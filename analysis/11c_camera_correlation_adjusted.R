# =============================================================================
# 11c_camera_correlation_adjusted.R
#   — fgsea 의 NES 팽창(유전자간 상관) 보정 + 데이터셋 일관성 + 조성 보정
# -----------------------------------------------------------------------------
# 진단 결과 요약 (11b):
#   D2  depth 단독 예측변수 → OXPHOS NES = -1.569 (관측 +2.392 와 반대 방향)
#       ⇒ depth 인공산물 가설 기각
#   D3  permutation null (dataset 내 셔플) → perm_sd = 1.18 ~ 2.27
#       ⇒ 이론값 ~1 의 2배. **fgsea 의 p 값 자체가 팽창**돼 있다.
#       원인: fgsea 는 '유전자'를 셔플하므로 유전자간 독립을 가정하는데,
#             실제 전사체는 pathway 단위로 강하게 공조절된다.
#             (Efron & Tibshirani 2007; Wu & Smyth 2012 CAMERA)
#   D4  SVA 는 이 상황의 해법이 아님 — 확산형(diffuse) 신호를 SV 가 흡수.
#       부호가 0 이 아니라 반대(-1.32)로 뒤집힌 것이 과보정의 징후.
#
# 따라서 이 스크립트는 SVA 대신
#   §1 CAMERA / FRY  — 유전자간 상관을 보정한 gene set 검정 (표준 해법)
#   §2 SV 진단       — SV 가 SLC7A7 자체를 흡수했는지 3줄로 확인
#   §3 LODO + 데이터셋별 NES — "한 study 가 끌고 있는가"
#   §4 조성 보정     — 아종 비율을 공변량으로 (Fig6 pro-infl I FDR=0.045 대응)
# 를 수행한다.  11_depth_sensitivity.R 실행 후 같은 세션에서 이어서 실행.
# =============================================================================

suppressPackageStartupMessages({ library(limma); library(edgeR); library(fgsea); library(dplyr) })

x0     <- as.numeric(dm[GENE, colnames(dge0$counts)])
design <- model.matrix(~ scale(x0) + donor_meta$depth + donor_meta$dataset)
colnames(design) <- make.names(colnames(design))
COEF <- colnames(design)[2]
v    <- voom(dge0, design, plot = FALSE)
idx  <- ids2indices(HALLMARK, rownames(v))

# =============================================================================
# §1. CAMERA (경쟁형, 유전자간 상관 보정) + FRY (회전검정)
# =============================================================================
message("[1] CAMERA — inter-gene correlation 을 데이터에서 추정")
cam <- camera(v, idx, design, contrast = COEF, inter.gene.cor = NA, sort = FALSE)
cam$pathway <- rownames(cam)
print(cam[TARGET_PATHWAYS, c("NGenes","Correlation","Direction","PValue","FDR")], digits = 3)
message("   ▶ Correlation 이 0.05 이상이면 fgsea 의 p 는 신뢰 불가 (팽창 확정)")

message("\n[1b] FRY — 표본 회전(rotation) 기반. 표본 상관구조를 보존")
fr <- fry(v, idx, design, contrast = COEF, sort = FALSE)
fr$pathway <- rownames(fr)
print(fr[TARGET_PATHWAYS, c("NGenes","Direction","PValue","FDR")], digits = 3)

cmp <- data.frame(
  pathway   = TARGET_PATHWAYS,
  fgsea_NES = as.numeric(obs_nes[TARGET_PATHWAYS]),
  camera_dir = cam[TARGET_PATHWAYS, "Direction"],
  camera_p   = cam[TARGET_PATHWAYS, "PValue"],
  camera_FDR = cam[TARGET_PATHWAYS, "FDR"],
  intergene_cor = cam[TARGET_PATHWAYS, "Correlation"],
  fry_p      = fr[TARGET_PATHWAYS, "PValue"],
  fry_FDR    = fr[TARGET_PATHWAYS, "FDR"]
)
write.csv(cmp, file.path(DIR_RESULTS, paste0("Fig6_GSEA_camera", SFX, ".csv")), row.names = FALSE)
print(cmp, digits = 3)
message("   ▶ 이 표의 camera_FDR / fry_FDR 을 논문에 보고하세요. fgsea padj 는 보고 금지.")

# 전체 50개 hallmark 도 저장 (본문 Figure 용)
cam_all <- camera(v, idx, design, contrast = COEF, inter.gene.cor = NA)
write.csv(data.frame(pathway = rownames(cam_all), cam_all),
          file.path(DIR_RESULTS, paste0("Fig6_GSEA_camera_all", SFX, ".csv")), row.names = FALSE)

# =============================================================================
# §2. SVA 가 생물학을 먹었는지 3줄 진단  (11b D4 해석용)
# =============================================================================
if (exists("svobj")) {
  message("\n[2] SV vs SLC7A7 상관 — |r| > 0.4 면 SV 가 신호를 흡수한 것 (D4 무효)")
  sv_cor <- apply(svobj$sv, 2, function(s) cor(s, scale(x0)[,1]))
  print(round(sv_cor, 3))
  message("    SV ~ dataset R²: ",
          paste(round(apply(svobj$sv, 2, function(s)
            summary(lm(s ~ donor_meta$dataset))$r.squared), 2), collapse = ", "))
  message("    SV ~ depth   r : ",
          paste(round(apply(svobj$sv, 2, function(s) cor(s, donor_meta$depth)), 2), collapse = ", "))
}

# =============================================================================
# §3. 데이터셋 일관성 — LODO + 데이터셋별 방향
# -----------------------------------------------------------------------------
#  Fig5 수송모듈이 9개 중 7개 음수였던 전례가 있으므로 반드시 확인.
# =============================================================================
message("\n[3] leave-one-dataset-out")
lodo <- lapply(levels(donor_meta$dataset), function(ds) {
  k <- donor_meta$dataset != ds
  if (sum(k) < 25) return(NULL)
  dgeK <- dge0[, k]; dgeK$samples$lib.size <- colSums(dgeK$counts)
  mdK  <- droplevels(donor_meta[k, ])
  dsn  <- model.matrix(~ scale(x0[k]) + mdK$depth + mdK$dataset)
  colnames(dsn) <- make.names(colnames(dsn))
  vK   <- voom(dgeK, dsn, plot = FALSE)
  cK   <- camera(vK, ids2indices(HALLMARK, rownames(vK)), dsn,
                 contrast = colnames(dsn)[2], inter.gene.cor = NA, sort = FALSE)
  data.frame(dropped = ds, n = sum(k),
             OXPHOS_dir = cK["HALLMARK_OXIDATIVE_PHOSPHORYLATION","Direction"],
             OXPHOS_p   = cK["HALLMARK_OXIDATIVE_PHOSPHORYLATION","PValue"],
             OXPHOS_FDR = cK["HALLMARK_OXIDATIVE_PHOSPHORYLATION","FDR"])
}) %>% bind_rows()
print(lodo, digits = 3)
write.csv(lodo, file.path(DIR_RESULTS, paste0("Supp_LODO_camera", SFX, ".csv")), row.names = FALSE)
message("   ▶ 어느 하나를 빼도 Up/유의 유지되어야 함. 특정 dataset 제거 시 무너지면 그 study 효과.")

# 데이터셋별 OXPHOS leading-edge 점수 방향 (donor 수 >= 8 인 dataset 만)
message("\n[3b] dataset 별 방향 (donor >= 8)")
oxg <- intersect(HALLMARK$HALLMARK_OXIDATIVE_PHOSPHORYLATION, rownames(dge0))
lcpm <- cpm(dge0, log = TRUE, prior.count = 1)
ox_score <- rowMeans(scale(t(lcpm[oxg, ])))          # donor 별 OXPHOS 평균 z
per_ds <- lapply(levels(donor_meta$dataset), function(ds) {
  k <- donor_meta$dataset == ds
  if (sum(k) < 8) return(NULL)
  m <- lm(ox_score[k] ~ scale(x0[k]) + donor_meta$depth[k])
  data.frame(dataset = ds, n = sum(k),
             beta = coef(m)[2], se = summary(m)$coef[2,2], p = summary(m)$coef[2,4])
}) %>% bind_rows()
print(per_ds, digits = 3)
write.csv(per_ds, file.path(DIR_RESULTS, paste0("Supp_OXPHOS_per_dataset", SFX, ".csv")), row.names = FALSE)

# =============================================================================
# §4. 조성 보정 — 아종 비율을 공변량으로
# -----------------------------------------------------------------------------
#  Fig6 에서 TAM-MG pro-infl I 비율이 SLC7A7 과 연관(β=+0.26, FDR=0.045).
#  "SLC7A7 높은 donor에 pro-infl I 이 많고 그 아종이 OXPHOS 가 높아서"
#  나온 신호인지, 세포 내부 상태인지 가른다. ← 생물학적으로 해석 가능한 검정
# =============================================================================
prop_file <- file.path(DIR_RESULTS, "Fig6_donor_subtype_proportions.csv")
if (file.exists(prop_file)) {
  message("\n[4] 아종 비율 보정")
  pr <- read.csv(prop_file)
  # donor × subtype wide 로 정리 (열 이름은 실제 파일에 맞게 조정)
  prw <- pr %>% tidyr::pivot_wider(id_cols = donor, names_from = subtype,
                                   values_from = prop, values_fill = 0) %>%
    as.data.frame()
  rownames(prw) <- prw$donor
  prw <- prw[colnames(dge0$counts), setdiff(colnames(prw), "donor"), drop = FALSE]
  prw <- prw[, colSums(is.na(prw)) == 0, drop = FALSE]
  # 마지막 열은 합=1 제약으로 중복 → 제거
  P <- as.matrix(prw)[, -ncol(prw), drop = FALSE]
  
  dsn2 <- model.matrix(~ scale(x0) + donor_meta$depth + donor_meta$dataset + P)
  colnames(dsn2) <- make.names(colnames(dsn2))
  v2   <- voom(dge0, dsn2, plot = FALSE)
  cam2 <- camera(v2, idx, dsn2, contrast = colnames(dsn2)[2],
                 inter.gene.cor = NA, sort = FALSE)
  out4 <- data.frame(pathway = TARGET_PATHWAYS,
                     dir_compadj = cam2[TARGET_PATHWAYS, "Direction"],
                     p_compadj   = cam2[TARGET_PATHWAYS, "PValue"],
                     FDR_compadj = cam2[TARGET_PATHWAYS, "FDR"])
  print(out4, digits = 3)
  write.csv(out4, file.path(DIR_RESULTS, paste0("Supp_GSEA_composition_adj", SFX, ".csv")),
            row.names = FALSE)
  message("   ▶ 조성 보정 후에도 OXPHOS 유의 = 세포 내부 상태 (강한 주장)")
  message("   ▶ 조성 보정 후 소실 = 아종 구성 효과 (그것도 정당한 결론, 서술만 바꾸면 됨)")
} else {
  message("\n[4] Fig6_donor_subtype_proportions.csv 없음 — 08 스크립트 먼저 실행")
}

message("\n=== 판정 ===")
message(" camera/fry FDR < 0.05 & LODO 안정 & 데이터셋별 방향 일치 → 보고 가능")
message(" camera 는 유의한데 dataset 별 방향이 갈리면 → 메타분석으로 서술 (Fig5 방식)")
message(" camera 도 ns → 단일세포 GSEA 주장 철회, 기술적(descriptive) 층으로 축소")
