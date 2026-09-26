# =============================================================================
# 07c_ucell_module_axes.R
#   UCell 모듈 축 분석 — "SLC7A7 은 수송·대사 축이 아니라 면역(IFN) 축에 있다"
# -----------------------------------------------------------------------------
# 07 / 07b 를 대체합니다.
#
# 설계 원칙
#   1) 모듈은 전부 MSigDB(Hallmark / Reactome / GO:BP)에서 프로그램으로 생성.
#      공개 set 이 없는 2개(Pan-myeloid, Homeostatic microglia)만 curated 이며,
#      spatial 분석(10_spatial.R MODULES)과 글자 그대로 동일하게 유지한다.
#   2) SLC7A7 은 모든 모듈에서 제거 (순환논리 차단).
#   3) transport / arg_enzyme 은 ★음성대조★ — 이들이 null 이어야 논지가 선다.
#   4) 검정은 transport 모듈을 폐기시킨 것과 동일한 수준으로:
#        lmer(depth + dataset 랜덤효과) + within-dataset 랜덤효과 메타 + I²
#        + depth 음성대조 + 발현매칭 200유전자 empirical null
#   5) UCell maxRank 는 2값 민감도 (1500 / 3000).
#
# 입력 : data/processed/microglia.rds
#        results/Supp_null_matched_genes_microglia.csv   (11_ 산출물, empirical null 용)
#        results/Fig6_pseudobulk_DEG_limma.csv           (ISG 유전자 수준 확인용)
# 출력 : results/ucell/*.csv , figures/ucell/*.pdf
# =============================================================================

source(here::here("R", "setup.R"))
suppressPackageStartupMessages({
  library(Seurat); library(UCell); library(msigdbr); library(lme4); library(lmerTest)
  library(dplyr); library(tidyr); library(tibble); library(purrr); library(ggplot2)
})

RES_U <- file.path(DIR_RESULTS, "ucell"); dir.create(RES_U, recursive = TRUE, showWarnings = FALSE)
FIG_U <- file.path(DIR_FIGURES, "ucell"); dir.create(FIG_U, recursive = TRUE, showWarnings = FALSE)

MIN_DONOR_CELLS   <- 50     # 08 / 11 과 동일 (donor 67명)
MIN_DATASET_DONOR <- 3
MAXRANK_SET       <- c(1500, 3000)
N_NULL            <- 200

obj <- readRDS(file.path(DIR_DATA_PROC, "microglia.rds"))
ds_col <- find_batch_col(obj); stopifnot(!is.na(ds_col))
message(sprintf("[0] microglia cells = %d | median nFeature = %.0f",
                ncol(obj), median(obj$nFeature_RNA)))
message("    ▶ maxRank 는 median nFeature 보다 충분히 커야 함. 위 값 확인 후 MAXRANK_SET 조정.")

UNIVERSE <- rownames(obj)

# =============================================================================
# §1. 모듈 정의 — MSigDB 에서 프로그램으로 생성
# =============================================================================
msig_get <- function(gs_name, coll, subcoll = NULL) {
  m <- tryCatch(msigdbr(species = "Homo sapiens", collection = coll, subcollection = subcoll),
                error = function(e) msigdbr(species = "Homo sapiens", category = coll,
                                            subcategory = subcoll))
  g <- unique(m$gene_symbol[m$gs_name == gs_name])
  if (!length(g)) warning("빈 gene set: ", gs_name)
  g
}

H <- local({
  m <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"),
                error = function(e) msigdbr(species = "Homo sapiens", category = "H"))
  split(m$gene_symbol, m$gs_name)
})
stopifnot(length(H) == 50)

# ── A. 1차 검정 축 (Hallmark, 경험적 유래) ───────────────────────────────────
MOD <- list(
  IFN_ALPHA    = H$HALLMARK_INTERFERON_ALPHA_RESPONSE,
  IFN_GAMMA    = H$HALLMARK_INTERFERON_GAMMA_RESPONSE,
  NFKB         = H$HALLMARK_TNFA_SIGNALING_VIA_NFKB,
  INFLAMMATORY = H$HALLMARK_INFLAMMATORY_RESPONSE,

  # ── B. ★음성대조★ — 기존 curation 함수 그대로 (공개 GO 3종 조합) ──────────
  TRANSPORT    = build_transporter_geneset(universe = UNIVERSE),
  ARG_ENZYME   = build_arg_enzyme_geneset(universe = UNIVERSE),

  # ── C. 미세환경 / 대사 대조 (Hallmark) ─────────────────────────────────────
  HYPOXIA      = H$HALLMARK_HYPOXIA,
  GLYCOLYSIS   = H$HALLMARK_GLYCOLYSIS,
  OXPHOS       = H$HALLMARK_OXIDATIVE_PHOSPHORYLATION,
  ROS          = H$HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY,

  # ── D. 교차검증 (다른 DB · 작은 set) ───────────────────────────────────────
  IFN_REACTOME = msig_get("REACTOME_INTERFERON_ALPHA_BETA_SIGNALING", "C2", "CP:REACTOME"),
  PPP_REACTOME = msig_get("REACTOME_PENTOSE_PHOSPHATE_PATHWAY",        "C2", "CP:REACTOME"),
  POLYAMINE    = msig_get("GOBP_POLYAMINE_BIOSYNTHETIC_PROCESS",       "C5", "GO:BP"),

  # ── E. 조성 / 정체성 — 공개 set 부재. spatial(10_spatial.R)과 동일한 curated ──
  PAN_MYELOID  = c("PTPRC","AIF1","CSF1R","C1QA","C1QB","C1QC","CD68",
                   "TYROBP","FCER1G","ITGAM","CD14","LAPTM5","MPEG1"),
  HOMEOSTATIC  = c("P2RY12","TMEM119","CX3CR1","OLFML3","SALL1")
)
CURATED <- c("PAN_MYELOID", "HOMEOSTATIC")   # Methods 에 curated 로 명시할 것

# self-gene 및 파트너 제거 + universe 교집합
MOD <- lapply(MOD, function(g) setdiff(intersect(unique(g), UNIVERSE), GENE_OF_INTEREST))
MOD <- MOD[lengths(MOD) >= 5]
print(data.frame(module = names(MOD), n_gene = lengths(MOD),
                 source = ifelse(names(MOD) %in% CURATED, "curated", "MSigDB")), row.names = FALSE)
write.csv(stack(MOD) %>% setNames(c("gene","module")),
          file.path(RES_U, "UC_module_genes.csv"), row.names = FALSE)

# ── 겹침 행렬 (IFN vs NFKB 주장을 위해 반드시 보고) ──────────────────────────
jac <- outer(seq_along(MOD), seq_along(MOD), Vectorize(function(i, j)
  length(intersect(MOD[[i]], MOD[[j]])) / length(union(MOD[[i]], MOD[[j]]))))
dimnames(jac) <- list(names(MOD), names(MOD))
write.csv(round(jac, 3), file.path(RES_U, "UC_module_jaccard.csv"))
message(sprintf("\n[1] IFN_ALPHA ∩ IFN_GAMMA = %d | IFN_GAMMA ∩ NFKB = %d | IFN_ALPHA ∩ NFKB = %d",
                length(intersect(MOD$IFN_ALPHA, MOD$IFN_GAMMA)),
                length(intersect(MOD$IFN_GAMMA, MOD$NFKB)),
                length(intersect(MOD$IFN_ALPHA, MOD$NFKB))))

# ── disjoint 버전 : 2개 이상 모듈에 속한 유전자를 전부 제거 ──────────────────
tally <- table(unlist(MOD))
shared <- names(tally)[tally > 1]
MOD_DISJ <- lapply(MOD, setdiff, y = shared)
MOD_DISJ <- MOD_DISJ[lengths(MOD_DISJ) >= 5]
names(MOD_DISJ) <- paste0(names(MOD_DISJ), "_disj")
message(sprintf("    공유 유전자 %d개 제거 → disjoint 모듈 %d개 유지",
                length(shared), length(MOD_DISJ)))

# =============================================================================
# §2. UCell 점수 (maxRank 민감도 2값) → donor 집계
# =============================================================================
md <- obj@meta.data
donor_col <- if ("donor_id" %in% colnames(md)) "donor_id" else "donor"

score_all <- function(mr) {
  message("\n[2] UCell maxRank = ", mr)
  score_auc(obj, c(MOD, MOD_DISJ), max_rank = mr)
}
U <- setNames(lapply(MAXRANK_SET, score_all), paste0("mr", MAXRANK_SET))
saveRDS(U, file.path(DIR_DATA_PROC, "all_path_score.rds"))
# U <- readRDS(file.path(DIR_DATA_PROC, "all_path_score.rds"))

donor_meta <- md %>%
  mutate(donor = as.character(.data[[donor_col]]),
         dataset = as.character(.data[[ds_col]])) %>%
  group_by(donor) %>%
  summarise(slc     = mean(FetchData(obj, vars = GENE_OF_INTEREST)[cur_group_rows(), 1]),
            depth   = mean(.data[[DEPTH_COVARIATE]]),
            n_cell  = n(),
            dataset = names(which.max(table(dataset))), .groups = "drop") %>%
  filter(n_cell >= MIN_DONOR_CELLS) %>% as.data.frame()
# dataset 당 donor 1명이면 제외
donor_meta <- donor_meta[donor_meta$dataset %in%
                           names(which(table(donor_meta$dataset) >= MIN_DATASET_DONOR)), ]
rownames(donor_meta) <- donor_meta$donor
message(sprintf("    donors = %d | datasets = %d | rho(SLC7A7, depth) = %.3f",
                nrow(donor_meta), n_distinct(donor_meta$dataset),
                cor(donor_meta$slc, donor_meta$depth, method = "spearman")))

agg_donor <- function(sc) {
  d <- data.frame(donor = as.character(md[[donor_col]][match(rownames(sc), rownames(md))]), sc,
                  check.names = FALSE)
  d %>% group_by(donor) %>% summarise(across(everything(), mean), .groups = "drop") %>%
    filter(donor %in% donor_meta$donor) %>% as.data.frame()
}
D <- lapply(U, agg_donor)
for (k in names(D)) rownames(D[[k]]) <- D[[k]]$donor

obj_uc <- AddMetaData(obj, U$mr3000)
saveRDS(obj_uc, file.path(DIR_DATA_PROC, "microglia_scored.rds"))


# =============================================================================
# §3. 검정 엔진 — transport 를 폐기시킨 것과 동일한 절차
# =============================================================================
#' @param y  donor 수준 모듈 점수 (named)
#' @param x  donor 수준 예측변수 (named) — SLC7A7 또는 대조 유전자
test_module <- function(y, x, mm = donor_meta) {
  d <- data.frame(y = as.numeric(y[mm$donor]), x = scale(as.numeric(x[mm$donor]))[, 1],
                  depth = mm$depth, dataset = factor(mm$dataset))
  d <- d[complete.cases(d), ]
  if (nrow(d) < 20 || sd(d$y) == 0) return(NULL)

  ## (a) 혼합모형 (primary)
  fm <- try(lmerTest::lmer(y ~ x + depth + (1 | dataset), data = d), silent = TRUE)
  if (inherits(fm, "try-error")) return(NULL)
  cf <- summary(fm)$coefficients

  ## (b) within-dataset 랜덤효과 메타 + I²  ← transport 를 죽인 검정
  per <- d %>% group_by(dataset) %>% filter(n() >= 4) %>%
    group_modify(~{
      s <- summary(lm(y ~ x + depth, data = .x))$coefficients
      if (!"x" %in% rownames(s)) return(data.frame())
      data.frame(beta = s["x", 1], se = s["x", 2], n = nrow(.x))
    }) %>% ungroup() %>% filter(is.finite(beta), is.finite(se), se > 0)

  meta <- list(beta = NA, se = NA, p = NA, I2 = NA, k = 0)
  if (nrow(per) >= 2) {
    w <- 1 / per$se^2; bf <- sum(w * per$beta) / sum(w)
    Q <- sum(w * (per$beta - bf)^2); dfq <- nrow(per) - 1
    tau2 <- max(0, (Q - dfq) / (sum(w) - sum(w^2) / sum(w)))     # DerSimonian–Laird
    wr <- 1 / (per$se^2 + tau2); br <- sum(wr * per$beta) / sum(wr); ser <- sqrt(1 / sum(wr))
    meta <- list(beta = br, se = ser, p = 2 * pnorm(-abs(br / ser)),
                 I2 = max(0, (Q - dfq) / Q) * 100, k = nrow(per))
  }
  data.frame(
    beta_lmer = cf["x", 1], se_lmer = cf["x", 2], p_lmer = cf["x", 5],
    beta_meta = meta$beta, se_meta = meta$se, p_meta = meta$p,
    I2 = meta$I2, k_dataset = meta$k, n_donor = nrow(d),
    n_pos_dataset = if (nrow(per)) sum(per$beta > 0) else NA
  )
}

slc_vec <- setNames(donor_meta$slc,   donor_meta$donor)
dep_vec <- setNames(donor_meta$depth, donor_meta$donor)

run_panel <- function(Dk, xvec, label, mrtag) {
  mods <- setdiff(names(Dk), "donor")
  out <- map_dfr(mods, function(m) {
    r <- test_module(setNames(Dk[[m]], Dk$donor), xvec)
    if (is.null(r)) return(NULL)
    cbind(data.frame(module = m, predictor = label, maxRank = mrtag,
                     n_gene = length(c(MOD, MOD_DISJ)[[m]]),
                     source = ifelse(m %in% CURATED, "curated", "MSigDB")), r)
  })
  out$fdr_lmer <- p.adjust(out$p_lmer, "BH")
  out$fdr_meta <- p.adjust(out$p_meta, "BH")
  out
}

res <- map_dfr(names(D), function(k) {
  mr <- sub("^mr", "", k)
  rbind(run_panel(D[[k]], slc_vec, "SLC7A7", mr),
        run_panel(D[[k]], dep_vec, "depth_negctrl", mr))   # ★ depth 음성대조
})
write.csv(res, file.path(RES_U, "UC_module_regression.csv"), row.names = FALSE)

main <- res %>% filter(predictor == "SLC7A7", maxRank == MAXRANK_SET[1],
                       !grepl("_disj$", module)) %>% arrange(p_meta)
cat("\n===== 1차 결과 (maxRank ", MAXRANK_SET[1], ", SLC7A7) =====\n", sep = "")
print(main %>% dplyr::select(module, n_gene, beta_lmer, p_lmer, fdr_lmer,
                      beta_meta, p_meta, fdr_meta, I2, k_dataset, n_pos_dataset),
      row.names = FALSE, digits = 3)
message("\n  ▶ 판정: transport / arg_enzyme 이 ns 이고 IFN 계열이 fdr_meta < 0.05 & I² 낮음 → 논지 성립")
message("  ▶ I² > 50 이면 '유의'해도 dataset 간 이질적이므로 메타 서술로 낮출 것")

# =============================================================================
# §4. 발현매칭 200유전자 empirical null  (11_ 과 동일 프레임)
# =============================================================================
mfile <- file.path(DIR_RESULTS, "Supp_null_matched_genes_microglia.csv")
if (file.exists(mfile)) {
  matched <- utils::read.csv(mfile)$gene
  matched <- intersect(matched, UNIVERSE)[seq_len(min(N_NULL, length(matched)))]
  message(sprintf("\n[4] empirical null: 매칭 대조 유전자 %d개", length(matched)))

  dat <- SeuratObject::LayerData(obj, assay = "RNA", layer = "data")
  cells_keep <- md[[donor_col]] %in% donor_meta$donor
  dfac <- factor(as.character(md[[donor_col]][cells_keep]), levels = donor_meta$donor)
  ind  <- Matrix::sparse.model.matrix(~ 0 + dfac); colnames(ind) <- levels(dfac)
  ncell <- Matrix::colSums(ind)
  gm <- sweep(as.matrix(dat[matched, cells_keep, drop = FALSE] %*% ind), 2, ncell, "/")

  Dk <- D[[1]]; mods <- setdiff(names(Dk), "donor"); mods <- mods[!grepl("_disj$", mods)]
  pbar <- txtProgressBar(min = 0, max = length(matched), style = 3)
  null_list <- vector("list", length(matched))
  for (i in seq_along(matched)) {
    xv <- setNames(gm[matched[i], donor_meta$donor], donor_meta$donor)
    null_list[[i]] <- map_dfr(mods, function(m) {
      r <- test_module(setNames(Dk[[m]], Dk$donor), xv)
      if (is.null(r)) return(NULL)
      data.frame(gene = matched[i], module = m, beta = r$beta_lmer, beta_meta = r$beta_meta)
    })
    setTxtProgressBar(pbar, i)
  }
  close(pbar)
  nulldf <- bind_rows(null_list)
  write.csv(nulldf, file.path(RES_U, "UC_null_distribution.csv"), row.names = FALSE)

  obs <- res %>% dplyr::filter(predictor == "SLC7A7", maxRank == MAXRANK_SET[1]) %>%
    dplyr::select(module, beta_lmer, beta_meta)
  emp <- nulldf %>% group_by(module) %>%
    summarise(n_null = n(), null_mean = mean(beta, na.rm = TRUE),
              null_sd = sd(beta, na.rm = TRUE), .groups = "drop") %>%
    inner_join(obs, by = "module") %>%
    rowwise() %>%
    mutate(z_vs_null = (beta_lmer - null_mean) / null_sd,
           p_emp_two = (1 + sum(abs(nulldf$beta[nulldf$module == module]) >= abs(beta_lmer),
                                na.rm = TRUE)) / (1 + n_null)) %>%
    ungroup() %>% mutate(p_emp_BH = p.adjust(p_emp_two, "BH")) %>% arrange(p_emp_two)
  print(as.data.frame(emp), row.names = FALSE, digits = 3)
  write.csv(emp, file.path(RES_U, "UC_empirical_p.csv"), row.names = FALSE)
} else {
  message("\n[4] 건너뜀 — ", basename(mfile), " 없음 (11_depth_sensitivity.R 먼저 실행)")
}

# =============================================================================
# §5. ISG 유전자 수준 확인 — 기존 DEG 표에서 바로 (spatial module control 과 동형)
# =============================================================================
dfile <- file.path(DIR_RESULTS, "Fig6_pseudobulk_DEG_limma.csv")
if (file.exists(dfile)) {
  deg <- utils::read.csv(dfile)
  ISG <- c("IFIT1","IFIT3","ISG15","OAS1","OAS2","OAS3","IFI6",
           "USP18","HERC6","CMPK2","RSAD2","STAT1","EIF2AK2","XAF1")   # spatial 과 동일
  gl <- deg[deg$gene %in% ISG, ]
  gl <- gl[order(gl$P.Value), ]
  print(gl, row.names = FALSE, digits = 3)
  write.csv(gl, file.path(RES_U, "UC_ISG_genelevel.csv"), row.names = FALSE)

  bg <- deg$t[is.finite(deg$t)]
  cat(sprintf("\n[5] ISG %d/%d detected | mean t = %.2f (배경 %.2f) | Wilcoxon p = %.3g\n",
              nrow(gl), length(ISG), mean(gl$t), mean(bg),
              wilcox.test(gl$t, bg)$p.value))
  if (requireNamespace("limma", quietly = TRUE)) {
    st <- setNames(deg$t, deg$gene); st <- st[is.finite(st)]
    print(limma::cameraPR(st, limma::ids2indices(list(ISG = ISG, PAN_MYELOID = MOD$PAN_MYELOID,
                                                      NFKB_EFF = MOD$NFKB), names(st))))
  }
  message("  ▶ 이 표가 세포 내재적 IFN 주장의 1차 근거. gene set 검정보다 상관 페널티가 작다.")
} else {
  message("\n[5] 건너뜀 — Fig6_pseudobulk_DEG_limma.csv 없음")
}

# =============================================================================
# §6. 그림
# =============================================================================
plt <- res %>% filter(maxRank == MAXRANK_SET[1], !grepl("_disj$", module)) %>%
  mutate(module = factor(module, levels = rev(main$module)),
         sig = fdr_meta < 0.05)
p1 <- ggplot(plt, aes(beta_meta, module, colour = predictor)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_errorbarh(aes(xmin = beta_meta - 1.96 * se_meta, xmax = beta_meta + 1.96 * se_meta),
                 height = .25, position = position_dodge(width = .6)) +
  geom_point(aes(shape = sig), size = 2.6, position = position_dodge(width = .6)) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), name = "FDR < 0.05") +
  scale_colour_manual(values = c(SLC7A7 = PALETTE_TWO[1], depth_negctrl = "grey55")) +
  labs(x = "random-effects meta β (within-dataset)", y = NULL,
       title = "UCell module ~ SLC7A7 (depth 음성대조 병기)") +
  theme_paper(base_size = 10)
ggsave(file.path(FIG_U, "UC_module_forest.pdf"), p1, width = 7.5, height = 5.5)

p2 <- plt %>% filter(predictor == "SLC7A7") %>%
  ggplot(aes(I2, module, fill = sig)) +
  geom_col(width = .65) + geom_vline(xintercept = 50, linetype = 2, colour = "firebrick") +
  scale_fill_manual(values = c(`FALSE` = "grey75", `TRUE` = PALETTE_TWO[1]), guide = "none") +
  labs(x = "I² (%) — dataset 간 이질성", y = NULL,
       title = "50%% 초과면 '일관된 효과'로 서술 불가") + theme_paper(base_size = 10)
ggsave(file.path(FIG_U, "UC_module_I2.pdf"), p2, width = 6.5, height = 5.5)

message("\ndone.")
message(" 본문 표   : results/ucell/UC_module_regression.csv")
message(" 본문 그림 : figures/ucell/UC_module_forest.pdf (+ UC_module_I2.pdf)")
message(" 핵심 근거 : results/ucell/UC_ISG_genelevel.csv , UC_empirical_p.csv")
message("\n※ 보고 원칙")
message("  · TRANSPORT / ARG_ENZYME 의 null 결과를 반드시 같은 표·같은 그림에 넣을 것 (음성대조)")
message("  · I² 를 모든 모듈에 대해 보고. IFN 만 보고하면 이중잣대.")
message("  · PAN_MYELOID / HOMEOSTATIC 은 curated 임을 Methods 에 명시 (공개 set 부재)")

