# =============================================================================
# 11b_null_diagnostics.R  —  empirical null 이 실패한 원인 진단
# -----------------------------------------------------------------------------
# 배경: microglia 에서 OXPHOS 관측 NES = +2.163 인데 매칭 200유전자 null 의
#       SD 가 2.344 (기대 ~1) 로 과대산포 → p_emp = 0.611 (통과 실패).
#       SLC3A2 편상관도 null 평균 t = +3.11 로, "아무 유전자나 SLC7A7 과
#       양의 연관을 보이는" 전역 donor 수준 공분산 구조가 존재함을 시사.
#
# 이 스크립트는 "OXPHOS 가 진짜 아닌가" vs "모델이 불안정한가" 를 가른다.
# 11_depth_sensitivity.R 을 먼저 끝까지 실행한 뒤 같은 세션에서 이어서 실행.
#   (dm, dge0, donor_meta, HALLMARK, run_donor_continuous, nes_of 를 재사용)
# =============================================================================

suppressPackageStartupMessages({library(sva); library(limma); library(edgeR); library(fgsea) })

OUT <- list()

# =============================================================================
# 진단 1. 모델 안정성 — 자유도와 표본
# =============================================================================
n_donor <- ncol(dge0$counts)
n_ds    <- nlevels(donor_meta$dataset)
df_res  <- n_donor - (2 + n_ds - 1)          # intercept + x + depth + (ds-1)
message(sprintf("[D1] donors=%d  datasets=%d  residual df=%d  genes=%d",
                n_donor, n_ds, df_res, nrow(dge0)))
message("     ▶ residual df < 40 이면 voom 가중치가 불안정 → null 과대산포의 주원인")
message("     ▶ genes=16272 는 myeloid DEG(4537)보다 3.6배 많음. filterByExpr 설정 확인 요망")

# =============================================================================
# 진단 2. depth 자체를 예측변수로 → "depth 가 OXPHOS 를 만드는가?"
# -----------------------------------------------------------------------------
#  이게 원래 묻고 싶었던 질문입니다. depth 단독으로 OXPHOS NES 가 크게 양수면
#  → depth 인공산물 확정. 0 근처면 → depth 는 범인이 아님.
# =============================================================================
run_with_x <- function(x, drop_depth = FALSE) {
  x <- scale(as.numeric(x))[, 1]
  design <- if (drop_depth) model.matrix(~ x + donor_meta$dataset)
            else            model.matrix(~ x + donor_meta$depth + donor_meta$dataset)
  colnames(design) <- make.names(colnames(design))
  v   <- voom(dge0, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design))
  tt  <- topTable(fit, coef = "x", number = Inf, sort.by = "none")
  st  <- setNames(tt$t, rownames(tt))
  st[setdiff(names(st), GENE)]
}

set.seed(GSEA_SEED)
nes_depth <- nes_of(run_with_x(donor_meta$depth, drop_depth = TRUE), seed = GSEA_SEED)
OUT$depth_as_predictor <- nes_depth
message("\n[D2] depth 를 예측변수로 넣었을 때의 NES")
print(round(nes_depth, 3))
message("     ▶ OXPHOS 가 +1.5 이상이면 depth 인공산물, 0 근처면 depth 무죄")

# 참고: donor당 세포수(n_cell)도 같이 확인
OUT$ncell_as_predictor <- nes_of(run_with_x(donor_meta$n_cell), seed = GSEA_SEED)
message("[D2b] n_cell 을 예측변수로:"); print(round(OUT$ncell_as_predictor, 3))

# =============================================================================
# 진단 3. 진짜 permutation null — "우연인가?"
# -----------------------------------------------------------------------------
#  SLC7A7 donor 값을 dataset 안에서 섞는다(배치구조 보존).
#  이 null 의 SD 는 이론상 ~1 이어야 한다. 실측이 1 근처면 모델은 정상이고
#  매칭유전자 null 의 SD 2.3 은 "실제 유전자들이 진짜 생물학을 갖기 때문"이다.
#  이 null 의 SD 도 2 를 넘으면 → 모델 자체가 불안정 (진단 4로).
# =============================================================================
N_PERM <- 200
x0 <- as.numeric(dm[GENE, colnames(dge0$counts)])
perm_mat <- matrix(NA_real_, N_PERM, length(TARGET_PATHWAYS),
                   dimnames = list(NULL, TARGET_PATHWAYS))
set.seed(GSEA_SEED)
pbar <- txtProgressBar(min = 0, max = N_PERM, style = 3)
for (i in seq_len(N_PERM)) {
  xp <- ave(x0, donor_meta$dataset, FUN = sample)   # dataset 내 섞기
  perm_mat[i, ] <- nes_of(run_with_x(xp), seed = GSEA_SEED + 10000 + i)
  setTxtProgressBar(pbar, i)
}
close(pbar)

perm_summary <- data.frame(
  pathway   = TARGET_PATHWAYS,
  observed  = as.numeric(obs_nes[TARGET_PATHWAYS]),
  perm_mean = apply(perm_mat, 2, mean, na.rm = TRUE),
  perm_sd   = apply(perm_mat, 2, sd,   na.rm = TRUE),
  p_perm_two = sapply(TARGET_PATHWAYS, function(p) {
    nl <- perm_mat[, p]; nl <- nl[is.finite(nl)]
    (1 + sum(abs(nl) >= abs(obs_nes[[p]]))) / (1 + length(nl))
  })
)
perm_summary$p_perm_BH <- p.adjust(perm_summary$p_perm_two, "BH")
message("\n[D3] permutation null (dataset 내 셔플)")
print(perm_summary, digits = 3)
message("     ▶ perm_sd 가 ~1 이면 모델 정상 / >2 면 모델 불안정")
write.csv(perm_summary,
          file.path(DIR_RESULTS, paste0("Supp_null_permutation", SFX, ".csv")),
          row.names = FALSE)

# =============================================================================
# 진단 4. 숨은 요인(surrogate variable) 보정 후 재실행
# -----------------------------------------------------------------------------
#  SLC3A2 편상관 null 평균 t = +3.11 은 "모든 유전자가 SLC7A7 과 양의 연관"
#  = 모델에 없는 전역 donor 요인(세포 조성·해리 스트레스·ambient RNA 등)이
#  존재한다는 뜻. svaseq 로 잡아 공변량에 추가한다.
# =============================================================================
mod0 <- model.matrix(~ donor_meta$depth + donor_meta$dataset)
mod1 <- model.matrix(~ scale(x0) + donor_meta$depth + donor_meta$dataset)
cpm_mat <- cpm(dge0, log = TRUE, prior.count = 1)
n_sv <- tryCatch(num.sv(cpm_mat, mod1, method = "be"), error = function(e) 2)
n_sv <- max(1, min(n_sv, 5))
message(sprintf("\n[D4] estimated surrogate variables: %d", n_sv))
svobj <- sva(cpm_mat, mod1, mod0, n.sv = n_sv)

run_with_sv <- function(x) {
  x <- scale(as.numeric(x))[, 1]
  design <- model.matrix(~ x + donor_meta$depth + donor_meta$dataset + svobj$sv)
  colnames(design) <- make.names(colnames(design))
  v   <- voom(dge0, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design))
  tt  <- topTable(fit, coef = "x", number = Inf, sort.by = "none")
  st  <- setNames(tt$t, rownames(tt)); st[setdiff(names(st), GENE)]
}

obs_sv <- nes_of(run_with_sv(x0), seed = GSEA_SEED)
message("[D4] SV 보정 후 관측 NES"); print(round(obs_sv, 3))

# SV 보정 상태에서 매칭유전자 null 재실행 (100회로 축약, 시간 절약)
N2 <- 100
mat2 <- matrix(NA_real_, N2, length(TARGET_PATHWAYS),
               dimnames = list(matched$gene[1:N2], TARGET_PATHWAYS))
pbar <- txtProgressBar(min = 0, max = N2, style = 3)
for (i in seq_len(N2)) {
  mat2[i, ] <- nes_of(run_with_sv(dm[matched$gene[i], colnames(dge0$counts)]),
                      seed = GSEA_SEED + 20000 + i)
  setTxtProgressBar(pbar, i)
}
close(pbar)

sv_summary <- data.frame(
  pathway  = TARGET_PATHWAYS,
  observed_SVadj = as.numeric(obs_sv[TARGET_PATHWAYS]),
  null_mean = apply(mat2, 2, mean, na.rm = TRUE),
  null_sd   = apply(mat2, 2, sd,   na.rm = TRUE),
  p_emp_two = sapply(TARGET_PATHWAYS, function(p) {
    nl <- mat2[, p]; nl <- nl[is.finite(nl)]
    (1 + sum(abs(nl) >= abs(obs_sv[[p]]))) / (1 + length(nl))
  })
)
sv_summary$p_emp_BH <- p.adjust(sv_summary$p_emp_two, "BH")
message("\n[D4] SV 보정 후 empirical null")
print(sv_summary, digits = 3)
write.csv(sv_summary,
          file.path(DIR_RESULTS, paste0("Supp_null_SVadjusted", SFX, ".csv")),
          row.names = FALSE)

# =============================================================================
# 진단 5. 공조절 유전자 배제 기준 민감도 (|rho| < 0.5 → 0.25)
# -----------------------------------------------------------------------------
#  현재 매칭 200개 중 |rho| 최대 0.498. EIF4B(0.481) 처럼 SLC7A7 과 꽤 상관된
#  유전자가 포함돼 null 을 부풀렸을 수 있다. 정직하게 두 기준 모두 보고할 것.
# =============================================================================
message("\n[D5] |rho| < 0.25 로 재선정 시 후보 수: ",
        sum(abs(gstat$rho_slc7a7) < 0.25, na.rm = TRUE), " genes")
message("     ▶ 200개 이상이면 EXCLUDE_COR_ABS <- 0.25 로 11_ 재실행하여 병기")

saveRDS(OUT, file.path(DIR_RESULTS, paste0("null_diagnostics", SFX, ".rds")))
message("\n=== 판정 가이드 ===")
message(" D2 OXPHOS |NES| > 1.5      → depth 인공산물. 결과 폐기")
message(" D3 perm_sd ~1 & p<0.05     → 모델 정상, 신호 실재. D1 매칭null은 '평범한 유전자보다 강한가' 라는 더 엄격한 질문")
message(" D3 perm_sd > 2             → 모델 불안정. D4 SV 보정 필수")
message(" D4 보정 후 p_emp < 0.05    → 통과. 논문에 SV 보정 버전으로 보고")
message(" 모두 실패                  → OXPHOS 주장 철회, bulk 중심 논문으로 축소")
