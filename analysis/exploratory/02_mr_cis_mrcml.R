# =============================================================================
# exploratory/02_mr_cis_mrcml.R   *** 논문 미사용 (탐색용 보존) ***
# SLC family cis-MR (TwoSampleMR + 민감도) + MRcML
#   exposure = eQTLGen cis-eQTL(OpenGWAS), outcome = 뇌종양/GBM GWAS
# 입력: data/processed/SLC_coxtest.csv    출력: results/MR_*.csv
# 원본: "MR통합 찐 최종 코드.R" (MR.R / MR자동화.R / MRcML자동화.R의 최종 통합본)
# -----------------------------------------------------------------------------
# ※ 결과의 유의성/일관성이 낮아 논문에는 포함하지 않았습니다.
# ※ OPENGWAS_JWT 는 .Renviron 에서 주입 (코드에 토큰 하드코딩 금지).
# =============================================================================
source(here::here("R", "setup.R"))
suppressPackageStartupMessages({
  library(TwoSampleMR); library(MRcML); library(MRPRESSO)
  library(biomaRt); library(jsonlite); library(org.Hs.eg.db)
})

Sys.setenv(OPENGWAS_JWT = get_opengwas_token())   # .Renviron 에서 읽음

# exposure 유전자 목록 (p<0.05 3군 모두 만족하는 SLC)
df_slc <- utils::read.csv(PATH_SLC_COXTEST)
df_filtered <- dplyr::filter(df_slc, overall_p <= 0.05, mut_p <= 0.05, wt_p <= 0.05)
ensembl_list <- stats::na.omit(df_filtered$ENSEMBL)

# outcome 선택 (config의 MR_OUTCOMES)
outcome  <- MR_OUTCOMES$gbm
mart     <- connect_mart()

# TwoSampleMR
mr_res_list <- run_mr_batch(ensembl_list, outcome$id, outcome$n, mart = mart)
mr_results  <- do.call(rbind, lapply(mr_res_list, function(x) x$summary))
het         <- do.call(rbind, lapply(mr_res_list, function(x) x$heterogeneity))
steiger     <- do.call(rbind, lapply(mr_res_list, function(x) x$steiger))
plt         <- do.call(rbind, lapply(mr_res_list, function(x) x$pleiotropy))
presso      <- do.call(rbind, lapply(mr_res_list, function(x) x$presso_summary))
loo         <- do.call(rbind, lapply(mr_res_list, function(x) x$leaveoneout))

write_result(mr_results, "MR_results_SLC_GBM.csv")
write_result(het,        "MR_heterogeneity.csv")
write_result(steiger,    "MR_steiger.csv")
write_result(plt,        "MR_pleiotropy.csv")
if (!is.null(presso)) write_result(presso, "MR_presso.csv")
write_result(loo,        "MR_leaveoneout.csv")

# MRcML
cml_list    <- run_cml_batch(ensembl_list, outcome$id, outcome$n, mart = mart)
cml_summary <- do.call(rbind, lapply(cml_list, function(x) x$summary))
write_result(cml_summary, "MRcML_summary.csv")
saveRDS(list(mr = mr_res_list, cml = cml_list), file.path(DIR_RESULTS, "MR_full.rds"))
message("완료: MR + MRcML (exploratory)")
