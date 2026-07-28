# =============================================================================
# 04_multicohort_validation.R           [Fig 3 — 논문의 절반]
# KM · time-dependent ROC · DCA · IDH/grade/age/MGMT 보정 다변량 Cox · IDH 층화
# (선행연구는 CGGA693만 사용 → 여기선 693+325 모두)
# =============================================================================
source(here::here("R", "setup.R"))
suppressPackageStartupMessages({ library(survival); library(survminer); library(timeROC) })

cohorts <- load_all_cohorts(); stopifnot(length(cohorts) > 0)
roc_tab <- list(); mv_tab <- list(); st_tab <- list(); dca_all <- list()

for (nm in names(cohorts)) {
  co <- cohorts[[nm]]
  if (!GENE_OF_INTEREST %in% rownames(co$expr)) { message(nm, ": 유전자 없음"); next }
  v <- cohort_gene(co)
  km <- try(km_cohort(co), silent = TRUE)
  if (!inherits(km, "try-error"))
    ggplot2::ggsave(file.path(DIR_FIGURES, paste0("Fig3_KM_", nm, ".pdf")),
                    survminer::arrange_ggsurvplots(list(km), print = FALSE), width = 5.5, height = 6.5)
  roc_tab[[nm]] <- timedep_roc(co, v, label = GENE_OF_INTEREST)
  r <- try(plot_timedep_roc(co, v), silent = TRUE)
  if (!inherits(r, "try-error"))
    ggplot2::ggsave(file.path(DIR_FIGURES, paste0("Fig3_ROC_", nm, ".pdf")), r, width = 5, height = 4.5)
  mv <- multivariable_cox(co$expr, co$clin, GENE_OF_INTEREST)
  if (!is.null(mv)) mv_tab[[nm]] <- dplyr::mutate(mv$table[mv$table$term == "gene", ],
                                                  cohort = nm, label = nm)
  st_tab[[nm]]  <- dplyr::mutate(stratified_cox(co$expr, co$clin, GENE_OF_INTEREST), cohort = nm)
  dca_all[[nm]] <- tryCatch(dca_cohort(co), error = function(e) NULL)
}

roc_df <- dplyr::bind_rows(roc_tab); write_result(roc_df, "Fig3_timeROC_AUC.csv")
mv_df  <- dplyr::bind_rows(mv_tab);  write_result(mv_df,  "Fig3_multivariable_cox.csv")
write_result(dplyr::bind_rows(st_tab), "Fig3_stratified_cox_idh.csv")
if (nrow(mv_df)) ggplot2::ggsave(file.path(DIR_FIGURES, "Fig3_forest_multicohort.pdf"),
  forest_multicohort(mv_df), width = 6, height = 0.5 * nrow(mv_df) + 2)
if (length(dca_all)) {
  dca_df <- dplyr::bind_rows(dca_all); write_result(dca_df, "Fig3_DCA.csv")
  ggplot2::ggsave(file.path(DIR_FIGURES, "Fig3_DCA.pdf"), plot_dca(dca_df), width = 9, height = 4)
}
print(roc_df); print(mv_df); message("완료: 다중 코호트 검증")
