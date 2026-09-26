# =============================================================================
# 03_lasso_cox_signature.R              [Fig 2]
# 수송체 gene set → LASSO Cox(10-fold CV) → 유전자 선정 → 단/다변량 Cox forest
# ★ SLC7A7 탈락 시 candidate-gene 설계로 보고(양쪽 결과 모두 산출)
# =============================================================================
source(here::here("R", "setup.R"))
suppressPackageStartupMessages({ library(glmnet); library(survival); library(survminer) })

co    <- load_cohort(TRAIN_COHORT)
stopifnot(ncol(co$expr) > 0, nrow(co$clin) > 0)     # ← 코호트 로드 검증
genes <- utils::read.csv(file.path(DIR_RESULTS, "transporter_geneset.csv"))$gene

las <- run_lasso_cox(co$expr, co$clin, genes)
write_result(las$selected, "LASSO_selected_genes.csv")
grDevices::pdf(file.path(DIR_FIGURES, "Fig2AB_lasso_path_cv.pdf"), width = 10, height = 5)
plot_lasso(las); grDevices::dev.off()
message(GENE_OF_INTEREST, if (GENE_OF_INTEREST %in% las$selected$gene) " 는 LASSO 선택됨."
        else " 는 LASSO 탈락 → candidate-gene 설계로 보고할 것.")

# plot_lasso(las)
uni <- univariate_cox(co$expr, co$clin, genes); write_result(uni, "univariate_cox_transporters.csv")
mv  <- multivariable_cox(co$expr, co$clin, GENE_OF_INTEREST)
if (!is.null(mv)) { print(mv$table); write_result(mv$table, "multivariable_cox_train.csv") }
write_result(stratified_cox(co$expr, co$clin, GENE_OF_INTEREST, by = "idh"),
             "stratified_cox_idh_train.csv")

top <- utils::head(uni, 15)
p <- ggplot2::ggplot(top, ggplot2::aes(HR, stats::reorder(gene, HR))) +
  ggplot2::geom_vline(xintercept = 1, color = "grey60") +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = lower, xmax = upper), height = .2) +
  ggplot2::geom_point(ggplot2::aes(color = FDR < 0.05), size = 2.6) +
  ggplot2::scale_color_manual(values = c(`TRUE` = "#D64B4B", `FALSE` = "grey60")) +
  ggplot2::scale_x_continuous(trans = "log2") +
  ggplot2::labs(x = "HR per SD", y = NULL, title = paste0(TRAIN_COHORT, " — transporter Cox")) +
  theme_paper()
p
ggplot2::ggsave(file.path(DIR_FIGURES, "Fig2C_forest_univariate.pdf"), p, width = 6, height = 6)
message("완료: LASSO/Cox signature")
