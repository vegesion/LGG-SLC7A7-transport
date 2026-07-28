# =============================================================================
# 05_bulk_clinical_gsea_decon.R         [Fig 4]
# 임상 상관 · 상관기반 GSEA · 면역 deconvolution(+myeloid 보정 Cox)
# =============================================================================
source(here::here("R", "setup.R"))
suppressPackageStartupMessages({ library(clusterProfiler); library(msigdbr); library(ggpubr);
  library(survival); library(estimate); library(xCell)})

co <- load_cohort(TRAIN_COHORT); expr <- co$expr; clin <- co$clin; v <- cohort_gene(co)

## (a) 임상 상관
cl <- data.frame(expr = v, idh = clin$idh, grade = clin$grade, codel = clin$codel)
for (vr in c("idh", "grade", "codel")) {
  d <- cl[!is.na(cl[[vr]]), ]; if (nrow(d) < 10 || dplyr::n_distinct(d[[vr]]) < 2) next
  p <- ggpubr::ggviolin(d, x = vr, y = "expr", fill = vr, add = "boxplot",
                        add.params = list(fill = "white", width = .12)) +
    ggpubr::stat_compare_means() + ggplot2::labs(x = NULL, y = paste(GENE_OF_INTEREST, "expression")) +
    theme_paper() + ggplot2::theme(legend.position = "none")
  ggplot2::ggsave(file.path(DIR_FIGURES, paste0("Fig4A_violin_", vr, ".pdf")), p, width = 4.2, height = 4.2)
}
write_result(cl, "Fig4_clinical_association.csv")

## (b) 상관 기반 GSEA
ranked <- correlation_ranked_list(expr, gene = GENE_OF_INTEREST, method = "spearman")
write_result(run_gsea_hallmark(ranked), "Fig4_GSEA_hallmark.csv", row.names = TRUE)
write_result(run_gsego(ranked, ont = "BP"), "Fig4_GSEA_GO_BP.csv")
write_result(run_gsego(ranked, ont = "MF"), "Fig4_GSEA_GO_MF.csv")
write_result(run_gsego(ranked, ont = "CC"), "Fig4_GSEA_GO_CC.csv")


## (c) Deconvolution + myeloid 보정
frac <- tryCatch(run_estimate(expr), error = function(e) { message("ESTIMATE 실패: ", conditionMessage(e)); NULL })
frac$sample <- gsub("\\.", "-", frac$sample)
xc <- run_xcell(expr)
if (!is.null(xc)) frac <- if (is.null(frac)) xc else dplyr::left_join(frac, xc, by = "sample")
if (!is.null(frac)) {
  write_result(frac, "Fig4_deconvolution_scores.csv")
  cors <- correlate_gene_fractions(expr, frac); write_result(cors, "Fig4_gene_vs_fractions.csv")
  print(utils::head(cors, 15))
  mye_cols <- c("Macrophages", "Macrophages M1", "Macrophages M2", "Monocytes",
                "ImmuneScore.x", "MicroenvironmentScore")     # 광범위 지표도 포함
  sens <- dplyr::bind_rows(lapply(intersect(mye_cols, names(frac)), function(cc) {
    a <- adjust_for_fraction(expr, clin, frac, fraction_col = cc)
    g <- a[a$term == "gene", ]
    data.frame(adjusted_for = cc, HR = g$HR, lower = g$lower, upper = g$upper, p = g$p)
  }))
  print(sens)
  write_result(sens, "Fig4_sensitivity_myeloid_adjustment.csv")
  for (mye in mye_cols){
    adj <- adjust_for_fraction(expr, clin, frac, fraction_col = mye)
    print(adj); write_result(adj, paste0("Fig4_cox_adjusted_for_", mye, ".csv"))
  }
  top <- utils::head(cors[order(-abs(cors$rho)), ], 20)
  ggplot2::ggsave(file.path(DIR_FIGURES, "Fig4C_deconvolution_corr.pdf"),
    ggplot2::ggplot(top, ggplot2::aes(rho, stats::reorder(cell_type, rho), fill = FDR < 0.05)) +
      ggplot2::geom_col() +
      ggplot2::scale_fill_manual(values = c(`TRUE` = "#D64B4B", `FALSE` = "grey75")) +
      ggplot2::labs(x = "Spearman rho", y = NULL) + theme_paper(), width = 6, height = 6)
}
message("완료: 임상/GSEA/deconvolution")
