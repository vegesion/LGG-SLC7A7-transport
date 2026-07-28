# =============================================================================
# 08_sc_myeloid_donor.R                  [Fig 6]
# 아종 UMAP · donor 수준 발현 vs 아종 비율 회귀(핵심) · MG vs MDM · pseudobulk DEG/GSEA
# =============================================================================
source(here::here("R", "setup.R"))
suppressPackageStartupMessages({
  library(Seurat); library(edgeR); library(limma); library(fgsea); library(ggpubr) })

mye <- readRDS(file.path(DIR_DATA_PROC, "myeloid.rds"))
mg  <- readRDS(file.path(DIR_DATA_PROC, "microglia.rds"))
sub_col <- intersect(c("annotation_level_4","celltype_original","annotation_level_3"),
                     colnames(mye@meta.data))[1]

ggplot2::ggsave(file.path(DIR_FIGURES, "Fig6A_umap_myeloid_subtype.pdf"),
  DimPlot(mye, group.by = sub_col, label = TRUE, raster = TRUE) + theme_paper(), width = 8, height = 6)

## donor 수준 발현 vs 아종 비율 (depth 보정 + dataset 혼합모형)
pr <- donor_subtype_proportion(mye, sub_col)
print(pr$result); write_result(pr$result, "Fig6_donor_subtype_proportion_regression.csv")
write_result(pr$data, "Fig6_donor_subtype_proportions.csv")
top_sub <- utils::head(pr$result$subtype[order(pr$result$p)], 4)
ggplot2::ggsave(file.path(DIR_FIGURES, "Fig6B_proportion_regression.pdf"),
  ggplot2::ggplot(dplyr::filter(pr$data, subtype %in% top_sub), ggplot2::aes(gene, prop)) +
    ggplot2::geom_point(size = 1.8, alpha = .8) + ggplot2::geom_smooth(method = "lm") +
    ggpubr::stat_cor(method = "spearman", size = 3) +
    ggplot2::facet_wrap(~subtype, scales = "free_y") +
    ggplot2::labs(x = paste(GENE_OF_INTEREST, "(donor mean)"), y = "Subtype proportion") +
    theme_paper(), width = 7.5, height = 6)

## microglia vs MDM
cmp <- compare_mg_mdm(mye); print(cmp$test); write_result(cmp$data, "Fig6_microglia_vs_MDM.csv")
ggplot2::ggsave(file.path(DIR_FIGURES, "Fig6C_mg_vs_mdm.pdf"),
  ggpubr::ggpaired(tidyr::pivot_longer(cmp$data, c(Microglia, MDM),
                                       names_to = "lineage", values_to = "expr"),
                   x = "lineage", y = "expr", id = "donor", line.color = "grey80") +
    ggpubr::stat_compare_means(paired = TRUE) +
    ggplot2::labs(x = NULL, y = paste(GENE_OF_INTEREST, "(donor mean)")) + theme_paper(),
  width = 4, height = 4.5)

## donor pseudobulk DEG (paired limma) + Hallmark GSEA
mg <- label_slc_status(mg)
pb <- run_pseudobulk_limma(mg); write_result(pb$result, "Fig6_pseudobulk_DEG.csv")
saveRDS(pb, file.path(DIR_DATA_PROC, "pb_fit.rds"))
gs <- run_fgsea_pseudobulk(pb$fit, coef = pb$coef); write_result(gs$table, "Fig6_pseudobulk_GSEA.csv")
print(utils::head(gs$table, 15))
ggplot2::ggsave(file.path(DIR_FIGURES, "Fig6D_pseudobulk_GSEA_NES.pdf"),
                plot_nes_bar(gs$table, top_n = 12), width = 7, height = 5)
message("완료: Fig6")
