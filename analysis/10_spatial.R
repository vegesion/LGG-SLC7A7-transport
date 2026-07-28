# =============================================================================
# 10_spatial.R                           [Fig 8 — 유일한 차별점]
# 공간 분포 · niche 점수 · spot 수준 SLC7A7 vs ASS1(합성 축) · 저산소 니치 검정
# 입력: data/raw/spatial/<sample>/ (spaceranger outs)
# =============================================================================
source(here::here("R", "setup.R"))
suppressPackageStartupMessages({ library(Seurat); library(ggpubr) })

samples <- load_spatial_all(); cor_all <- list(); niche_all <- list()
for (nm in names(samples)) {
  obj <- add_spatial_modules(samples[[nm]])
  if (!GENE_OF_INTEREST %in% rownames(obj)) { message(nm, ": 유전자 없음"); next }
  ggplot2::ggsave(file.path(DIR_FIGURES, paste0("Fig8_spatial_", nm, "_gene.pdf")),
                  plot_spatial_gene(obj), width = 5.5, height = 5)
  for (s in intersect(c("hypoxia_score","myeloid_score"), names(obj@meta.data)))
    ggplot2::ggsave(file.path(DIR_FIGURES, paste0("Fig8_spatial_", nm, "_", s, ".pdf")),
                    Seurat::SpatialFeaturePlot(obj, features = s), width = 5.5, height = 5)
  cor_all[[nm]] <- spatial_spot_correlation(obj)
  nt <- spatial_niche_test(obj); niche_all[[nm]] <- nt$data
  message(nm, " | 저산소 니치 p = ", signif(nt$test$p.value, 3))
  samples[[nm]] <- obj
}
cor_df <- dplyr::bind_rows(cor_all); write_result(cor_df, "Fig8_spot_correlations.csv"); print(cor_df)
if (length(niche_all)) {
  nd <- dplyr::bind_rows(niche_all); write_result(nd, "Fig8_niche_expression.csv")
  ggplot2::ggsave(file.path(DIR_FIGURES, "Fig8_niche_boxplot.pdf"),
    ggpubr::ggboxplot(nd, x = "group", y = "expr", fill = "group",
                      palette = PALETTE_TWO, facet.by = "sample") +
      ggpubr::stat_compare_means() +
      ggplot2::labs(x = NULL, y = paste(GENE_OF_INTEREST, "(spot)")) + theme_paper(),
    width = 8, height = 4.5)
}
saveRDS(samples, file.path(DIR_DATA_PROC, "spatial_scored.rds"))
message("완료: Fig8")
