# =============================================================================
# 07_sc_atlas_donor.R                    [Fig 5]
# atlas 개관(donor 수준): UMAP · 세포유형별 발현 · AUCell(수송) · SLC3A2 공발현
# ★ dataset 랜덤효과 + within-dataset 메타분석 (ICC 높음 대응)
# =============================================================================
source(here::here("R", "setup.R"))
suppressPackageStartupMessages({ library(Seurat); library(UCell); library(ggpubr) })

seu <- readRDS(file.path(DIR_DATA_PROC, "seurat_gbmap.rds"))
ggplot2::ggsave(file.path(DIR_FIGURES, "Fig5A_umap_celltype.pdf"),
  DimPlot(seu, group.by = "cell_type", label = TRUE, raster = TRUE) + theme_paper(), width = 7, height = 6)
ggplot2::ggsave(file.path(DIR_FIGURES, "Fig5B_umap_gene.pdf"),
  FeaturePlot(seu, features = GENE_OF_INTEREST, raster = TRUE) + theme_paper(), width = 6, height = 5.5)

## 세포유형별 — donor 수준 요약
df <- data.frame(expr = FetchData(seu, vars = GENE_OF_INTEREST)[, 1],
                 ct = as.character(seu$cell_type), donor = seu$donor_id)
donor_ct <- df %>% dplyr::group_by(ct, donor) %>%
  dplyr::summarise(m = mean(expr), n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n >= SC_DONOR_MIN_CELLS)
write_result(donor_ct, "Fig5_donor_celltype_expression.csv")
ggplot2::ggsave(file.path(DIR_FIGURES, "Fig5C_celltype_donor.pdf"),
  ggpubr::ggboxplot(donor_ct, x = "ct", y = "m", add = "jitter") +
    ggplot2::labs(x = NULL, y = paste(GENE_OF_INTEREST, "(donor mean)")) + theme_paper() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)), width = 7, height = 5)

## AUCell (depth 강건) + donor 연속 검정
tr_set <- build_transporter_geneset(universe = rownames(seu))
auc <- score_auc(seu, list(transport = tr_set,
                           arg_enzyme = build_arg_enzyme_geneset(universe = rownames(seu))))
seu$AUC_transport <- auc[, "transport"]; seu$AUC_argenzyme <- auc[, "arg_enzyme"]

r <- donor_continuous_test(seu, "AUC_transport")
print(r$coef); if (!is.null(r$mixed_coef)) print(r$mixed_coef)
cat(sprintf("within-dataset meta: beta=%.4g, se=%.4g, p=%.3g (%s datasets)\n",
            r$meta$pooled, r$meta$se, r$meta$p, r$meta$n_dataset))
write_result(r$data, "Fig5_donor_AUC_transport.csv")
write_result(r$meta$per_dataset, "Fig5_within_dataset_meta.csv")

## SLC3A2 공발현 (화학량론) — donor 수준 + depth 보정
part <- data.frame(slc7a7 = FetchData(seu, vars = GENE_OF_INTEREST)[, 1],
                   slc3a2 = if (GENE_PARTNER %in% rownames(seu)) FetchData(seu, vars = GENE_PARTNER)[, 1] else NA,
                   depth = seu[[DEPTH_COVARIATE]][, 1], donor = seu$donor_id) %>%
  dplyr::group_by(donor) %>%
  dplyr::summarise(dplyr::across(c(slc7a7, slc3a2, depth), mean), n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n >= SC_DONOR_MIN_CELLS)
print(summary(stats::lm(slc3a2 ~ slc7a7 + depth, data = part))$coefficients)
write_result(part, "Fig5_donor_SLC3A2_coexpression.csv")
ggplot2::ggsave(file.path(DIR_FIGURES, "Fig5D_SLC3A2_donor.pdf"),
  ggplot2::ggplot(part, ggplot2::aes(slc7a7, slc3a2)) + ggplot2::geom_point(size = 2, alpha = .8) +
    ggplot2::geom_smooth(method = "lm") + ggpubr::stat_cor(method = "spearman") +
    ggplot2::labs(x = paste(GENE_OF_INTEREST, "(donor)"), y = paste(GENE_PARTNER, "(donor)")) +
    theme_paper(), width = 4.6, height = 4.4)

saveRDS(seu, file.path(DIR_DATA_PROC, "seurat_gbmap_scored.rds"))
saveRDS(auc, file.path(DIR_DATA_PROC, "score_auc.rds"))
message("완료: Fig5")




"=================="

# add ---------------------------------------------------------------------

#12개 중 7개 양수, 5개 음수(Mathewson2021은 −0.614). 이건 "일관되게 재현됨"이 아니라 이질적입니다.
#최소 donor 수를 8~10으로 올리고 랜덤효과 메타 + I²로 다시 계산하세요.

per <- r$meta$per_dataset
per <- per[per$n >= 8, ]                       # 소표본 dataset 제외
w  <- 1/per$se^2; b <- sum(w*per$beta)/sum(w)
Q  <- sum(w*(per$beta - b)^2); df <- nrow(per) - 1
tau2 <- max(0, (Q - df)/(sum(w) - sum(w^2)/sum(w)))   # DerSimonian–Laird
wr <- 1/(per$se^2 + tau2); br <- sum(wr*per$beta)/sum(wr); ser <- sqrt(1/sum(wr))
c(beta = br, se = ser, p = 2*pnorm(-abs(br/ser)), I2 = max(0,(Q-df)/Q)*100, k = nrow(per))


# 순환논리 점검 -----------------------------------------------------------------

tr_noself <- setdiff(tr_set, GENE_OF_INTEREST)
auc2 <- score_auc(seu, list(transport_ns = tr_noself))
seu$AUC_transport_ns <- auc2[, "transport_ns"]
r2 <- donor_continuous_test(seu, "AUC_transport_ns")
print(r2$mixed_coef)     # 이 값이 진짜 효과
