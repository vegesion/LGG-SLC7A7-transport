# =============================================================================
# 12_figures.R
# 논문 제출용 "그림 틀" 생성 (축/제목/텍스트 제거, tick만) → PPT에서 라벨 재작성
#   (a) Volcano  (b) myeloid SLC7A7 violin  (c) Vln/Feature "틀"
# 입력: results/SLC7A7_DEG_microglia_pseudobulk.csv, data/processed/*.rds
# 출력: figures/*.tiff / *.emf
# 원본: volcano.R, seurat그래프그리는코드.R, GBmap 분석.R(violin)
# =============================================================================
source(here::here("R", "setup.R"))
suppressPackageStartupMessages({
  library(readr); library(ggrepel); library(ggpubr); library(Seurat); library(devEMF)
})

# ── (a) Volcano ──────────────────────────────────────────────────────────────
fc_cut <- 0.15; p_cut <- 0.05
df <- readr::read_csv(file.path(DIR_RESULTS, "Fig6_pseudobulk_DEG_limma.csv"), show_col_types = FALSE)
colnames(df)[1] <- "gene"
df <- dplyr::filter(df, !grepl("^ENSG", gene))
df$mlog10P <- -log10(df$adj.P.Val)
df$mlog10P[is.infinite(df$mlog10P)] <- 300
df$Significance <- "Not significant"
df$Significance[df$logFC >=  fc_cut & df$adj.P.Val < p_cut] <- "UP"
df$Significance[df$logFC <= -fc_cut & df$adj.P.Val < p_cut] <- "DOWN"
df$Significance <- factor(df$Significance, levels = c("UP", "DOWN", "Not significant"))

volcano <- ggplot(df, aes(logFC, mlog10P, color = Significance)) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_color_manual(values = c(UP = "#D64B4B", DOWN = "#4B6FD6", "Not significant" = "grey70")) +
  geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = 5, linewidth = 0.6) +
  geom_hline(yintercept = -log10(p_cut), linetype = 5, linewidth = 0.6) +
  labs(x = "log2(FC)", y = "-log10(adj.P)") +
  theme_bw()

save_tiff(volcano, "Fig6_pseudobulk_volcano.tiff", width = 7, height = 7)
save_emf(volcano, "Fig6_pseudobulk_volcano.emf")




# OVA (GO enrichment analysis) --------------------------------------------


# ── (b) myeloid 3종 SLC7A7 violin (통계표 포함) ──────────────────────────────
myeloid <- readRDS(file.path(DIR_DATA_PROC, "myeloid.rds"))   # 05에서 저장
vln_df <- FetchData(myeloid, vars = c(GENE_OF_INTEREST, "cell_type"))
names(vln_df)[1] <- "expr"
summary_tbl <- vln_df %>% dplyr::group_by(cell_type) %>%
  dplyr::summarise(n = dplyr::n(), mean = mean(expr), median = median(expr),
                   sd = sd(expr), IQR = IQR(expr))
write_result(summary_tbl, "myeloid_SLC7A7_summary.csv", row.names = FALSE)

vln <- ggplot(vln_df, aes(cell_type, expr)) +
  geom_violin(trim = FALSE) + geom_boxplot(width = 0.12, outlier.shape = NA) +
  stat_compare_means(comparisons = list(
    c("microglial cell", "macrophage"), c("monocyte", "microglial cell"),
    c("macrophage", "monocyte")), method = "wilcox.test") +
  theme_paper() + labs(x = NULL, y = paste(GENE_OF_INTEREST, "expression"))
save_tiff(vln, "violin_myeloid_SLC7A7.tiff", width = 7, height = 7)
message("완료: figures (volcano + violin, 틀 포함)")
