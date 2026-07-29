# =============================================================================
# 09_trajectory_cellchat.R               [Fig 7]
# Monocle3(SLC7A7 연속값 투영) · Slingshot 교차검증 · CellChat(donor high/low)
# ※ prolif 클러스터는 trajectory 에서 제외 (세포주기 신호가 활성화 연속체를 교란).
#    Methods 에 근거 명시. 비율/CCC 분석에는 포함.
# =============================================================================
source(here::here("R", "setup.R"))
suppressPackageStartupMessages({
  library(Seurat); library(monocle3); library(slingshot); library(SingleCellExperiment)
  library(CellChat) })

mg <- readRDS(file.path(DIR_DATA_PROC, "microglia.rds"))
set.seed(SUBSAMPLE_SEED)
sm <- mg
if ("annotation_level_4" %in% colnames(sm@meta.data))
  sm <- subset(sm, subset = annotation_level_4 != "TAM-MG prolif")
cnt <- as(GetAssayData(sm, assay = "RNA", layer = "counts"), "dgCMatrix")

cds <- new_cell_data_set(cnt, cell_metadata = sm@meta.data,
        gene_metadata = data.frame(gene_short_name = rownames(cnt), row.names = rownames(cnt)))
cds <- preprocess_cds(cds, num_dim = 30)
dn  <- table(colData(cds)$donor_id); cds <- cds[, colData(cds)$donor_id %in% names(dn[dn >= 20])]
colData(cds)$donor_id <- droplevels(as.factor(colData(cds)$donor_id))
cds <- align_cds(cds, alignment_group = "donor_id",
                 alignment_k = min(20, min(table(colData(cds)$donor_id)) - 1))
cds <- reduce_dimension(cds, preprocess_method = "Aligned")
cds <- cluster_cells(cds, resolution = 1e-4)
cds <- learn_graph(cds, use_partition = length(unique(partitions(cds))) > 1, close_loop = FALSE)

hm <- intersect(c("P2RY12","TMEM119","CX3CR1","SALL1","MEF2C"), rownames(cds))
cds$homeo <- colMeans(as.matrix(exprs(cds)[hm, , drop = FALSE]))
root_cl <- names(sort(tapply(cds$homeo, monocle3::clusters(cds), mean), decreasing = TRUE))[1]
cds <- order_cells(cds, root_cells = colnames(cds)[monocle3::clusters(cds) == root_cl])
colData(cds)$pseudotime <- pseudotime(cds)
saveRDS(cds, file.path(DIR_DATA_PROC, "cds_microglia.rds"))

# ★ 연속값 투영 (이분 금지)
ggplot2::ggsave(file.path(DIR_FIGURES, "Fig7A_monocle3_gene_continuous.pdf"),
  plot_cells(cds, genes = GENE_OF_INTEREST, cell_size = .6, label_cell_groups = FALSE,
             label_leaves = FALSE, label_branch_points = FALSE, label_roots = FALSE) +
    ggplot2::scale_color_viridis_c(), width = 6, height = 5)
ggplot2::ggsave(file.path(DIR_FIGURES, "Fig7B_monocle3_pseudotime.pdf"),
  plot_cells(cds, color_cells_by = "pseudotime", cell_size = .6, label_cell_groups = FALSE,
             label_leaves = FALSE, label_branch_points = FALSE, label_roots = FALSE),
  width = 6, height = 5)
ggplot2::ggsave(file.path(DIR_FIGURES, "Fig7D_monocle3_genebypseudotime.pdf"),
                plot_genes_in_pseudotime(cell_size = 1.5, cds[c("SLC7A7"), ], min_expr = 0.1))

# donor 수준 상관 (pseudoreplication 회피)
pt <- data.frame(pt = colData(cds)$pseudotime,
                 expr = as.numeric(exprs(cds)[GENE_OF_INTEREST, ]),
                 donor = colData(cds)$donor_id) %>%
  dplyr::filter(is.finite(pt)) %>% dplyr::group_by(donor) %>%
  dplyr::summarise(pt = mean(pt), expr = mean(expr), .groups = "drop")
print(stats::cor.test(pt$pt, pt$expr, method = "spearman"))
write_result(pt, "Fig7_donor_pseudotime_vs_gene.csv")

## Slingshot 교차검증
sce <- as.SingleCellExperiment(NormalizeData(CreateSeuratObject(cnt, meta.data = sm@meta.data)))
reducedDims(sce)$UMAP <- Embeddings(sm, "umap")[colnames(sce), ]
sce$cluster <- sm$annotation_level_4[colnames(sce)]
sce <- slingshot(sce, clusterLabels = "cluster", reducedDim = "UMAP", start.clus = "TAM-MG aging sig")
saveRDS(slingPseudotime(sce), file.path(DIR_DATA_PROC, "slingshot_pseudotime.rds"))

## CellChat: donor high/low (세포 단위 이분 아님)
seu <- readRDS(file.path(DIR_DATA_PROC, "seurat_gbmap.rds"))
dexp <- data.frame(donor = seu$donor_id, e = FetchData(seu, vars = GENE_OF_INTEREST)[, 1]) %>%
  dplyr::group_by(donor) %>% dplyr::summarise(e = mean(e), .groups = "drop")
hi <- dexp$donor[dexp$e > stats::median(dexp$e)]
seu$donor_group <- ifelse(seu$donor_id %in% hi, "SLC7A7_high", "SLC7A7_low")
set.seed(SUBSAMPLE_SEED)
cells <- unlist(lapply(c("SLC7A7_high","SLC7A7_low"), function(g) {
  cl <- colnames(seu)[seu$donor_group == g]; sample(cl, min(CELLCHAT_SUBSAMPLE, length(cl))) }))
seu_cc <- subset(seu, cells = cells)
# CellChat 실행 전에 — 두 군 모두에서 최소 세포 수를 만족하는 cell type만 유지
min_cells <- 50
tab <- table(seu_cc$cell_type, seu_cc$donor_group)
keep_ct <- rownames(tab)[apply(tab, 1, min) >= min_cells]
message("유지(", length(keep_ct), "): ", paste(keep_ct, collapse = ", "))
message("제외: ", paste(setdiff(rownames(tab), keep_ct), collapse = ", "))
print(tab[setdiff(rownames(tab), keep_ct), , drop = FALSE])   # 몇 개였는지 기록용

seu_cc <- subset(seu_cc, subset = cell_type %in% keep_ct)

# ★ 핵심: character 가 아니라 '레벨을 고정한 factor' 로 지정 → 두 객체 차원이 강제로 일치
seu_cc$cellchat_group <- factor(as.character(seu_cc$cell_type), levels = keep_ct)

cc_hi <- run_cellchat(subset(seu_cc, subset = donor_group == "SLC7A7_high"))
cc_lo <- run_cellchat(subset(seu_cc, subset = donor_group == "SLC7A7_low"))
saveRDS(cc_hi, file.path(DIR_RESULTS, "cellchat_donor_high.rds"))
saveRDS(cc_lo, file.path(DIR_RESULTS, "cellchat_donor_low.rds"))
common.pathways <- intersect(cc_lo@netP$pathways,
                             cc_hi@netP$pathways)
merged <- mergeCellChat(list(SLC7A7_low = cc_lo, SLC7A7_high = cc_hi),
                        add.names = c("SLC7A7_low","SLC7A7_high"), cell.prefix = TRUE)
grDevices::pdf(file.path(DIR_FIGURES, "Fig7C_cellchat_comparison.pdf"), width = 10, height = 6)
print(compareInteractions(merged, show.legend = FALSE, group = c(1, 2)))
print(rankNet(merged, mode = "comparison", stacked = TRUE, do.stat = TRUE, signaling = common.pathways,
              color.use = c("#4B6FD6", "#D64B4B")))
grDevices::dev.off()
message("완료: Fig7")
