# =============================================================================
# R/gsea.R  ---  GSEA 관련 함수 / Gene set enrichment
# -----------------------------------------------------------------------------
# 세 가지 GSEA 방식(원본 TCGA GSEA.R / GSEA_LGG.R / GSEA_pseudobulk.R)을 통합:
#  (1) correlation 기반 (bulk): 특정 유전자와의 spearman 상관 순위 → GSEA/gseGO
#  (2) DEG logFC 기반 (single-cell): FindMarkers avg_log2FC 순위 → Hallmark GSEA
#  (3) limma-voom t 기반 (pseudobulk): topTable t 순위 → fgsea
# 세 경로 모두 동일 Hallmark(H) gene set을 써야 상호 비교가 가능합니다.
# =============================================================================

# ── Hallmark gene set (TERM2GENE / list) ─────────────────────────────────────
get_hallmark_t2g <- function() {
  msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::select(gs_name, gene_symbol)
}
get_hallmark_list <- function() {
  h <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
  split(h$gene_symbol, h$gs_name)
}

get_gobp_list <- function() {
  h <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
  split(h$gene_symbol, h$gs_name)
}

get_gocc_list <- function() {
  h <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC")
  split(h$gene_symbol, h$gs_name)
}

get_gomf_list <- function() {
  h <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF")
  split(h$gene_symbol, h$gs_name)
}
# ── 발현행렬에서 특정 유전자와의 상관 순위 벡터 / correlation ranking ─────────
correlation_ranked_list <- function(expr_mat, gene = GENE_OF_INTEREST, method = "spearman") {
  gene_vec <- as.numeric(expr_mat[gene, ])
  cor_vec  <- apply(expr_mat, 1, function(x) stats::cor(x, gene_vec, method = method))
  sort(cor_vec, decreasing = TRUE)
}

# ── (1) correlation → Hallmark GSEA (clusterProfiler) ────────────────────────
# 반환: 정렬된 결과 data.frame. 원본 GSEA 객체는 attr(x, "gsea_obj") 로 접근
# (gseaplot2 등 시각화에 필요).
run_gsea_hallmark <- function(ranked_list, seed = TRUE) {
  obj <- clusterProfiler::GSEA(
    geneList = ranked_list, TERM2GENE = get_hallmark_t2g(),
    pvalueCutoff = GSEA_PVAL, minGSSize = GSEA_MIN_SIZE, maxGSSize = GSEA_MAX_SIZE,
    seed = seed, verbose = FALSE
  )
  out <- dplyr::arrange(as.data.frame(obj@result), dplyr::desc(NES))
  attr(out, "gsea_obj") <- obj
  out
}

# ── (1b) correlation → gseGO (SYMBOL keyType) ────────────────────────────────
run_gsego <- function(ranked_list, ont = "BP") {
  res <- clusterProfiler::gseGO(
    geneList = ranked_list, OrgDb = org.Hs.eg.db::org.Hs.eg.db, ont = ont,
    keyType = "SYMBOL", minGSSize = GSEA_MIN_SIZE, maxGSSize = GSEA_MAX_SIZE,
    pvalueCutoff = GSEA_PVAL, pAdjustMethod = "BH", verbose = FALSE,
    seed = TRUE, eps = 1e-10
  )
  dplyr::arrange(as.data.frame(res), dplyr::desc(NES))
}

# ── (2) DEG(FindMarkers) avg_log2FC → Hallmark GSEA ──────────────────────────
ranked_list_from_deg <- function(deg, logfc_col = "avg_log2FC", drop_prefix = c("^MT", "^ENSG")) {
  keep <- rep(TRUE, nrow(deg))
  for (p in drop_prefix) keep <- keep & !grepl(p, rownames(deg))
  deg <- deg[keep, , drop = FALSE]
  gl <- deg[[logfc_col]]; names(gl) <- rownames(deg)
  sort(gl, decreasing = TRUE)
}

# ── (3) limma-voom fit → fgsea (t 통계량 순위) ───────────────────────────────
run_fgsea_pseudobulk <- function(fit, coef = paste0(GENE_OF_INTEREST, "_statusPositive"),
                                 seed = GSEA_SEED) {
  res_df <- limma::topTable(fit, coef = coef, number = Inf) %>%
    tibble::rownames_to_column("gene")
  ranks <- stats::setNames(res_df$t, res_df$gene)
  ranks <- sort(ranks, decreasing = TRUE)
  set.seed(seed)
  gsea_pb <- fgsea::fgsea(pathways = get_hallmark_list(), stats = ranks,
                          minSize = 10, maxSize = GSEA_MAX_SIZE, eps = 0)
  list(ranks = ranks,
       table = gsea_pb[order(gsea_pb$pval), c("pathway", "NES", "pval", "padj", "size", "ES")])
}

# ── NES bar plot (유의 pathway) / NES barplot of significant pathways ────────
plot_nes_bar <- function(gsea_tbl, padj_cut = 0.05, top_n = 10, blank = FALSE) {
  nes_plot <- gsea_tbl %>% dplyr::filter(padj <= padj_cut) %>%
    dplyr::arrange(dplyr::desc(abs(NES))) %>% utils::head(top_n)
  p <- ggplot2::ggplot(nes_plot, ggplot2::aes(x = stats::reorder(pathway, NES), y = NES, fill = NES > 0)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = c("TRUE" = "#D64B4B", "FALSE" = "#4B6FD6")) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "NES", title = "Hallmark pathway enrichment") +
    theme_paper(base_size = 14) + ggplot2::theme(legend.position = "none")
  if (blank) p <- p + theme_blank_frame()
  p
}
