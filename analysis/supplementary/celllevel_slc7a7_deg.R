# =============================================================================
# 06_scrna_slc7a7_deg.R
# microglia에서 SLC7A7 +/- 그룹 DEG (Wilcoxon, 전체 유전자) + 정리 + cell-level GSEA
# 입력: data/processed/microglia.rds   출력: results/SLC7A7_DEG_microglia_allgenes.csv
# 원본: GBmap 분석.R (DEG/GSEA 부분)
# =============================================================================
source(here::here("R", "setup.R"))
suppressPackageStartupMessages({ library(Seurat); library(clusterProfiler); library(msigdbr) })

microglia <- readRDS(file.path(DIR_DATA_PROC, "microglia.rds"))
microglia <- label_slc_status(microglia)                       # SLC7A7_status

deg <- find_slc_deg(microglia, all_genes = TRUE)               # 모든 유전자 p 반환
write_result(tibble::rownames_to_column(deg, "gene"),
             paste0(GENE_OF_INTEREST, "_DEG_microglia_allgenes.csv"), row.names = FALSE)

deg_sig   <- dplyr::filter(deg, p_val_adj < 0.05)
deg_clean <- clean_deg_genes(deg_sig, drop_self = GENE_OF_INTEREST)
message("유의 DEG: ", nrow(deg_sig), " → 정리 후: ", nrow(deg_clean))

# cell-level GSEA (참고용; 주 분석은 07 pseudobulk 기반)
ranked <- ranked_list_from_deg(deg, logfc_col = "avg_log2FC")
gsea   <- clusterProfiler::GSEA(ranked, TERM2GENE = get_hallmark_t2g(),
                                pvalueCutoff = GSEA_PVAL, seed = TRUE, verbose = FALSE)
write_result(dplyr::arrange(as.data.frame(gsea@result), dplyr::desc(NES)),
             paste0("GSEA_microglia_celllevel_", GENE_OF_INTEREST, ".csv"), row.names = FALSE)
message("완료: microglia SLC7A7 DEG + GSEA")
