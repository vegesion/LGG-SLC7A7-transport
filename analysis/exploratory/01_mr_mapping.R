# =============================================================================
# exploratory/01_mr_mapping.R
# LGG Cox 결과에서 SLC family만 추출 → SYMBOL을 ENSEMBL/ENTREZID로 매핑 → SLC_coxtest.csv
# 입력: data/processed/LGG_all_cox_zscore.csv (02번 산출)
# 출력: data/processed/SLC_coxtest.csv  (MR exposure 목록)
# 원본: "MR위한 SLC mapping 임시 R 파일.R", "MR위한 cox_filtered mapping 임시 R 파일.R"
# =============================================================================
source(here::here("R", "setup.R"))
suppressPackageStartupMessages({ library(org.Hs.eg.db); library(clusterProfiler) })

df <- utils::read.csv(PATH_LGG_COX_ALL)
df_slc <- df[grepl(GENE_FAMILY, df$gene), ]
df_slc <- df_slc[order(df_slc$overall_p), ]
df_slc$ENTREZID <- map_gene_ids(df_slc$gene, from = "SYMBOL", to = "ENTREZID")
df_slc$ENSEMBL  <- map_gene_ids(df_slc$gene, from = "SYMBOL", to = "ENSEMBL")

utils::write.csv(df_slc, PATH_SLC_COXTEST, row.names = FALSE)
message("완료: SLC_coxtest.csv (", nrow(df_slc), " genes)")
