# =============================================================================
# 05_scrna_preprocess_gbmap.R
# GBmap core atlas(h5ad) → Seurat v5 객체(seu) — BPCells 온디스크(RAM 안전)
#   + obs 메타데이터 전체 + UMAP + LogNormalize + myeloid subset
# 입력: data/raw/<GBmap>.h5ad          출력: data/processed/{seurat_gbmap,microglia,mac,myeloid}.rds
# 원본: GBmap 분석.R (전처리 부분)
# -----------------------------------------------------------------------------
# ★ 대용량 atlas(수십만 세포)를 RAM에 다 올리지 않도록 BPCells 온디스크 행렬 사용.
#   최초 실행 시 h5ad → data/processed/gbmap_bpcells 로 1회 변환하고, 이후 재사용.
#   ※ 이 BPCells 폴더는 지우면 안 됨 — 저장된 seu 가 이 폴더를 참조함.
# =============================================================================
source(here::here("R", "setup.R"))
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(hdf5r); library(BPCells)
  library(AnnotationDbi); library(org.Hs.eg.db)
})

BP    <- file.path(DIR_DATA_PROC, "gbmap_bpcells")          # 온디스크 행렬 폴더
h5ad  <- require_data(basename(PATH_H5AD))                  # 원본 h5ad 경로

# open_matrix_dir + 방향정렬 + SYMBOL + obs 전체 + UMAP + NormalizeData 를 한 번에
seu <- load_gbmap_bpcells(BP, h5ad)
stopifnot(GENE_OF_INTEREST %in% rownames(seu))
print(colnames(seu@meta.data))
print(table(seu$cell_type))

# 주요 subset (원본과 동일 정의)
microglia <- subset(seu, subset = cell_type == "microglial cell")
mac       <- subset(seu, subset = cell_type == "macrophage")
myeloid   <- subset(seu, subset = cell_type %in% c("macrophage", "microglial cell", "monocyte"))

saveRDS(seu,       file.path(DIR_DATA_PROC, "seurat_gbmap.rds"))
saveRDS(microglia, file.path(DIR_DATA_PROC, "microglia.rds"))
saveRDS(mac,       file.path(DIR_DATA_PROC, "mac.rds"))
saveRDS(myeloid,   file.path(DIR_DATA_PROC, "myeloid.rds"))
message("완료: GBmap 전처리 | microglia ", ncol(microglia), " cells")
