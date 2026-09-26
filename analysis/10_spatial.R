# ==============================================================================
# 01_load_visium_h5ad.R
# GBM-space (De Jong et al., bioRxiv 2025) 10X Visium — h5ad 불러오기
#
# 핵심 주의사항 (README 기준):
#   - X 행렬에 gene expression 말고도 3종류의 feature가 섞여 있음.
#     var$feature_type ∈ {"Gene Expression",
#                         "Cell state abundances",       # cell2location deconvolution
#                         "Histopath annotation overlap",# IvyGAP (일부 section만)
#                         "Spatial niche abundances"}    # NMF niche
#   - 따라서 raw count 쓰려면 반드시 feature_type == "Gene Expression" 으로 subset.
#   - obs index = "<barcode>_<sample ID>" → 파일 간 colname 충돌 없음.
# ==============================================================================

## ---- 0. 패키지 --------------------------------------------------------------
# BiocManager::install(c("zellkonverter", "SingleCellExperiment", "rhdf5", "HDF5Array"))
# install.packages(c("Matrix", "dplyr"))
source(here::here("R", "setup.R"))
suppressPackageStartupMessages({
  library(zellkonverter)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(rhdf5)
  library(Matrix)
})
BiocManager::install("zellkonverter")

## ---- 1. 경로 설정 -----------------------------------------------------------
DATA_DIR <- "C:\\Users\\helis\\Desktop\\LGG-SLC7A7-transport\\data\\raw\\spatial_data_visium"   # <-- 수정
OUT_DIR  <- "data/processed/visium"                  # 파일별 RDS 저장 위치
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

h5ad_files <- list.files(DATA_DIR, pattern = "\\.h5ad$", full.names = TRUE)
names(h5ad_files) <- sub("\\.h5ad$", "", basename(h5ad_files))

cat("발견된 h5ad 파일:", length(h5ad_files), "개\n")
print(names(h5ad_files))


## ---- 2. 통째로 읽기 전에 훑어보기 -------------------------------------------

# 2-1. var의 feature type 컬럼 이름을 자동 탐지 (feature_types / feature_type)
detect_ft_key <- function(path) {
  on.exit(rhdf5::h5closeAll(), add = TRUE)
  vars <- h5ls(path)
  vars <- vars$name[vars$group == "/var"]
  hit  <- intersect(c("feature_types", "feature_type"), vars)
  if (length(hit) == 0) NA_character_ else hit[1]
}

FT_KEY <- detect_ft_key(h5ad_files[1])
cat("feature type 컬럼:", FT_KEY, "\n")   # -> "feature_types"


# 2-2. 파일 전체를 안 읽고 feature_types 분포만 확인
peek_feature_types <- function(path, key = FT_KEY) {
  on.exit(rhdf5::h5closeAll(), add = TRUE)
  tryCatch({
    # anndata >= 0.8 : categorical은 categories/codes 그룹
    cats  <- as.character(h5read(path, paste0("var/", key, "/categories")))
    codes <- as.integer(h5read(path, paste0("var/", key, "/codes")))
    table(cats[codes + 1L])
  }, error = function(e) {
    tryCatch(table(as.character(h5read(path, paste0("var/", key)))),
             error = function(e2) NA)
  })
}

peek_feature_types(h5ad_files[1])
# 예상: Gene Expression ~36600, 나머지 100여 개가 cell state / niche / histopath

# 2-3. 30개 파일이 같은 feature 구성인지 (README상 histopath는 일부 section만)
ft_summary <- lapply(h5ad_files, peek_feature_types)
do.call(rbind, lapply(ft_summary, function(x) {
  if (all(is.na(x))) return(NULL)
  as.data.frame(t(as.matrix(x)))
}))

# 2-4. uns/spatial 안에 뭐가 들어있는지 (H&E 이미지 크기 가늠용)
h5ls(h5ad_files[1], recursive = 4)[
  grepl("^/uns", h5ls(h5ad_files[1], recursive = 4)$group), ]


## ---- 3. 로더 함수 -----------------------------------------------------------
FT_MAP <- c(
  "Cell state abundances"        = "cell_state",
  "Spatial niche abundances"     = "niche",
  "Histopath annotation overlap" = "histopath"
)

# uns(H&E 이미지)는 기본적으로 건너뜀 — 30개 파일이면 이것만으로 수 GB.
# 나중에 이미지 필요하면 keep_uns = TRUE 또는 schard 사용.
.read_h5ad <- function(path, keep_uns = FALSE) {
  args <- list(file = path, X_name = "counts", reader = "R", verbose = FALSE)
  if (!keep_uns && "uns" %in% names(formals(zellkonverter::readH5AD))) {
    args$uns <- FALSE
  }
  do.call(zellkonverter::readH5AD, args)
}

load_visium_h5ad <- function(path,
                             gene_feature   = "Gene Expression",
                             ft_key         = FT_KEY,
                             filter_tissue  = TRUE,
                             keep_uns       = FALSE,
                             verbose        = TRUE) {
  
  # reader = "R" : rhdf5 네이티브. Windows에서 basilisk/conda 안 깔아도 됨.
  #                실패하면 args에서 reader = "python" 으로 교체.
  sce <- .read_h5ad(path, keep_uns = keep_uns)
  
  ft <- as.character(rowData(sce)[[ft_key]])
  if (is.null(ft)) stop("var$", ft_key, " 없음: ", basename(path))
  
  # 3-1. gene expression만 메인 assay로
  sce_gene <- sce[ft == gene_feature, ]
  assayNames(sce_gene) <- "counts"
  
  # 3-2. 나머지 feature block은 altExp로 (count 아님 → assay명 'value')
  for (blk in setdiff(unique(ft), gene_feature)) {
    sub_sce <- sce[ft == blk, ]
    assayNames(sub_sce) <- "value"
    nm_alt <- if (blk %in% names(FT_MAP)) FT_MAP[[blk]] else make.names(blk)
    altExp(sce_gene, nm_alt) <- sub_sce
  }
  
  # 3-3. spatial 좌표 (obsm/spatial → reducedDim)
  rd <- reducedDimNames(sce_gene)
  sp <- grep("spatial", rd, ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(sp)) {
    colnames(reducedDim(sce_gene, sp)) <- c("x", "y")
    sce_gene$spatial_x <- reducedDim(sce_gene, sp)[, 1]
    sce_gene$spatial_y <- reducedDim(sce_gene, sp)[, 2]
  } else {
    warning("spatial 좌표 못 찾음: ", basename(path))
  }
  
  # 3-4. tissue 밖 spot 제거
  if (filter_tissue && "in_tissue" %in% colnames(colData(sce_gene))) {
    n0 <- ncol(sce_gene)
    sce_gene <- sce_gene[, sce_gene$in_tissue == 1]
    if (verbose && ncol(sce_gene) < n0) {
      cat(sprintf("  in_tissue 필터: %d -> %d spots\n", n0, ncol(sce_gene)))
    }
  }
  
  # 3-5. 파일명 파싱
  #   AT3-BRA5-FO-1_1   -> tumor=AT3, site=BRA5-FO-1,  rep=1
  #   AT5-BRA-5-FO-1_1  -> tumor=AT5, site=BRA-5-FO-1, rep=1
  nm <- sub("\\.h5ad$", "", basename(path))
  sce_gene$file_id   <- nm
  sce_gene$tumor_id  <- sub("^(AT[0-9]+).*$", "\\1", nm)
  sce_gene$replicate <- sub("^.*_([0-9]+)$", "\\1", nm)
  sce_gene$site      <- sub("_[0-9]+$", "", sub("^AT[0-9]+-", "", nm))
  
  if (verbose) {
    cat(sprintf("  %s | genes: %d | spots: %d | altExp: %s\n",
                nm, nrow(sce_gene), ncol(sce_gene),
                paste(altExpNames(sce_gene), collapse = ", ")))
  }
  sce_gene
}


## ---- 4. 한 개로 먼저 검증 ---------------------------------------------------
sce1 <- load_visium_h5ad(h5ad_files[[1]])

dim(sce1)                                   # 유전자 x spot
assayNames(sce1); altExpNames(sce1)
head(colData(sce1))
head(reducedDim(sce1, "spatial"))



"SLC7A7" %in% rownames(sce1)                # 타겟 존재 확인

# counts가 진짜 정수인지 (FLOAT로 저장돼 있으므로 확인 필요)
cx <- counts(sce1)[, 1:min(50, ncol(sce1))]
all(cx@x == round(cx@x))                    # TRUE 여야 정상
rm(cx)

# cell2location cell state 목록 — microglia/TAM substate 라벨 확인
rownames(altExp(sce1, "cell_state"))
# NMF spatial niche 목록
rownames(altExp(sce1, "niche"))
# histopath (있는 경우)
if ("histopath" %in% altExpNames(sce1)) rownames(altExp(sce1, "histopath"))


## ---- 5. 전체 배치 로드 (파일별 RDS 저장) ------------------------------------
log_fail <- character(0)

for (i in seq_along(h5ad_files)) {
  nm  <- names(h5ad_files)[i]
  out <- file.path(OUT_DIR, paste0(nm, "_sce.rds"))
  if (file.exists(out)) { cat("[skip] ", nm, "\n"); next }
  cat(sprintf("[%d/%d] %s\n", i, length(h5ad_files), nm))
  
  ok <- tryCatch({
    x <- load_visium_h5ad(h5ad_files[[i]])
    saveRDS(x, out); rm(x); gc(verbose = FALSE); TRUE
  }, error = function(e) { message("  실패: ", conditionMessage(e)); FALSE })
  
  if (!ok) log_fail <- c(log_fail, nm)
}

cat("실패:", if (length(log_fail) == 0) "없음" else paste(log_fail, collapse = ", "), "\n")


## ---- 6. 다시 불러오기 / 합치기 ----------------------------------------------
rds_files <- list.files(OUT_DIR, pattern = "_sce\\.rds$", full.names = TRUE)
names(rds_files) <- sub("_sce\\.rds$", "", basename(rds_files))
sce_list <- lapply(rds_files, readRDS)

common_genes <- Reduce(intersect, lapply(sce_list, rownames))
cat("공통 유전자:", length(common_genes), "\n")
rm(sce_list)
sce_all <- do.call(cbind, lapply(sce_list, function(x) {
  y <- x[common_genes, ]; altExps(y) <- NULL; y
}))
rm(sce_list)
sce_all
saveRDS(sce_all, file.path(OUT_DIR, "visium_all_merged_sce.rds"))

sce_all <- readRDS(file.path(OUT_DIR, "visium_all_merged_sce.rds"))

# ==============================================================================
# 옵션 B: H&E 이미지까지 필요하면 (Seurat spatial 객체)
# ------------------------------------------------------------------------------
remotes::install_github("cellgeni/schard")

srt <- schard::h5ad2seurat_spatial(h5ad_files[[1]])
ft  <- srt[["Spatial"]]@meta.data$feature_types      # 복수형 주의
srt <- subset(srt, features = rownames(srt)[ft == "Gene Expression"])
srt <- NormalizeData(
  srt,
  assay = DefaultAssay(srt)
)
Seurat::SpatialFeaturePlot(srt, features = "SLC7A7")

# ==============================================================================
"================================================="

## ============================================================================
## 7. 패키지 · 상수
## ============================================================================
suppressPackageStartupMessages({
  library(scuttle); library(dplyr); library(tidyr); library(tibble)
  library(purrr);   library(ggplot2); library(ggrepel); library(patchwork)
})

theme_set(theme_bw(base_size = 11))
FIG_DIR <- "figures/spatial"; dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
RES_DIR <- "results/spatial";  dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)

TARGET_GENE <- "SLC7A7"

TAM_STATES <- c(
  "Resident-TAMs", "Resident BAM TAMs",
  "Pro-inflammatory TAMs", "Anti-inflammatory TAMs",
  "Interferon TAMs", "Angiogenic TAMs",
  "Stress-response TAMs", "Proliferative TAMs",
  "Astrocyte-like TAMs", "RTN1+ TAMs"
)
MYELOID_STATES <- c(TAM_STATES, "Monocytes", "Dendritic cells")

# 같이 뽑아둘 marker (sanity check + figure용)
GENES_KEEP <- c(
  TARGET_GENE,
  "PTPRC", "AIF1", "CSF1R", "CD68",            # pan-myeloid
  "P2RY12", "TMEM119", "CX3CR1",               # homeostatic microglia
  "CD163", "MRC1", "MSR1", "SPP1", "APOE",     # M2/anti-infl
  "TNF", "IL1B", "NFKBIA", "CCL3", "CCL4",     # pro-infl / TNFA-NFKB
  "ISG15", "IFI6", "MX1", "STAT1", "IFIT3",    # IFN response
  "VEGFA", "HIF1A", "CA9",                     # hypoxia/angiogenesis
  "GFAP", "EGFR", "SOX2", "MKI67"              # tumor/AC
)


## ============================================================================
## 8. spot 단위 tidy table 만들기 (이후 모든 분석의 기반)
## ============================================================================
# altExp의 라벨은 공백/괄호/+ 가 섞여 있어 컬럼명으로 쓰기 나쁘다.
# 안전한 컬럼명 <-> 원 라벨 매핑표를 먼저 만들어 둔다.
make_key <- function(labels, prefix) {
  tibble(label = labels, col = paste0(prefix, make.names(labels)))
}

sce_ref  <- readRDS(rds_files[1])
CS_KEY   <- make_key(rownames(altExp(sce_ref, "cell_state")), "cs_")
NIC_KEY  <- make_key(rownames(altExp(sce_ref, "niche")),      "nic_")
HP_KEY   <- if ("histopath" %in% altExpNames(sce_ref))
  make_key(rownames(altExp(sce_ref, "histopath")), "hp_") else NULL
rm(sce_ref)

TAM_COLS <- CS_KEY$col[CS_KEY$label %in% TAM_STATES]
MYE_COLS <- CS_KEY$col[CS_KEY$label %in% MYELOID_STATES]
# Ambiguous / Undefined 는 주 분석에서 제외 (민감도 분석 때만 사용)
JUNK_COLS <- CS_KEY$col[grepl("^(Ambiguous|Undefined)", CS_KEY$label)]
BONA_COLS <- setdiff(CS_KEY$col, JUNK_COLS)


build_spot_df <- function(sce, genes = GENES_KEEP) {
  # 빈 spot 제거 후 library size 정규화
  sce <- sce[, colSums(counts(sce)) > 0]
  sce <- scuttle::logNormCounts(sce)
  
  g    <- intersect(genes, rownames(sce))
  expr <- as.data.frame(as.matrix(t(logcounts(sce)[g, , drop = FALSE])))
  names(expr) <- paste0("g_", make.names(g))
  
  grab <- function(name, key) {
    if (!name %in% altExpNames(sce)) return(NULL)
    m <- as.data.frame(as.matrix(t(assay(altExp(sce, name), "value"))))
    names(m) <- key$col[match(names(m), key$label)]
    m
  }
  cs  <- grab("cell_state", CS_KEY)
  nic <- grab("niche",      NIC_KEY)
  hp  <- if (!is.null(HP_KEY)) grab("histopath", HP_KEY) else NULL
  
  meta <- as.data.frame(colData(sce)) %>%
    transmute(file_id, tumor_id, site, replicate,
              spatial_x, spatial_y,
              n_counts   = colSums(counts(sce)),
              n_detected = colSums(counts(sce) > 0))
  
  out <- bind_cols(tibble(spot_id = colnames(sce)), meta, expr, cs, nic, hp)
  
  # 조성 요약값
  out$tam_total  <- rowSums(out[, intersect(TAM_COLS, names(out)), drop = FALSE])
  out$mye_total  <- rowSums(out[, intersect(MYE_COLS, names(out)), drop = FALSE])
  out$cell_total <- rowSums(out[, intersect(BONA_COLS, names(out)), drop = FALSE])
  out$tam_frac   <- out$tam_total / pmax(out$cell_total, 1e-6)
  
  # dominant niche (argmax)
  nic_cols <- intersect(NIC_KEY$col, names(out))
  out$dom_niche <- NIC_KEY$label[match(nic_cols[max.col(as.matrix(out[, nic_cols]),
                                                        ties.method = "first")],
                                       NIC_KEY$col)]
  out
}

spot_df <- map_dfr(seq_along(rds_files), function(i) {
  cat(sprintf("[%d/%d] %s\n", i, length(rds_files), names(rds_files)[i]))
  x <- readRDS(rds_files[i]); d <- build_spot_df(x); rm(x); gc(verbose = FALSE); d
})

dim(spot_df)
saveRDS(spot_df, file.path(RES_DIR, "spot_df.rds"))
spot_df <- readRDS(file.path(RES_DIR, "spot_df.rds"))

## ============================================================================
## 9. QC — 이 데이터로 SLC7A7 얘기를 할 수 있는지부터 확인
## ============================================================================
# 9-1. donor(tumor) 수 = 실질적 반복 수. 여기서 통계 설계가 결정된다.
spot_df %>% dplyr::count(tumor_id, file_id) %>% dplyr::count(tumor_id, name = "n_sections")

# 9-2. SLC7A7 검출률. 너무 낮으면(<10~15%) spot-level 상관은 힘이 없다.
qc_det <- spot_df %>%
  group_by(tumor_id, file_id) %>%
  summarise(n_spot     = n(),
            det_rate   = mean(g_SLC7A7 > 0),
            mean_expr  = mean(g_SLC7A7),
            med_counts = median(n_counts),
            mean_tam   = mean(tam_frac), .groups = "drop") %>%
  arrange(det_rate)
print(qc_det, n = 40)

p_qc <- ggplot(qc_det, aes(reorder(file_id, det_rate), det_rate, fill = tumor_id)) +
  geom_col() + coord_flip() +
  labs(x = NULL, y = "SLC7A7 detection rate (spots)", fill = "Tumour")
ggsave(file.path(FIG_DIR, "QC_SLC7A7_detection.pdf"), p_qc, width = 7, height = 7)

# 9-3. SLC7A7 vs total TAM — 교란의 크기를 눈으로 확인
p_conf <- spot_df %>%
  ggplot(aes(tam_frac, g_SLC7A7)) +
  geom_point(size = .2, alpha = .15) +
  geom_smooth(method = "lm", se = FALSE, colour = "firebrick") +
  facet_wrap(~ tumor_id) +
  labs(x = "TAM fraction (cell2location)", y = "SLC7A7 (logcounts)",
       subtitle = "이 관계가 강할수록 보정 없는 상관 결과는 전부 자명해진다")
ggsave(file.path(FIG_DIR, "QC_SLC7A7_vs_TAMfraction.pdf"), p_conf, width = 8, height = 6)

## ============================================================================
## 9-4. 진단 — 어떤 컬럼이, 어느 section에서 비어 있나
## ============================================================================
cs_present <- intersect(CS_KEY$col, names(spot_df))

# (a) spot_df 전체에서 NA가 있는 cell state 컬럼
na_report <- spot_df %>%
  summarise(across(all_of(cs_present), ~ sum(is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "col", values_to = "n_na") %>%
  filter(n_na > 0) %>%
  left_join(CS_KEY, by = "col") %>%
  arrange(desc(n_na))
print(na_report, n = 50)

# (b) 어느 section에 없는지 — tumor별로 갈리면 donor-matched reference 문제
if (nrow(na_report) > 0) {
  spot_df %>%
    group_by(tumor_id, file_id) %>%
    summarise(across(all_of(na_report$col), ~ mean(is.na(.x))), .groups = "drop") %>%
    pivot_longer(-c(tumor_id, file_id), names_to = "col", values_to = "frac_na") %>%
    filter(frac_na > 0) %>%
    left_join(CS_KEY, by = "col") %>%
    dplyr::count(tumor_id, label) %>%
    print(n = 100)
}

# (c) 모든 section에 공통으로 존재하는 컬럼만 남긴다
COMPLETE_COLS <- cs_present[colSums(is.na(spot_df[, cs_present, drop = FALSE])) == 0]
cat(sprintf("cell state: 전체 %d개 중 모든 section 공통 %d개\n",
            length(cs_present), length(COMPLETE_COLS)))

BONA_COLS <- intersect(BONA_COLS, COMPLETE_COLS)
TAM_COLS  <- intersect(TAM_COLS,  COMPLETE_COLS)

# TAM state가 빠진 게 있으면 tam_total 재계산 필요 (section별 정의가 달라지므로)
dropped_tam <- setdiff(CS_KEY$label[CS_KEY$label %in% TAM_STATES],
                       CS_KEY$label[match(TAM_COLS, CS_KEY$col)])
if (length(dropped_tam) > 0) {
  warning("일부 section에 없는 TAM state: ", paste(dropped_tam, collapse = ", "))
  spot_df$tam_total  <- rowSums(spot_df[, TAM_COLS, drop = FALSE])
  spot_df$cell_total <- rowSums(spot_df[, BONA_COLS, drop = FALSE])
  spot_df$tam_frac   <- spot_df$tam_total / pmax(spot_df$cell_total, 1e-6)
}


## ============================================================================
## 10. 분석 1 — total TAM 보정 partial correlation  [수정판]
## ============================================================================

# Spearman partial correlation 닫힌 형태:
#   r_xy.z = (r_xy - r_xz * r_yz) / sqrt((1 - r_xz^2)(1 - r_yz^2))
section_pcor <- function(d, state_cols, target = "g_SLC7A7", ctrl = "tam_total") {
  
  cols <- unique(c(target, ctrl, intersect(state_cols, names(d))))
  if (!all(c(target, ctrl) %in% cols)) return(tibble())
  
  M <- as.matrix(d[, cols, drop = FALSE])
  storage.mode(M) <- "double"
  
  # 핵심 수정: sd가 NA인 열(전부 NA)도 함께 제거.
  # isTRUE()로 감싸서 논리 첨자에 NA가 절대 들어가지 않게 한다.
  ok <- vapply(seq_len(ncol(M)), function(j) {
    s <- stats::sd(M[, j], na.rm = TRUE)
    isTRUE(is.finite(s) && s > 0)
  }, logical(1))
  M <- M[, ok, drop = FALSE]
  
  if (!all(c(target, ctrl) %in% colnames(M))) return(tibble())
  
  R  <- suppressWarnings(stats::cor(M, method = "spearman",
                                    use = "pairwise.complete.obs"))
  st <- setdiff(colnames(M), c(target, ctrl))
  st <- st[!is.na(st) & st %in% colnames(R)]          # 2차 방어
  if (length(st) == 0) return(tibble())
  
  rxy <- R[target, st]
  rxz <- as.numeric(R[target, ctrl])
  ryz <- R[ctrl, st]
  
  den <- sqrt(pmax((1 - rxz^2) * (1 - ryz^2), .Machine$double.eps))
  
  tibble(col     = st,
         rho_raw = as.numeric(rxy),
         rho_adj = as.numeric((rxy - rxz * ryz) / den),
         n_spot  = nrow(M))
}

cor_sec <- spot_df %>%
  group_by(tumor_id, file_id) %>%
  group_modify(~ section_pcor(.x, intersect(BONA_COLS, names(.x)))) %>%
  ungroup()

stopifnot(nrow(cor_sec) > 0)

cor_sec <- cor_sec %>%
  left_join(CS_KEY, by = "col") %>%
  mutate(is_TAM = label %in% TAM_STATES)

# 몇 개 section이 실제로 계산됐는지 확인 (조용히 빠진 게 없는지)
cor_sec %>% distinct(tumor_id, file_id) %>% count(tumor_id)
cat("계산된 section 수:", n_distinct(cor_sec$file_id),
    "/ 전체:", n_distinct(spot_df$file_id), "\n")


# 10-1. tumor 수준으로 한 단계 올린 뒤 집계 (nested 구조 반영)
cor_tumor <- cor_sec %>%
  group_by(label, is_TAM, tumor_id) %>%
  summarise(rho_adj = mean(rho_adj, na.rm = TRUE), .groups = "drop")

meta_cor <- cor_tumor %>%
  group_by(label, is_TAM) %>%
  summarise(mean_rho = mean(rho_adj, na.rm = TRUE),
            sd_rho   = sd(rho_adj,  na.rm = TRUE),
            n_tumor  = sum(is.finite(rho_adj)),
            n_pos    = sum(rho_adj > 0, na.rm = TRUE),
            .groups  = "drop") %>%
  mutate(se = sd_rho / sqrt(n_tumor),
         t  = mean_rho / se,
         p  = 2 * pt(-abs(t), df = n_tumor - 1),
         fdr = p.adjust(p, "BH")) %>%
  arrange(desc(mean_rho))

print(meta_cor %>% filter(is_TAM), n = 20)
write.csv(meta_cor, file.path(RES_DIR, "partial_cor_SLC7A7_cellstate.csv"), row.names = FALSE)

# 10-2. Figure: TAM substate forest plot (보정 전/후 대비)
fp_dat <- cor_sec %>%
  filter(is_TAM) %>%
  pivot_longer(c(rho_raw, rho_adj), names_to = "type", values_to = "rho") %>%
  mutate(type = factor(type, c("rho_raw", "rho_adj"),
                       c("Unadjusted", "Adjusted for total TAM")))

p_forest <- ggplot(fp_dat, aes(rho, reorder(label, rho), colour = tumor_id)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_point(position = position_jitter(height = .15), size = 1.4, alpha = .8) +
  stat_summary(fun = mean, geom = "point", colour = "black",
               shape = 18, size = 3.2, aes(group = 1)) +
  facet_wrap(~ type) +
  labs(x = expression(Spearman~rho~"(spot-level, per section)"),
       y = NULL, colour = "Tumour",
       title = "SLC7A7 expression vs TAM substate abundance")
ggsave(file.path(FIG_DIR, "F1_TAMsubstate_partial_correlation.pdf"),
       p_forest, width = 10, height = 5.5)

# 10-3. Figure: section × cell state heatmap
hm_dat <- cor_sec %>% filter(is_TAM)
p_hm <- ggplot(hm_dat, aes(file_id, reorder(label, rho_adj), fill = rho_adj)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = expression(rho[adj])) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = 7))
ggsave(file.path(FIG_DIR, "F2_TAMsubstate_heatmap.pdf"), p_hm, width = 10, height = 4.5)


## ============================================================================
## 11. 분석 2 — TAM 내부 조성 (compositional). 여기가 GBmap 결과와 직접 대응
## ============================================================================
# TAM이 충분히 있는 spot만. 임계값은 데이터 보고 조정.
TAM_MIN <- quantile(spot_df$tam_total, 0.60)
cat("TAM_MIN =", round(TAM_MIN, 3), "\n")

tam_df <- spot_df %>%
  filter(tam_total >= TAM_MIN) %>%
  mutate(across(all_of(TAM_COLS), ~ .x / tam_total, .names = "p_{.col}"))
P_COLS <- paste0("p_", TAM_COLS)

# 11-1. section 내 SLC7A7 3분위 → TAM 조성 비교
tam_ter <- tam_df %>%
  group_by(file_id) %>%
  mutate(slc_tert = factor(
    ntile(g_SLC7A7, 3),
    levels = 1:3,
    labels = c("Low", "Mid", "High")
  )
  ) %>%
  ungroup() %>%
  filter(!is.na(slc_tert))

comp_sec <- tam_ter %>%
  group_by(tumor_id, file_id, slc_tert) %>%
  summarise(across(all_of(P_COLS), mean), .groups = "drop") %>%
  pivot_longer(all_of(P_COLS), names_to = "col", values_to = "prop") %>%
  mutate(col = sub("^p_", "", col)) %>%
  left_join(CS_KEY, by = "col")

# High vs Low 짝지어 비교 (section이 반복단위)
delta <- comp_sec %>%
  filter(slc_tert %in% c("Low", "High")) %>%
  pivot_wider(names_from = slc_tert, values_from = prop) %>%
  mutate(diff = High - Low) %>%
  group_by(label) %>%
  summarise(mean_diff = mean(diff, na.rm = TRUE),
            n_sec     = n(),
            n_up      = sum(diff > 0, na.rm = TRUE),
            p         = tryCatch(wilcox.test(High, Low, paired = TRUE)$p.value,
                                 error = function(e) NA_real_),
            .groups = "drop") %>%
  mutate(fdr = p.adjust(p, "BH")) %>%
  arrange(desc(mean_diff))
print(delta, n = 20)
write.csv(delta, file.path(RES_DIR, "TAMcomposition_SLC7A7high_vs_low.csv"), row.names = FALSE)

# 11-2. Figure: stacked bar + paired slope
p_stack <- comp_sec %>%
  group_by(slc_tert, label) %>% summarise(prop = mean(prop), .groups = "drop") %>%
  ggplot(aes(slc_tert, prop, fill = label)) +
  geom_col(width = .7) +
  labs(x = "SLC7A7 tertile (within section)", y = "Fraction of TAM compartment",
       fill = "TAM substate")

p_slope <- comp_sec %>%
  filter(slc_tert %in% c("Low", "High")) %>%
  ggplot(aes(slc_tert, prop, group = file_id, colour = tumor_id)) +
  geom_line(alpha = .5) + geom_point(size = 1) +
  facet_wrap(~ label, scales = "free_y", ncol = 5) +
  labs(x = NULL, y = "Fraction of TAM compartment", colour = "Tumour")

ggsave(file.path(FIG_DIR, "F3_TAMcomposition.pdf"),
       p_stack / p_slope + plot_layout(heights = c(1, 2)), width = 12, height = 10)


## ============================================================================
## 12. 분석 3 — niche 수준 (Immune (TAMs) vs Immune (resident))
## ============================================================================
niche_sec <- spot_df %>%
  group_by(tumor_id, file_id, dom_niche) %>%
  summarise(mean_slc = mean(g_SLC7A7), n = n(), .groups = "drop") %>%
  filter(n >= 20)                       # spot 너무 적은 niche 제외

niche_summary <- niche_sec %>%
  group_by(dom_niche) %>%
  summarise(mean_slc = mean(mean_slc), sd = sd(mean_slc),
            n_sec = n(), .groups = "drop") %>%
  arrange(desc(mean_slc))
print(niche_summary, n = 20)

# 두 immune niche 직접 대조 (같은 section 내 짝비교)
imm_pair <- niche_sec %>%
  filter(dom_niche %in% c("Immune (TAMs)", "Immune (resident)")) %>%
  select(tumor_id, file_id, dom_niche, mean_slc) %>%
  pivot_wider(names_from = dom_niche, values_from = mean_slc) %>%
  filter(complete.cases(.))
if (nrow(imm_pair) >= 5) print(wilcox.test(imm_pair[["Immune (TAMs)"]],
                                           imm_pair[["Immune (resident)"]], paired = TRUE))

p_niche <- ggplot(niche_sec, aes(reorder(dom_niche, mean_slc), mean_slc)) +
  geom_boxplot(outlier.shape = NA, fill = "grey92") +
  geom_jitter(aes(colour = tumor_id), width = .15, size = 1.6, alpha = .85) +
  coord_flip() +
  labs(x = NULL, y = "Mean SLC7A7 (logcounts) per section", colour = "Tumour",
       title = "SLC7A7 across NMF spatial niches")
ggsave(file.path(FIG_DIR, "F4_SLC7A7_by_niche.pdf"), p_niche, width = 8, height = 5)


## ============================================================================
## 13. 분석 4 — IvyGAP histopath 영역 (있는 section만)
## ============================================================================
if (!is.null(HP_KEY) && any(HP_KEY$col %in% names(spot_df))) {
  hp_cols <- intersect(HP_KEY$col, names(spot_df))
  hp_df <- spot_df %>%
    filter(rowSums(across(all_of(hp_cols))) > 0)     # annotation 있는 spot만
  
  hp_sec <- hp_df %>%
    group_by(tumor_id, file_id) %>%
    group_modify(~ {
      d <- .x
      map_dfr(hp_cols, function(cc) tibble(
        col = cc,
        rho = suppressWarnings(cor(d$g_SLC7A7, d[[cc]], method = "spearman")),
        rho_adj = { R <- suppressWarnings(cor(cbind(d$g_SLC7A7, d$tam_total, d[[cc]]),
                                              method = "spearman"))
        (R[1,3] - R[1,2]*R[2,3]) / sqrt((1-R[1,2]^2)*(1-R[2,3]^2)) }))
    }) %>% ungroup() %>% left_join(HP_KEY, by = "col")
  
  p_hp <- ggplot(hp_sec, aes(rho_adj, reorder(label, rho_adj), colour = tumor_id)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_point(size = 2, alpha = .85) +
    labs(x = expression(rho[adj]~"(TAM-adjusted)"), y = NULL, colour = "Tumour",
         title = "SLC7A7 vs IvyGAP histopathology annotation")
  ggsave(file.path(FIG_DIR, "F5_SLC7A7_by_histopath.pdf"), p_hp, width = 8, height = 4.5)
  write.csv(hp_sec, file.path(RES_DIR, "SLC7A7_histopath.csv"), row.names = FALSE)
} else {
  message("histopath feature가 이 subset에 없음 — section 13 건너뜀")
}


## ============================================================================
## 14. Figure — 대표 section 공간 지도
## ============================================================================
# Visium 좌표는 이미지 좌표계(위→아래로 y 증가) → scale_y_reverse() 필요
plot_spatial <- function(d, col, title = NULL, pal = "magma", trans = "identity") {
  ggplot(d, aes(spatial_x, spatial_y, colour = .data[[col]])) +
    geom_point(size = .8) + coord_fixed() + scale_y_reverse() +
    scale_colour_viridis_c(option = pal, trans = trans) +
    labs(title = title, colour = NULL) +
    theme_void(base_size = 10) +
    theme(plot.title = element_text(size = 10, hjust = .5))
}
# library(purrr)
# SLC7A7 검출률이 가장 높은 상위 5개 section을 대표로
rep_id <- qc_det$file_id[order(qc_det$det_rate, decreasing = TRUE)[1:5]]
# print(CS_KEY,100)
show_states <- c(
  "Hypoxic 1 (cell state)",
  "Hypoxic 2 (cell state)",
  "Endothelial (capillary)",
  "Pericytes 1",
  "Pro-inflammatory TAMs",
  "Anti-inflammatory TAMs",
  "Angiogenic TAMs",
  "Resident-TAMs"
)
"#####################################################################"
"#####################################################################"
"#####################################################################"
"#####################################################################"
"#####################################################################"
"#####################################################################"


# h5ad 파일명 -> file_id 매핑
h5ad_file_id <- tools::file_path_sans_ext(basename(h5ad_files))

# file_id가 실제로 모두 대응되는지 확인
stopifnot(all(rep_id %in% h5ad_file_id))

for (rid in rep_id) {
  
  # --------------------------------------------------
  # 1. spot_df에서 해당 replicate 추출
  # --------------------------------------------------
  d1 <- spot_df %>%
    dplyr::filter(file_id == rid)
  
  # --------------------------------------------------
  # 2. 기존 spatial plots
  # --------------------------------------------------
  panels <- list(
    plot_spatial(d1, "g_SLC7A7", "SLC7A7")
  )
  
  # --------------------------------------------------
  # 3. 해당 file_id의 h5ad 찾기
  # --------------------------------------------------
  h5ad_idx <- match(rid, h5ad_file_id)
  
  if (is.na(h5ad_idx)) {
    stop("h5ad 파일을 찾을 수 없습니다: ", rid)
  }
  
  # --------------------------------------------------
  # 4. h5ad -> Seurat spatial object
  # --------------------------------------------------
  srt <- schard::h5ad2seurat_spatial(
    h5ad_files[[h5ad_idx]]
  )
  
  # Gene Expression만 유지
  ft <- srt[["Spatial"]]@meta.data$feature_types
  
  srt <- subset(
    srt,
    features = rownames(srt)[ft == "Gene Expression"]
  )
  
  # Normalize
  srt <- NormalizeData(
    srt,
    assay = DefaultAssay(srt)
  )
  
  # SLC7A7 존재 확인
  if (!"SLC7A7" %in% rownames(srt)) {
    stop("SLC7A7이 h5ad 파일에 없습니다: ", rid)
  }
  
  # --------------------------------------------------
  # 5. Seurat SpatialFeaturePlot
  # --------------------------------------------------
  p_srt <- Seurat::SpatialFeaturePlot(
    srt,
    features = "SLC7A7"
  ) +
    ggplot2::ggtitle("SLC7A7")
  
  
  panels <- c(panels, purrr::map(show_states, function(state) {
        
        col <- CS_KEY$col[CS_KEY$label == state]
        
        # state에 대응하는 column이 없으면 명확하게 에러
        if (length(col) != 1 || is.na(col)) {
          stop(
            "CS_KEY에서 state를 찾을 수 없습니다: ",
            state)
        }
        plot_spatial(d1, col, state)})
  )
  
  
  # --------------------------------------------------
  # 6. Dominant niche map
  # --------------------------------------------------
  p_niche_map <- ggplot(
    d1,
    aes(spatial_x, spatial_y, colour = dom_niche)
  ) +
    geom_point(size = .8) +
    coord_fixed() +
    scale_y_reverse() +
    labs(
      title = paste("Dominant niche:", rid),
      colour = NULL
    ) +
    theme_void(base_size = 10)
  
  # --------------------------------------------------
  # 7. 전체 패널 합치기
  # --------------------------------------------------
  p_map <- patchwork::wrap_plots(
    c(panels, list(p_srt, p_niche_map)),
    ncol = 4)
  # --------------------------------------------------
  # 8. 저장
  # --------------------------------------------------
  ggsave(
    file.path(
      FIG_DIR,
      paste0("F6_spatialmap_", rid, ".pdf")
    ),
    p_map,
    width = 15,
    height = 8
  )
}
## ============================================================================
## 15. 분석 5 — niche 단위 pseudobulk DE (기존 파이프라인과 동일한 틀)
## ============================================================================
# spot-level 상관이 약하게 나올 경우의 대안이자 보강.
# (section × dominant niche) 로 raw count를 합쳐서 DESeq2로 넘긴다.
#
# make_pseudobulk <- function(rds_paths, group_col = "dom_niche") {
#   mats <- list(); metas <- list()
#   for (p in rds_paths) {
#     sce <- readRDS(p)
#     d   <- build_spot_df(sce)
#     grp <- paste(d$file_id, d[[group_col]], sep = "__")
#     m   <- counts(sce)[, match(d$spot_id, colnames(sce))]
#     agg <- t(rowsum(t(as.matrix(m)), group = grp))     # 유전자 x 그룹
#     mats[[p]]  <- agg
#     metas[[p]] <- distinct(tibble(sample = grp, file_id = d$file_id,
#                                   tumor_id = d$tumor_id, group = d[[group_col]]))
#     rm(sce, m, agg); gc(verbose = FALSE)
#   }
#   genes <- Reduce(intersect, lapply(mats, rownames))
#   list(counts = round(do.call(cbind, lapply(mats, function(x) x[genes, ]))),
#        meta   = bind_rows(metas))
# }
#
# pb <- make_pseudobulk(rds_files)
# dds <- DESeq2::DESeqDataSetFromMatrix(pb$counts, pb$meta, ~ tumor_id + group)


## ============================================================================
## 16. 결과 한 줄 요약
## ============================================================================
cat("\n=== TAM substate, TAM 총량 보정 후 (tumor-level) ===\n")
print(meta_cor %>% filter(is_TAM) %>% select(label, mean_rho, n_pos, n_tumor, p, fdr))
cat("\n=== SLC7A7 high vs low spot의 TAM 조성 변화 ===\n")
print(delta %>% select(label, mean_diff, n_up, n_sec, p, fdr))



# ==============================================================================
# 03_spatial_revision.R
#   A. 컬럼 세트 감사 (section 10 / 11 불일치 해소)
#   B. 조성 분석 재계산 — TAM 밀도 매칭 + CLR + tumor 단위 검정
#   C. dropout 대응: 이웃 평활 + 검출 여부(binary) 기반 재분석
#   D. 본론: 공간 pseudobulk DE + Hallmark GSEA  <-- 여기가 실제 figure 후보
# ==============================================================================
source(here::here("R", "setup.R"))
BiocManager::install("apeglm")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(purrr); library(ggplot2)
  library(Matrix); library(SingleCellExperiment)
})

## ============================================================================
## A. 컬럼 세트 감사 — 어느 분석이 어떤 state 집합을 썼는지 확정
## ============================================================================
cs_present <- intersect(CS_KEY$col, names(spot_df))

n_tumor_present <- spot_df %>%
  group_by(tumor_id) %>%
  summarise(across(all_of(cs_present), ~ !all(is.na(.x))), .groups = "drop") %>%
  summarise(across(all_of(cs_present), sum)) %>%
  unlist()

audit <- spot_df %>%
  summarise(across(all_of(cs_present), ~ mean(is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "col", values_to = "frac_na") %>%
  left_join(CS_KEY, by = "col") %>%
  mutate(n_tumor_present = n_tumor_present[col]) %>%
  arrange(desc(frac_na))
print(audit %>% filter(label %in% TAM_STATES), n = 20)

# 전체 교집합을 강요하지 말고, state별 n_tumor를 가변으로 두고 보고한다.
# 단 TAM 조성(합=1) 계산에는 '모든 tumor 공통' 집합이 필요하므로 두 세트를 분리.
TAM_ALL      <- CS_KEY$col[CS_KEY$label %in% TAM_STATES]
TAM_ALL      <- intersect(TAM_ALL, cs_present)
TAM_COMMON   <- audit$col[audit$col %in% TAM_ALL & audit$frac_na == 0]
cat("TAM state: 전체", length(TAM_ALL), "/ 모든 section 공통", length(TAM_COMMON), "\n")
cat("공통에서 빠진 것:",
    paste(CS_KEY$label[match(setdiff(TAM_ALL, TAM_COMMON), CS_KEY$col)], collapse=", "), "\n")

# 조성 분석은 TAM_COMMON 으로 통일하고 tam_total 재계산
spot_df$tam_total <- rowSums(spot_df[, TAM_COMMON, drop = FALSE])
spot_df$tam_frac  <- spot_df$tam_total / pmax(spot_df$cell_total, 1e-6)


## ============================================================================
## B. 조성 재분석 — TAM 밀도로 매칭 + CLR + tumor 단위
## ============================================================================
clr <- function(p, eps = 1e-4) { p <- p + eps; log(p) - mean(log(p)) }

# section 내에서 TAM 밀도 10분위로 층화 → 층 안에서만 SLC7A7 high/low 비교.
# 이러면 "SLC7A7 높은 spot = TAM 많은 spot" 교란이 제거된다.
matched_comp <- function(d, tam_cols, n_bin = 10, min_n = 40) {
  empty_out <- tibble(col = character(), d_raw = double(), d_clr = double())
  
  d <- d %>% filter(tam_total > 0) %>% mutate(bin = ntile(tam_total, n_bin))
  if (nrow(d) == 0) return(empty_out)
  
  res <- map_dfr(unique(d$bin), function(b) {
    s <- d %>% filter(bin == b)
    if (nrow(s) < min_n || sd(s$g_SLC7A7) == 0) return(NULL)
    s <- s %>% mutate(grp = ifelse(g_SLC7A7 > median(g_SLC7A7), "High", "Low"))
    if (n_distinct(s$grp) < 2) return(NULL)
    P <- as.matrix(s[, tam_cols, drop = FALSE]) / s$tam_total
    Pc <- t(apply(P, 1, clr))
    hi <- s$grp == "High"
    tibble(col   = tam_cols,
           d_raw = colMeans(P[hi, , drop=FALSE])  - colMeans(P[!hi, , drop=FALSE]),
           d_clr = colMeans(Pc[hi, , drop=FALSE]) - colMeans(Pc[!hi, , drop=FALSE]),
           w     = nrow(s))
  })
  
  if (nrow(res) == 0) return(empty_out)   # <- 핵심: 모든 bin이 스킵된 경우
  
  res %>%
    group_by(col) %>%
    summarise(d_raw = weighted.mean(d_raw, w), d_clr = weighted.mean(d_clr, w), .groups = "drop")
}

comp_sec2 <- spot_df %>%
  group_by(tumor_id, file_id) %>%
  group_modify(~ matched_comp(.x, TAM_COMMON)) %>%
  ungroup()

# section -> tumor -> n=12 검정
comp_tumor <- comp_sec2 %>%
  group_by(tumor_id, col) %>%
  summarise(across(c(d_raw, d_clr), mean), .groups = "drop")

comp_res <- comp_tumor %>%
  group_by(col) %>%
  summarise(mean_raw = mean(d_raw), mean_clr = mean(d_clr),
            sd_clr = sd(d_clr), n_tumor = n(), n_up = sum(d_clr > 0),
            p = tryCatch(wilcox.test(d_clr, mu = 0)$p.value, error=function(e) NA),
            .groups = "drop") %>%
  mutate(fdr = p.adjust(p, "BH")) %>%
  left_join(CS_KEY, by = "col") %>%
  arrange(desc(mean_clr))

print(comp_res %>% select(label, mean_raw, mean_clr, n_up, n_tumor, p, fdr), n = 20)
write.csv(comp_res, file.path(RES_DIR, "TAMcomposition_matched_tumorlevel.csv"), row.names = FALSE)

# 이전(section 단위, 비매칭) 결과와 나란히 비교 — 얼마나 달라지는지 보고용
p_compare <- comp_res %>%
  select(label, matched_tumor = mean_raw) %>%
  left_join(read.csv(file.path(RES_DIR, "TAMcomposition_SLC7A7high_vs_low.csv")) %>%
              select(label, naive_section = mean_diff), by = "label") %>%
  pivot_longer(-label) %>%
  ggplot(aes(value, reorder(label, value), colour = name)) +
  geom_vline(xintercept = 0, linetype = 2) + geom_point(size = 2.5) +
  labs(x = "Δ fraction (SLC7A7 high − low)", y = NULL, colour = NULL,
       title = "밀도 매칭 + tumor 단위 집계 전후 비교")
ggsave(file.path(FIG_DIR, "S1_composition_naive_vs_matched.pdf"), p_compare, width = 8, height = 4)


## ============================================================================
## C. dropout 대응 — 이웃 평활 & 검출 여부
## ============================================================================
# install.packages("FNN")
library(FNN)

smooth_knn <- function(d, col = "g_SLC7A7", k = 6) {
  xy <- as.matrix(d[, c("spatial_x", "spatial_y")])
  nn <- FNN::get.knn(xy, k = k)$nn.index
  v  <- d[[col]]
  rowMeans(cbind(v, matrix(v[nn], nrow = nrow(d))))   # 자기 자신 + 6-이웃
}

spot_df <- spot_df %>%
  group_by(file_id) %>%
  mutate(slc_smooth = smooth_knn(cur_data_all()),
         slc_det    = as.integer(g_SLC7A7 > 0)) %>%
  ungroup()

# 검출률 (감쇠가 얼마나 심한지의 지표)
spot_df %>% group_by(tumor_id) %>%
  summarise(det = mean(slc_det), .groups="drop") %>% print(n = 20)

# 평활값으로 partial correlation 재계산 → 감쇠가 원인이었는지 판정
cor_smooth <- spot_df %>%
  group_by(tumor_id, file_id) %>%
  group_modify(~ section_pcor(.x, intersect(BONA_COLS, names(.x)),
                              target = "slc_smooth", ctrl = "tam_total")) %>%
  ungroup() %>% left_join(CS_KEY, by = "col") %>%
  mutate(is_TAM = label %in% TAM_STATES)

cor_smooth %>% filter(is_TAM) %>%
  group_by(label, tumor_id) %>% summarise(r = mean(rho_adj), .groups="drop") %>%
  group_by(label) %>%
  summarise(mean_rho = mean(r), n_pos = sum(r > 0), n = n(),
            p = tryCatch(wilcox.test(r, mu=0)$p.value, error=function(e) NA)) %>%
  arrange(desc(mean_rho)) %>% print(n = 20)
# 평활 후 |rho|가 눈에 띄게 커지면 -> 원래 결과는 dropout 감쇠 탓
# 그대로면 -> 진짜 null. 이 판정을 논문에 명시할 것.


## ============================================================================
## D. 공간 pseudobulk DE + GSEA  ← 실제로 figure가 될 부분
## ============================================================================
# 설계: section 안에서 TAM 밀도 층화 후 SLC7A7 high/low 두 그룹으로 나누고
#       raw count를 합쳐 (section × group) pseudobulk 생성.
#       ~ file_id + group  (section을 blocking → GBmap의 donor blocking과 동형)

assign_groups <- function(d, n_bin = 5, min_n = 50) {
  d %>%
    filter(tam_total > 0) %>%
    mutate(bin = ntile(tam_total, n_bin)) %>%
    group_by(bin) %>%
    filter(n() >= min_n, stats::sd(slc_smooth) > 0) %>%
    mutate(group = case_when(
      slc_smooth >= quantile(slc_smooth, 0.70) ~ "high",
      slc_smooth <= quantile(slc_smooth, 0.30) ~ "low",
      TRUE ~ NA_character_
    )) %>%
    ungroup() %>%
    filter(!is.na(group)) %>%
    select(spot_id, tumor_id, group)
}

spot_groups <- spot_df %>%
  group_by(file_id) %>%
  group_modify(~ assign_groups(.x)) %>%
  ungroup()

spot_groups %>% dplyr::count(group)

make_spatial_pb <- function(rds_paths, groups) {
  mats <- list(); metas <- list()
  for (p in rds_paths) {
    sce <- readRDS(p)
    fid <- as.character(colData(sce)$file_id[1])
    g   <- groups %>% filter(file_id == fid)
    if (nrow(g) < 100) { rm(sce); gc(verbose=FALSE); next }
    
    idx  <- match(g$spot_id, colnames(sce)); keep <- !is.na(idx)
    lab  <- factor(paste(fid, g$group[keep], sep = "__"))
    
    # 추가: 매칭 후 group이 1종류뿐이면 pseudobulk 비교 불가 -> skip
    if (nlevels(lab) < 2) {
      cat("x", sep = "")   # 스킵된 파일 표시(선택)
      rm(sce); gc(verbose = FALSE)
      next
    }
    
    m    <- counts(sce)[, idx[keep], drop = FALSE]
    D    <- Matrix::sparse.model.matrix(~ 0 + lab)
    agg  <- as.matrix(m %*% D); colnames(agg) <- levels(lab)
    mats[[fid]]  <- agg
    metas[[fid]] <- tibble(sample = colnames(agg), file_id = fid,
                           tumor_id = g$tumor_id[1],
                           group = sub(".*__", "", colnames(agg)))
    rm(sce, m, D, agg); gc(verbose = FALSE)
    cat(".", sep = "")
  }
  genes <- Reduce(intersect, lapply(mats, rownames))
  list(counts = round(do.call(cbind, lapply(mats, function(x) x[genes, , drop=FALSE]))),
       meta   = bind_rows(metas))
}

pb <- make_spatial_pb(rds_files, spot_groups)
dim(pb$counts); table(pb$meta$group)
saveRDS(pb, file.path(RES_DIR, "spatial_pseudobulk.rds"))
pb <- readRDS(file.path(RES_DIR, "spatial_pseudobulk.rds"))

# --- DESeq2 -----------------------------------------------------------------
suppressPackageStartupMessages(library(DESeq2))
meta <- pb$meta %>% column_to_rownames("sample")
meta <- meta[colnames(pb$counts), ]
meta$group   <- factor(meta$group, c("low", "high"))
meta$file_id <- factor(meta$file_id)

# section 짝이 완전한 것만 (high/low 둘 다 있는 section)
ok_sec <- names(which(table(meta$file_id) == 2))
sel    <- meta$file_id %in% ok_sec
dds <- DESeqDataSetFromMatrix(pb$counts[, sel], meta[sel, ], ~ file_id + group)
dds <- dds[rowSums(counts(dds) >= 10) >= 0.2 * ncol(dds), ]
dds <- DESeq(dds)
## ---- 1. 통계량(순위용) + shrunken LFC(보고용) 을 각각 받아서 합친다 --------
res_stat <- results(dds, name = "group_high_vs_low")          # stat 있음
res_shr  <- lfcShrink(dds, coef = "group_high_vs_low", type = "apeglm")

res <- as.data.frame(res_stat) %>%
  rownames_to_column("gene") %>%
  select(gene, baseMean, lfc_raw = log2FoldChange, stat, pvalue, padj) %>%
  left_join(as.data.frame(res_shr) %>% rownames_to_column("gene") %>%
              select(gene, lfc_shrunk = log2FoldChange, lfcSE_shrunk = lfcSE),
            by = "gene") %>%
  arrange(padj)

write.csv(res, file.path(RES_DIR, "spatial_pseudobulk_DE.csv"), row.names = FALSE)

## ---- 2. 순위 벡터 (data mask 사고 방지를 위해 $ 로 명시적 접근) ------------
 


# --- Hallmark GSEA: GBmap 결과가 공간에서도 재현되는가 -----------------------
suppressPackageStartupMessages({ library(fgsea); library(msigdbr) })
hall <- msigdbr(species = "Homo sapiens", collection = "H")
paths <- split(hall$gene_symbol, hall$gs_name)

SEED            <- 42
INTER_GENE_COR  <- 0.01 
HALLMARK <- get_hallmark_list()
GOBP <- get_gobp_list()
GOMF <- get_gomf_list()
GOCC <- get_gocc_list()

run_gsea_pair <- function(stat_vec, term) {
  stat_vec <- sort(stat_vec[is.finite(stat_vec)], decreasing = TRUE)
  set.seed(SEED)
  fg <- fgsea::fgsea(pathways = term, stats = stat_vec,
                     minSize = 10, maxSize = 500, eps = 0, nPermSimple = 100000)
  idx <- limma::ids2indices(term, names(stat_vec))
  idx <- idx[vapply(idx, length, 1L) >= 10]
  cam <- limma::cameraPR(stat_vec, idx, inter.gene.cor = INTER_GENE_COR,
                         use.ranks = FALSE)
  cam$pathway <- rownames(cam)
  out <- fg %>%
    dplyr::select(pathway, NES, pval, padj, size) %>%
    dplyr::left_join(cam %>% dplyr::select(pathway, Direction,
                                           camera_p = PValue, camera_FDR = FDR),
                     by = "pathway") %>%
    dplyr::arrange(pval)
  out
}
stat <- res %>% filter(!is.na(stat)) %>% { setNames(.$stat, .$gene) }

gs_hall <- run_gsea_pair(stat, HALLMARK) %>% arrange(padj)
gs_bp <- run_gsea_pair(stat, GOBP) %>% arrange(padj)
gs_mf <- run_gsea_pair(stat, GOMF) %>% arrange(padj)
gs_CC <- run_gsea_pair(stat, GOCC) %>% arrange(padj)

write_result(gs_hall, "Fig8_spatial_gsea_HALL.csv")
write_result(gs_bp, "Fig8_spatial_gsea_BP.csv")
write_result(gs_mf, "Fig8_spatial_gsea_MF.csv")
write_result(gs_CC, "Fig8_spatial_gsea_CC.csv")

# GBmap에서 나온 3개 경로만 직접 확인
key <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB",
         "HALLMARK_APOPTOSIS",
         "HALLMARK_INTERFERON_ALPHA_RESPONSE",
         "HALLMARK_INTERFERON_GAMMA_RESPONSE",
         "HALLMARK_HYPOXIA", "HALLMARK_GLYCOLYSIS")
print(gs %>% filter(pathway %in% key) %>% select(pathway, NES, pval, padj))

p_gsea <- gs %>% filter(padj < 0.25) %>%
  mutate(pathway = gsub("HALLMARK_", "", pathway),
         key = pathway %in% gsub("HALLMARK_", "", key)) %>%
  ggplot(aes(NES, reorder(pathway, NES), fill = key)) +
  geom_col() + scale_fill_manual(values = c("grey70", "firebrick"), guide = "none") +
  geom_vline(xintercept = 0) +
  labs(x = "NES (SLC7A7-high vs low spots)", y = NULL,
       title = "Spatial pseudobulk GSEA — Hallmark")
ggsave(file.path(FIG_DIR, "F7_spatial_GSEA.pdf"), p_gsea, width = 8, height = 6)


"======================================================"


  # 10_spatial.R 맨 아래에 이어붙일 것.
  #
  # 목적: 공간 pseudobulk에서 나온 IFN 상향이 단순히 "SLC7A7-high spot에
  #       면역세포가 더 많아서" 생긴 조성 효과가 아님을 보이는 통제 분석.
  #       논문에서는 이 표/그림 하나만 보고하면 충분하다.
  # ==============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(purrr); library(ggplot2)
})

## ---- 0. res 에 stat 컬럼 보장 (apeglm 결과만 있는 경우 대비) ---------------
if (!"stat" %in% names(res)) {
  res <- as.data.frame(results(dds, name = "group_high_vs_low")) %>%
    rownames_to_column("gene") %>%
    select(gene, baseMean, stat, pvalue, padj) %>%
    left_join(res %>% select(gene, lfc_shrunk = log2FoldChange), by = "gene")
}

## ---- 1. 모듈 정의 ----------------------------------------------------------
# 주의: Mx1은 C57BL/6 계열에서 기능 결실이므로 마우스 실험과 비교할 때는 제외.
#       여기는 인간 조직이라 포함해도 되지만 일관성 위해 뺐다.
MODULES <- list(
  "ISG (interferon)"      = c("IFIT1","IFIT3","ISG15","OAS1","OAS2","OAS3","IFI6",
                              "USP18","HERC6","CMPK2","RSAD2","STAT1","EIF2AK2","XAF1"),
  "Pan-myeloid"           = c("PTPRC","AIF1","CSF1R","C1QA","C1QB","C1QC","CD68",
                              "TYROBP","FCER1G","ITGAM","CD14","LAPTM5","MPEG1"),
  "Homeostatic microglia" = c("P2RY12","TMEM119","CX3CR1","OLFML3","SALL1"),
  "M2 / anti-inflammatory"= c("CD163","MRC1","MSR1","TREM2","STAB1","SPP1","APOE"),
  "NF-kB effectors"       = c("TNF","TNFAIP3","CCL2","CXCL2","ICAM1","RELB",
                              "NFKBIA","IL1B","CD83"),
  "Hypoxia"               = c("SLC2A1","NDRG1","ADM","ANGPTL4","CA9","VEGFA",
                              "HK2","P4HA1","ENO1")
)

mod_df <- imap_dfr(MODULES, function(g, nm) {
  res %>% filter(gene %in% g, !is.na(stat)) %>% mutate(module = nm)
}) %>% mutate(module = factor(module, names(MODULES)))

bg <- res$stat[!is.na(res$stat)]        # 전체 유전자 배경분포

## ---- 2. 모듈 요약 + 배경 대비 검정 -----------------------------------------
mod_summary <- mod_df %>%
  group_by(module) %>%
  summarise(n_detected  = n(),
            mean_stat   = mean(stat),
            median_stat = median(stat),
            mean_lfc    = mean(lfc_shrunk, na.rm = TRUE),
            n_padj05    = sum(padj < 0.05, na.rm = TRUE),
            p_vs_bg     = wilcox.test(stat, bg)$p.value,
            .groups = "drop") %>%
  mutate(fdr_vs_bg = p.adjust(p_vs_bg, "BH")) %>%
  arrange(desc(mean_stat))

print(as.data.frame(mod_summary), digits = 3)
cat("\n배경(전체 유전자) mean stat =", round(mean(bg), 3), "\n\n")

## ---- 3. 핵심 검정: ISG vs 범-myeloid 직접 비교 ------------------------------
isg <- mod_df$stat[mod_df$module == "ISG (interferon)"]
mye <- mod_df$stat[mod_df$module == "Pan-myeloid"]

wt <- wilcox.test(isg, mye)
tt <- t.test(isg, mye)
cat("ISG vs Pan-myeloid\n",
    " mean stat : ", round(mean(isg), 2), " vs ", round(mean(mye), 2), "\n",
    " padj<0.05 : ", sum(mod_df$padj[mod_df$module=="ISG (interferon)"] < 0.05, na.rm=TRUE),
    "/", length(isg), " vs ",
    sum(mod_df$padj[mod_df$module=="Pan-myeloid"] < 0.05, na.rm=TRUE), "/", length(mye), "\n",
    " Wilcoxon p: ", format.pval(wt$p.value, digits = 3), "\n",
    " Welch t  p: ", format.pval(tt$p.value, digits = 3), "\n\n", sep = "")

## ---- 4. cameraPR — 유전자 간 상관을 보정한 competitive test -----------------
# reviewer가 "gene set 안 유전자들은 서로 상관이 있어 p값이 부풀려진다"고
# 지적할 수 있으므로 같이 보고하면 방어가 쉽다.
if (requireNamespace("limma", quietly = TRUE)) {
  st  <- setNames(res$stat, res$gene); st <- st[!is.na(st)]
  idx <- limma::ids2indices(MODULES, names(st))
  cam <- limma::cameraPR(st, idx) %>% rownames_to_column("module")
  print(cam, digits = 3)
  mod_summary <- left_join(mod_summary,
                           cam %>% select(module, cam_dir = Direction,
                                          cam_p = PValue, cam_fdr = FDR),
                           by = "module")
}

write.csv(mod_summary, file.path(RES_DIR, "spatial_module_control.csv"), row.names = FALSE)

## ---- 5. Figure (supplementary) ---------------------------------------------
p_mod <- ggplot(mod_df, aes(stat, module)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_boxplot(outlier.shape = NA, fill = "grey93", width = .6) +
  geom_point(aes(colour = padj < 0.05),
             position = position_jitter(height = .18, seed = 1),
             size = 1.9, alpha = .9) +
  scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick"),
                      name = "padj < 0.05", na.translate = FALSE) +
  scale_y_discrete(limits = rev(levels(mod_df$module))) +
  labs(x = "Wald statistic (SLC7A7-high vs low spots)", y = NULL,
       title = "Marker module response in spatial pseudobulk",
       subtitle = "IFN 상향은 범-myeloid 조성 증가로 설명되지 않음")

ggsave(file.path(FIG_DIR, "S2_module_control.pdf"), p_mod, width = 8, height = 4.8)
p_mod

