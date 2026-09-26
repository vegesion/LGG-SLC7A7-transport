# =============================================================================
# R/scrna.R  ---  단일세포(GBmap) 분석 함수 / single-cell helpers
# -----------------------------------------------------------------------------
# 원본 GBmap 분석.R / patient_pseudo_bulk.R / cellchat.R 의 로직을 함수화.
# 대상: GBmap core atlas의 myeloid(특히 microglia)에서 SLC7A7 +/- 그룹 비교.
# =============================================================================

# ── h5ad 를 R에서 직접 희소행렬로 읽기 / read h5ad as sparse WITHOUT Python ───
# reticulate/anndata 버전에 따라 sparse X 가 dense 로 변환되며 RAM 이 폭발하는 문제를
# 원천 차단합니다(파이썬 미경유). hdf5r + Matrix 로 CSR/CSC 를 그대로 dgCMatrix(gene×cell)
# 로 구성합니다. obs(메타데이터)의 categorical 인코딩(categories/codes)도 처리.
# ── h5ad 크기만 먼저 확인 (데이터 미읽음, 메모리 거의 0) / peek h5ad size ─────
# X 를 읽지 않고 shape 와 nnz(비영요소 수)만 조회해 예상 메모리를 알려줍니다.
# 반드시 read_h5ad_sparse() 전에 이걸로 규모를 먼저 확인하세요.
h5ad_info <- function(path) {
  stopifnot(requireNamespace("hdf5r", quietly = TRUE))
  h5 <- hdf5r::H5File$new(path, mode = "r"); on.exit(h5$close_all(), add = TRUE)
  Xg    <- h5[["X"]]
  enc   <- hdf5r::h5attr(Xg, "encoding-type")
  shape <- hdf5r::h5attr(Xg, "shape")                 # c(n_obs, n_var)
  nnz   <- Xg[["data"]]$dims                          # 길이만(읽지 않음)
  est_gb <- nnz * (8 + 4) / 1e9                        # dgCMatrix: x(double)+i(int)
  message(sprintf(
    "h5ad: %s\n  cells = %s, genes = %s, encoding = %s\n  nnz = %s  →  예상 sparse 메모리 ~ %.1f GB",
    basename(path), format(shape[1], big.mark = ","), format(shape[2], big.mark = ","),
    enc, format(nnz, big.mark = ","), est_gb))
  invisible(list(cells = shape[1], genes = shape[2], nnz = nnz, est_gb = est_gb, encoding = enc))
}

read_h5ad_sparse <- function(path) {
  stopifnot(requireNamespace("hdf5r", quietly = TRUE),
            requireNamespace("Matrix", quietly = TRUE))
  h5 <- hdf5r::H5File$new(path, mode = "r")
  on.exit(h5$close_all(), add = TRUE)

  # --- X (sparse) → genes(var) × cells(obs) dgCMatrix ---
  Xg    <- h5[["X"]]
  enc   <- hdf5r::h5attr(Xg, "encoding-type")            # "csr_matrix" / "csc_matrix"
  shape <- hdf5r::h5attr(Xg, "shape")                    # c(n_obs, n_var)
  dat   <- Xg[["data"]]$read()
  idx   <- Xg[["indices"]]$read()
  ptr   <- Xg[["indptr"]]$read()
  n_obs <- shape[1]; n_var <- shape[2]
  if (identical(enc, "csr_matrix")) {
    # CSR(obs-major): p=ptr(길이 n_obs+1), i=var index → CSC 로 dims=(n_var,n_obs)=gene×cell
    counts <- Matrix::sparseMatrix(i = idx + 1L, p = ptr, x = dat,
                                   dims = c(n_var, n_obs), index1 = TRUE)
  } else {
    # CSC(var-major): obs×var 로 만든 뒤 전치
    X <- Matrix::sparseMatrix(i = idx + 1L, p = ptr, x = dat,
                              dims = c(n_obs, n_var), index1 = TRUE)
    counts <- Matrix::t(X)
  }

  # --- 이름 (var/_index = genes, obs/_index = cells) ---
  read_index <- function(grp) {
    ix <- hdf5r::h5attr(h5[[grp]], "_index")
    h5[[paste0(grp, "/", ix)]]$read()
  }
  rownames(counts) <- read_index("var")
  colnames(counts) <- read_index("obs")

  # --- obs 메타데이터: read_obs_all() 로 일괄 처리(길이 가드 포함) ---
  meta <- read_obs_all(path, cells = colnames(counts))

  # --- UMAP ---
  umap <- NULL
  if ("obsm" %in% names(h5) && "X_umap" %in% names(h5[["obsm"]])) {
    um <- h5[["obsm/X_umap"]]$read()
    if (nrow(um) != n_obs) um <- t(um)               # (2×n_obs) 로 저장된 경우 전치
    rownames(um) <- colnames(counts)
    umap <- um
  }
  list(counts = counts, meta = meta, umap = umap)
}

# ── h5ad → Seurat 객체 / build Seurat with UMAP (메모리 안전) ─────────────────
# use_python = FALSE(기본): hdf5r 직접 읽기(권장, RAM 안전).
# use_python = TRUE       : 예전처럼 anndata::read_h5ad 사용(호환용, 대용량은 위험).
load_gbmap_seurat <- function(h5ad_path = PATH_H5AD, use_python = FALSE) {
  if (use_python) {
    adata  <- anndata::read_h5ad(h5ad_path)
    counts <- Matrix::t(adata$X)
    gene_names <- reticulate::py_to_r(adata$var_names$tolist())
    dimnames(counts) <- list(gene_names, rownames(reticulate::py_to_r(adata$obs)))
    meta <- adata$obs
    umap <- as.matrix(adata$obsm[["X_umap"]])
  } else {
    dat    <- read_h5ad_sparse(h5ad_path)            # 파이썬 미경유, sparse 유지
    counts <- dat$counts; meta <- dat$meta; umap <- dat$umap
  }

  seu <- Seurat::CreateSeuratObject(counts = counts, meta.data = meta)
  if (!is.null(umap)) {
    rownames(umap) <- colnames(seu)
    seu[["umap"]] <- Seurat::CreateDimReducObject(embeddings = umap, key = "UMAP_")
  }
  seu <- Seurat::NormalizeData(seu)

  # ENSEMBL rowname → SYMBOL (매핑 실패 시 원본 ENSEMBL 유지)
  gene_map <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
    keys = rownames(seu), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  gene_map[is.na(gene_map)] <- rownames(seu)[is.na(gene_map)]
  rownames(seu@assays$RNA@layers$counts) <- gene_map
  rownames(seu@assays$RNA) <- make.unique(gene_map)
  seu
}

# ── obs 메타데이터 전체 읽기 / read ALL obs columns (길이 가드 포함) ──────────
# 범주형(categories/codes), nullable(values/mask, AnnData 0.10), 일반 배열 모두 처리.
# cells 인자를 주면 그 순서/길이에 맞추고, 세포 수와 길이가 다른 컬럼은 건너뜁니다.
read_obs_all <- function(path, cells = NULL) {
  stopifnot(requireNamespace("hdf5r", quietly = TRUE))
  h5  <- hdf5r::H5File$new(path, "r"); on.exit(h5$close_all(), add = TRUE)
  obs <- h5[["obs"]]
  idx <- hdf5r::h5attr(obs, "_index")
  file_cells <- obs[[idx]]$read()
  n <- length(file_cells)
  cols <- setdiff(names(obs), c(idx, "_index", "__categories"))
  out <- data.frame(row.names = file_cells)
  skipped <- character(0)
  for (nm in cols) {
    nd <- obs[[nm]]; v <- NULL
    if (inherits(nd, "H5Group")) {
      ch <- names(nd)
      if (all(c("categories", "codes") %in% ch)) {
        cats  <- nd[["categories"]]$read(); codes <- as.integer(nd[["codes"]]$read())
        v <- factor(cats[codes + 1L], levels = cats)
      } else if (all(c("values", "mask") %in% ch)) {          # nullable
        v <- nd[["values"]]$read(); v[as.logical(nd[["mask"]]$read())] <- NA
      }
    } else if (inherits(nd, "H5D")) {
      v <- tryCatch(nd$read(), error = function(e) NULL)
    }
    if (is.null(v) || length(v) != n) { skipped <- c(skipped, nm); next }
    out[[nm]] <- v
  }
  if (length(skipped))
    message("read_obs_all: 길이 불일치로 건너뜀 (세포 ", n, "): ", paste(skipped, collapse = ", "))
  if (!is.null(cells)) out <- out[cells, , drop = FALSE]        # seu 세포 순서에 맞추기
  out
}

# ── BPCells 온디스크 행렬 + h5ad 메타로 Seurat v5 구성 (RAM 안전) ─────────────
# bp_dir: write_matrix_dir 로 만든 폴더. 없으면 h5ad 에서 만들어 둔다.
# load_gbmap_bpcells <- function(bp_dir, h5ad_path = PATH_H5AD, normalize = TRUE) {
#   stopifnot(requireNamespace("BPCells", quietly = TRUE))
#   if (!dir.exists(bp_dir)) {
#     message("[BPCells] 온디스크 행렬 생성(최초 1회) → ", bp_dir)
#     mat <- BPCells::open_matrix_anndata_hdf5(h5ad_path)
#     BPCells::write_matrix_dir(mat, bp_dir)
#   }
#   counts <- BPCells::open_matrix_dir(bp_dir)
#   if (nrow(counts) > ncol(counts)) counts <- t(counts)          # genes×cells (전치는 공짜)
# 
#   if (any(grepl("^ENSG", utils::head(rownames(counts), 50)))) { # ENSEMBL → SYMBOL
#     g   <- rownames(counts)
#     sym <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = g,
#              column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
#     sym[is.na(sym)] <- g[is.na(sym)]
#     rownames(counts) <- make.unique(sym)
#   }
# 
#   seu  <- Seurat::CreateSeuratObject(counts = counts)
#   meta <- read_obs_all(h5ad_path, cells = colnames(seu))        # obs 전체
#   seu  <- Seurat::AddMetaData(seu, meta)
# 
#   h5 <- hdf5r::H5File$new(h5ad_path, "r")                       # UMAP
#   if ("obsm" %in% names(h5) && "X_umap" %in% names(h5[["obsm"]])) {
#     um <- h5[["obsm/X_umap"]]$read()
#     if (nrow(um) != ncol(seu)) um <- t(um)
#     rownames(um) <- colnames(seu); colnames(um) <- c("UMAP_1", "UMAP_2")
#     seu[["umap"]] <- Seurat::CreateDimReducObject(embeddings = um, key = "UMAP_", assay = "RNA")
#   }
#   h5$close_all()
# 
#   if (normalize) {seu <- Seurat::NormalizeData(seu); seu <- FindVariableFeatures(seu); seu <- ScaleData(seu); seu <- RunPCA(seu, npcs = 50, verbose = FALSE)}      # CellChat/FetchData 에 필요
#   seu
# }


h5_read_vec <- function(h5, path) {
  if (!h5$exists(path)) return(NULL)
  o <- h5[[path]]
  if (inherits(o, "H5Group")) {
    cats  <- as.character(o[["categories"]]$read())
    codes <- o[["codes"]]$read()
    out <- cats[codes + 1L]; out[codes < 0] <- NA_character_
    out
  } else as.vector(o$read())
}

REQUIRED_OBS <- c("donor_id", "author", "suspension_type", "assay",
                  "annotation_level_3", "annotation_level_4")

load_gbmap_bpcells <- function(bp_root, h5ad_path = PATH_H5AD,
                               group      = "raw/X",
                               normalize  = TRUE,
                               run_pca    = FALSE,
                               as_uint32  = TRUE,
                               reductions = c("X_umap", "X_scANVI", "X_scanvi_emb",
                                              "X_pca", "X_scvi")) {
  stopifnot(requireNamespace("BPCells", quietly = TRUE),
            requireNamespace("hdf5r",   quietly = TRUE))
  
  ## === phase 1 : h5 에서 이름 정보만 읽고 즉시 닫기 ==========================
  info <- local({
    h5 <- hdf5r::H5File$new(h5ad_path, "r")
    on.exit(h5$close_all(), add = TRUE)
    if (!h5$exists(group)) stop("h5ad 에 '", group, "' 없음.")
    vp <- if (startsWith(group, "raw/") && h5$exists("raw/var/_index")) "raw/var" else "var"
    list(
      var_prefix  = vp,
      gene_id     = h5_read_vec(h5, paste0(vp, "/_index")),
      gene_name   = h5_read_vec(h5, paste0(vp, "/feature_name")),
      cell_id     = h5_read_vec(h5, "obs/_index"),
      n_filtered  = if (h5$exists("var/feature_is_filtered"))
        sum(h5[["var/feature_is_filtered"]]$read()) else NA_integer_,
      obsm_names  = if ("obsm" %in% names(h5)) names(h5[["obsm"]]) else character(0)
    )
  })
  n_var <- length(info$gene_id); n_obs <- length(info$cell_id)
  message(sprintf("[h5ad] var = %s | n_var = %d | n_obs = %d | feature_is_filtered = %s",
                  info$var_prefix, n_var, n_obs, info$n_filtered))
  message("[obsm] ", paste(info$obsm_names, collapse = ", "))
  
  ## === phase 2 : BPCells 온디스크 행렬 =======================================
  bp_dir <- file.path(bp_root, gsub("[^A-Za-z0-9]", "_", group))
  if (!dir.exists(bp_dir)) {
    message("[BPCells] 생성 (", group, ") → ", bp_dir, "  ※ 수십 분")
    dir.create(bp_root, recursive = TRUE, showWarnings = FALSE)
    mat <- BPCells::open_matrix_anndata_hdf5(h5ad_path, group = group)
    if (as_uint32) mat <- BPCells::convert_matrix_type(mat, "uint32_t")
    BPCells::write_matrix_dir(mat, bp_dir)
  } else {
    message("[BPCells] 기존 캐시 사용: ", bp_dir)
  }
  counts <- BPCells::open_matrix_dir(bp_dir)
  
  d <- as.integer(dim(counts))
  if (identical(d, c(n_obs, n_var))) {
    counts <- t(counts)
  } else if (!identical(d, c(n_var, n_obs))) {
    stop(sprintf("행렬 %d x %d 가 n_var=%d / n_obs=%d 와 불일치", d[1], d[2], n_var, n_obs))
  }
  
  sym <- if (!is.null(info$gene_name) && length(info$gene_name) == n_var) {
    ifelse(is.na(info$gene_name) | info$gene_name == "", info$gene_id, info$gene_name)
  } else { warning("feature_name 없음 → ENSEMBL ID 유지"); info$gene_id }
  
  dup <- unique(sym[duplicated(sym)])
  if (length(dup)) message("[genes] 중복 심볼: ", paste(dup, collapse = ", "))
  rownames(counts) <- make.unique(as.character(sym))
  colnames(counts) <- as.character(info$cell_id)
  
  seu <- Seurat::CreateSeuratObject(counts = counts)
  
  # Seurat 이 밑줄을 대시로 치환하므로 실제 바뀐 이름 확인
  changed <- setdiff(make.unique(as.character(sym)), rownames(seu))
  if (length(changed))
    message(sprintf("[genes] Seurat 이 이름 변경한 유전자 %d 개 (예: %s)",
                    length(changed), paste(utils::head(changed, 5), collapse = ", ")))
  
  ## === phase 3 : obs 메타데이터 (여기서 h5 핸들 안 들고 있어야 함) ===========
  meta <- read_obs_all(h5ad_path, cells = colnames(seu))
  if (is.null(rownames(meta))) rownames(meta) <- as.character(info$cell_id)
  stopifnot(identical(rownames(meta), colnames(seu)))
  seu <- Seurat::AddMetaData(seu, meta)
  
  miss <- setdiff(REQUIRED_OBS, colnames(seu@meta.data))
  if (length(miss))
    warning("필수 obs 컬럼 누락: ", paste(miss, collapse = ", "),
            " — read_obs_all 이 건너뛰었는지 확인할 것")
  
  ## === phase 4 : obsm 임베딩 (h5 재오픈) =====================================
  emb <- local({
    tgt <- intersect(reductions, info$obsm_names)
    if (!length(tgt)) return(list())
    h5 <- hdf5r::H5File$new(h5ad_path, "r")
    on.exit(h5$close_all(), add = TRUE)
    stats::setNames(lapply(tgt, function(nm) h5[[paste0("obsm/", nm)]]$read()), tgt)
  })
  for (nm in names(emb)) {
    em <- emb[[nm]]
    if (nrow(em) != ncol(seu)) em <- t(em)
    if (nrow(em) != ncol(seu)) next
    red <- gsub("^X_", "", nm); key <- paste0(red, "_")
    rownames(em) <- colnames(seu); colnames(em) <- paste0(key, seq_len(ncol(em)))
    seu[[red]] <- Seurat::CreateDimReducObject(embeddings = em, key = key, assay = "RNA")
  }
  
  ## === phase 5 : 정수 검증 → 정규화 ==========================================
  smp <- SeuratObject::LayerData(seu, layer = "counts")[
    seq_len(min(2000, nrow(seu))), seq_len(min(200, ncol(seu)))]
  smp <- as.matrix(smp)                    # ★ BPCells lazy subset → base matrix
  v <- smp[smp > 0]
  is_raw <- length(v) > 0 && all(abs(v - round(v)) < 1e-8)
  message(sprintf("[counts] 정수 = %s | 예시: %s", is_raw,
                  paste(utils::head(sort(unique(v)), 6), collapse = ", ")))
  if (!is_raw) stop("counts 가 정수 아님 — group 인자 확인 (raw/X).")
  
  if (normalize) seu <- Seurat::NormalizeData(seu, verbose = FALSE)
  if (run_pca) {
    seu <- Seurat::FindVariableFeatures(seu, verbose = FALSE)
    seu <- Seurat::ScaleData(seu, verbose = FALSE)      # ~5 GB dense, 배치 미보정
    seu <- Seurat::RunPCA(seu, npcs = 50, verbose = FALSE)
  }
  
  message(sprintf("[done] %d genes x %d cells | layer = %s", nrow(seu), ncol(seu), group))
  seu
}

# ── SLC7A7 발현 기준 +/- 라벨 부여 / label Positive/Negative by expression ───
label_slc_status <- function(obj, gene = GENE_OF_INTEREST,
                             threshold = SCRNA_POS_THRESHOLD,
                             col = paste0(GENE_OF_INTEREST, "_status")) {
  expr <- Seurat::FetchData(obj, vars = gene)
  obj[[col]] <- ifelse(expr[[gene]] > threshold, "Positive", "Negative")
  obj
}

# ── SLC7A7 +/- DEG (Wilcoxon, 전체 유전자) / FindMarkers pos vs neg ──────────
find_slc_deg <- function(obj, status_col = paste0(GENE_OF_INTEREST, "_status"),
                         all_genes = TRUE) {
  Seurat::Idents(obj) <- status_col
  args <- list(object = obj, ident.1 = "Positive", ident.2 = "Negative",
               test.use = "wilcox", verbose = TRUE)
  if (all_genes) { args$min.pct <- 0; args$logfc.threshold <- 0; args$return.thresh <- 1 }
  else           { args$min.pct <- 0.1; args$logfc.threshold <- 0.25 }
  do.call(Seurat::FindMarkers, args)
}

# ── 환자별 pseudobulk + paired limma-voom / patient pseudobulk DE ────────────
# 원본 patient_pseudo_bulk.R: donor를 covariate로 넣은 paired design, filterByExpr,
# voom+limma, coef = <GENE>_statusPositive.
run_pseudobulk_limma <- function(obj, status_col = paste0(GENE_OF_INTEREST, "_status"),
                                 min_cells = PB_MIN_CELLS, min_count = PB_FILTER_MIN_COUNT,
                                 pos_level = "Positive", neg_level = "Negative") {
  status_vec <- obj[[status_col]][, 1]

  # 1. donor×status 조합 최소 세포 수 필터 (pseudobulk 노이즈 방지)
  cc <- obj@meta.data %>%
    dplyr::group_by(donor_id, .data[[status_col]]) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop")
  valid <- dplyr::filter(cc, n >= min_cells)
  keep  <- paste(obj$donor_id, status_vec) %in% paste(valid$donor_id, valid[[status_col]])
  obj_f <- subset(obj, cells = colnames(obj)[keep])

  # 2. pseudobulk (raw count 합산).
  #    ★ Seurat AggregateExpression 은 버전에 따라 "layer matched by multiple actual arguments"
  #      버그가 있어 쓰지 않는다. 대신 희소 지시행렬(indicator) 곱으로 직접 합산한다.
  #      - donor_id 에 '_'/'-' 가 있어도 안전(구분자 '|@|' 사용)
  #      - BPCells 온디스크 행렬과도 그대로 동작(결과 genes×groups 는 작음)
  counts <- Seurat::GetAssayData(obj_f, assay = "RNA", layer = "counts")   # genes × cells
  grp    <- factor(paste(obj_f$donor_id, obj_f[[status_col]][, 1], sep = "|@|"))
  ind    <- Matrix::sparse.model.matrix(~ 0 + grp)                         # cells × groups
  colnames(ind) <- levels(grp)
  pb     <- as.matrix(counts %*% ind)                                      # genes × groups
  parts  <- do.call(rbind, strsplit(colnames(pb), "|@|", fixed = TRUE))
  coldata <- data.frame(sample = colnames(pb),
                        donor_id = parts[, 1], status = parts[, 2],
                        stringsAsFactors = FALSE)

  # 3. paired donor만 (양쪽 status 모두 있는 donor)
  donor_n <- table(coldata$donor_id)
  paired  <- names(donor_n[donor_n == 2])
  if (length(paired) < 2)
    stop("paired donor 가 2명 미만입니다 (paired = ", length(paired),
         "). min_cells(", min_cells, ") 를 낮추거나 status 라벨을 확인하세요.", call. = FALSE)
  sel    <- coldata$donor_id %in% paired
  pb_p   <- pb[, sel, drop = FALSE]
  cold_p <- coldata[sel, , drop = FALSE]

  # ★ 방향 고정: Positive vs Negative (coef 가 항상 statusPositive) + donor 는 present-only factor
  #   (서브셋 후 남는 unused level 로 인한 rank-deficient design → 계수 NA/통계 왜곡 방지)
  cold_p$status   <- factor(cold_p$status, levels = c(neg_level, pos_level))
  cold_p$donor_id <- factor(cold_p$donor_id)
  if (anyNA(cold_p$status))
    stop("status 값이 '", neg_level, "'/'", pos_level, "' 와 다릅니다: ",
         paste(unique(coldata$status), collapse = ", "), call. = FALSE)

  # 4. filterByExpr → voom → limma (donor 를 covariate 로 하는 paired design)
  dge   <- edgeR::DGEList(counts = pb_p)
  keepg <- edgeR::filterByExpr(dge, group = cold_p$status, min.count = min_count)
  dge   <- edgeR::calcNormFactors(edgeR::DGEList(counts = pb_p[keepg, , drop = FALSE]))
  design <- stats::model.matrix(~ donor_id + status, data = cold_p)
  v     <- limma::voom(dge, design)
  fit   <- limma::eBayes(limma::lmFit(v, design))

  coef_name <- paste0("status", pos_level)                       # = "statusPositive"
  res  <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
  list(fit = fit, coef = coef_name, paired_donors = paired, coldata = cold_p,
       n_genes = sum(keepg),
       result = tibble::rownames_to_column(as.data.frame(res), "gene"))
}

# ── cellchat용 rigorous cohort 구성 / build capped paired cohort ─────────────
# paired donor로 제한 + cell_type×donor 당 최대 N개 세포로 다운샘플.
build_cellchat_cohort <- function(seu, paired_donors, cap = CELLCHAT_CAP_PER_CT,
                                  seed = SUBSAMPLE_SEED) {
  mother <- subset(seu, donor_id %in% paired_donors)
  set.seed(seed)
  mother$cell_id <- rownames(mother@meta.data)
  keep <- mother@meta.data %>%
    dplyr::group_by(cell_type, donor_id) %>%
    dplyr::group_modify(~ dplyr::slice_sample(.x, n = min(nrow(.x), cap))) %>%
    dplyr::pull(cell_id)
  full <- subset(mother, cells = keep)
  label_slc_status(full)
}

# ── CellChat 객체 1개 실행 / run one CellChat pipeline ───────────────────────
run_cellchat <- function(seurat_obj, group_col = "cellchat_group", nboot = 100) {
  data.input <- as(Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = "data"), "dgCMatrix")  # BPCells→sparse(dgCMatrix)
  cc <- CellChat::createCellChat(object = data.input, meta = seurat_obj@meta.data, group.by = group_col)
  cc@DB <- CellChat::CellChatDB.human
  cc <- CellChat::subsetData(cc)
  cc <- CellChat::identifyOverExpressedGenes(cc)
  cc <- CellChat::identifyOverExpressedInteractions(cc)
  cc <- CellChat::computeCommunProb(cc, population.size = FALSE, nboot = nboot)
  cc <- CellChat::filterCommunication(cc, min.cells = 10)
  cc <- CellChat::computeCommunProbPathway(cc)
  cc <- CellChat::aggregateNet(cc)
  CellChat::netAnalysis_computeCentrality(cc)
}

# ── microglia를 SLC7A7 상태로 나눈 pos/neg 배경공유 객체 쌍 만들기 ───────────
# 비-microglia 배경은 동일하게 유지하고, microglia만 +/- 로 나눠 두 객체 생성.
make_pos_neg_objects <- function(seurat_full, status_col = paste0(GENE_OF_INTEREST, "_status")) {
  seurat_full$cellchat_group <- as.character(seurat_full$cell_type)
  is_mg <- seurat_full$cell_type == "microglial cell"
  seurat_full$cellchat_group[is_mg] <- "microglia"
  cells_pos <- colnames(seurat_full)[!(is_mg & seurat_full[[status_col]][, 1] == "Negative")]
  cells_neg <- colnames(seurat_full)[!(is_mg & seurat_full[[status_col]][, 1] == "Positive")]
  list(obj_pos = subset(seurat_full, cells = cells_pos),
       obj_neg = subset(seurat_full, cells = cells_neg),
       seurat_full = seurat_full)
}
