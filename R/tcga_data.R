# =============================================================================
# R/tcga_data.R  ---  TCGA 다운로드 및 발현행렬 준비 / download & prepare
# -----------------------------------------------------------------------------
# 원본 TCGA.R / TCGA_GBM.R 의 다운로드·전처리 로직을 project 인자 하나로 통합.
# (TCGA-LGG, TCGA-GBM 이 완전히 동일한 절차였음)
# =============================================================================

# ── GDC에서 STAR-counts 다운로드 후 SummarizedExperiment 반환 ────────────────
# directory: 원본 파일 보관 폴더(기본 data/raw/GDCdata). 이미 받아둔 파일이 있으면
# GDCdownload 가 재다운로드하지 않고 그대로 사용합니다.
download_tcga_expression <- function(project = c("TCGA-LGG", "TCGA-GBM"),
                                     directory       = DIR_GDC,
                                     method          = TCGA_DOWNLOAD_METHOD,
                                     files_per_chunk = TCGA_FILES_PER_CHUNK,
                                     max_retries     = 3) {
  project <- match.arg(project)
  # 다운로드 안정화: 타임아웃 확대 + libcurl(부분 다운로드/트렁케이션 완화)
  old_to <- getOption("timeout"); on.exit(options(timeout = old_to), add = TRUE)
  options(timeout = max(10000, old_to))
  dm <- getOption("download.file.method")
  if (is.null(dm) || dm %in% c("default", "internal", "wininet")) {
    options(download.file.method = "libcurl")
  }

  query <- TCGAbiolinks::GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = "Primary Tumor"
  )

  # 안정화:
  #  - method="client": gdc-client(공식 도구) 사용. 재개 가능·대용량 안정적. 트렁케이션에 강함.
  #  - method="api"   : files.per.chunk 로 덩어리를 잘게(작을수록 백신/프록시 컷오프 회피).
  #  두 경우 모두 실패 시 최대 max_retries 회 재시도.
  ok <- FALSE
  for (attempt in seq_len(max_retries)) {
    ok <- tryCatch({
      if (identical(method, "client")) {
        TCGAbiolinks::GDCdownload(query, directory = directory, method = "client")
      } else {
        TCGAbiolinks::GDCdownload(query, directory = directory,
                                  method = "api", files.per.chunk = files_per_chunk)
      }
      TRUE
    }, error = function(e) {
      message(sprintf("  [다운로드 실패 %d/%d] %s", attempt, max_retries, conditionMessage(e)))
      FALSE
    })
    if (ok) break
    if (attempt < max_retries) { message("  5초 후 재시도..."); Sys.sleep(5) }
  }
  if (!ok) stop(
    "GDCdownload 가 ", max_retries, "회 모두 실패했습니다.\n",
    "  1) 부분 다운로드 삭제: unlink(DIR_GDC, recursive = TRUE)\n",
    "  2) api 방식이면 config 의 TCGA_FILES_PER_CHUNK 를 1 로 낮추기\n",
    "  3) 백신 실시간검사 잠시 끄기 / 다른 네트워크(핫스팟) 시도\n",
    "  4) 그래도 안 되면 config 의 TCGA_DOWNLOAD_METHOD 를 \"client\" 로 (gdc-client 설치 필요)",
    call. = FALSE)

  TCGAbiolinks::GDCprepare(query, directory = directory)
}

# ── 캐시 우선 로더 / cache-first loader ──────────────────────────────────────
# 준비된 {tag}_counts_clinical.rds 가 있으면 그걸 읽고(다운로드 X),
# 없을 때만 GDC에서 받아 준비한 뒤 .rds 로 저장합니다. force=TRUE 면 강제 재다운로드.
get_tcga_data <- function(project = c("TCGA-LGG", "TCGA-GBM"), force = FALSE) {
  project <- match.arg(project)
  tag <- sub("TCGA-", "", project)
  rds <- file.path(DIR_DATA_PROC, paste0(tag, "_counts_clinical.rds"))
  if (!force && file.exists(rds)) {
    message("[cache] 로컬 저장본 사용 → ", rds, "  (GDC 재다운로드 안 함)")
    return(readRDS(rds))
  }
  message("[download] GDC에서 ", project, " 다운로드/준비 (최초 1회) ...")
  t0  <- Sys.time()
  out <- prepare_count_matrix(download_tcga_expression(project))

  # ── 다운로드 provenance 기록 (일시 등) ──────────────────────────────────────
  out$download_info <- list(
    project              = project,
    downloaded_at        = format(t0, "%Y-%m-%d %H:%M:%S %Z"),
    n_samples            = ncol(out$counts),
    n_genes              = nrow(out$counts),
    tcgabiolinks_version = as.character(utils::packageVersion("TCGAbiolinks")),
    gdc_directory        = normalizePath(DIR_GDC, mustWork = FALSE)
  )
  saveRDS(out, rds)

  # 사람이 읽는 provenance 텍스트도 남김 (data/processed/{tag}_download_info.txt)
  info_txt <- file.path(DIR_DATA_PROC, paste0(tag, "_download_info.txt"))
  writeLines(c(
    paste("project:              ", out$download_info$project),
    paste("downloaded_at:        ", out$download_info$downloaded_at),
    paste("n_samples:            ", out$download_info$n_samples),
    paste("n_genes:              ", out$download_info$n_genes),
    paste("TCGAbiolinks_version: ", out$download_info$tcgabiolinks_version)
  ), info_txt)

  message("[saved] ", rds, "\n[log]   ", info_txt,
          " (다운로드 일시: ", out$download_info$downloaded_at, ")")
  out
}

# ── SE → (count_matrix_symbol, clinical_data) / raw counts + SYMBOL rownames ─
# Primary Tumor 필터, ENSEMBL 버전 접미사 제거, SYMBOL 매핑(중복 첫 항목 유지).
prepare_count_matrix <- function(se) {
  clinical <- SummarizedExperiment::colData(se)
  count_matrix <- SummarizedExperiment::assay(se)

  primary_idx  <- clinical$sample_type == "Primary Tumor"
  count_matrix <- count_matrix[, primary_idx]
  clinical     <- as.data.frame(clinical[primary_idx, ])

  rownames(count_matrix) <- sub("\\..*$", "", rownames(count_matrix))  # ENSG버전 제거

  gene_map <- AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = rownames(count_matrix), keytype = "ENSEMBL", columns = "SYMBOL"
  )
  gene_map <- gene_map[!is.na(gene_map$SYMBOL), ]
  gene_map <- gene_map[!duplicated(gene_map$SYMBOL), ]

  count_matrix_symbol <- count_matrix[gene_map$ENSEMBL, ]
  rownames(count_matrix_symbol) <- gene_map$SYMBOL

  list(counts = count_matrix_symbol, clinical = clinical)
}

# ── raw counts → logCPM 정규화 / edgeR log-CPM normalization ─────────────────
# 사용자 원본(egdeR.R)과 동일: edgeR::cpm(counts, log = TRUE, prior.count = 1).
# 입력 counts: SYMBOL(rowname) × sample 정수 count 행렬 (prepare_count_matrix 산출).
# 반환: data.frame (첫 열 Gene = SYMBOL, 이후 sample별 logCPM) — load_expr_matrix 와 호환.
normalize_logcpm <- function(counts, prior_count = 1) {
  logcpm <- edgeR::cpm(counts, log = TRUE, prior.count = prior_count)
  out <- as.data.frame(logcpm)
  tibble::rownames_to_column(out, var = "Gene")
}

# 이미 저장해둔 rawcounts .txt(첫 행 헤더 + 유전자명 열)를 count 행렬로 되읽기.
# (원본 egdeR.R 의 read.table → header/gene 처리 흐름을 함수화)
read_rawcounts_txt <- function(path) {
  df <- utils::read.table(path)
  names(df) <- df[1, ]; df <- df[-1, ]
  gene <- df[[1]]; df <- df[, -1, drop = FALSE]
  df[] <- lapply(df, as.numeric); rownames(df) <- gene
  as.matrix(df)
}

# ── 정규화 발현행렬을 clinical에 붙이기 / attach gene logCPM to clinical ──────
# normalized_counts_final: load_expr_matrix() 결과 (Gene + samples, sample 순서 = clinical 순서 가정)
attach_gene_to_clinical <- function(clinical_df, expr_df, gene = GENE_OF_INTEREST) {
  row <- expr_df$Gene == gene
  clinical_df[[paste0(gene, "_logCPM")]] <- as.numeric(expr_df[row, ])[-1]
  clinical_df[[paste0(gene, "_group")]] <- ifelse(
    clinical_df[[paste0(gene, "_logCPM")]] >
      stats::median(clinical_df[[paste0(gene, "_logCPM")]], na.rm = TRUE),
    "High", "Low"
  )
  clinical_df
}
