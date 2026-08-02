# =============================================================================
# R/utils.R  ---  공통 유틸리티 / Shared utilities
# =============================================================================

# ── 발현행렬 로드 헬퍼 / Load a gene x sample expression matrix ───────────────
# 첫 열을 유전자명(Gene)으로 잡고, 샘플명의 '.'을 '-'로 되돌립니다.
# (read.csv 가 barcode의 '-'를 '.'로 바꾸는 문제를 원복)
load_expr_matrix <- function(path, restore_dash = TRUE) {
  m <- utils::read.csv(path, check.names = FALSE)
  names(m)[1] <- "Gene"
  if (restore_dash) {
    names(m)[-1] <- gsub("\\.", "-", names(m)[-1])
  }
  m
}

# 발현 data.frame(Gene + samples) → numeric matrix (rownames = Gene)
expr_df_to_matrix <- function(expr_df) {
  mat <- as.matrix(expr_df[, -1, drop = FALSE])
  rownames(mat) <- expr_df$Gene
  mode(mat) <- "numeric"
  mat
}

# ── ENSEMBL → SYMBOL 등 ID 매핑 / gene ID mapping via org.Hs.eg.db ────────────
map_gene_ids <- function(keys, from = "SYMBOL", to = "ENSEMBL") {
  AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = keys, column = to, keytype = from, multiVals = "first"
  )
}

# ── DEG 유전자명 정리 / drop MT-, ENSG, ribosomal pseudogenes ────────────────
clean_deg_genes <- function(df, drop_self = NULL) {
  keep <- !grepl("^MT-?", rownames(df)) &
          !grepl("^ENSG",  rownames(df)) &
          !grepl("^RP[LS][0-9]", rownames(df))
  out <- df[keep, , drop = FALSE]
  if (!is.null(drop_self)) out <- out[rownames(out) != drop_self, , drop = FALSE]
  out
}

# ── 논문용 ggplot 테마 / Publication theme ───────────────────────────────────
theme_paper <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text  = ggplot2::element_text(color = "black"),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.line  = ggplot2::element_line(linewidth = 0.5),
      legend.position = "right"
    )
}

# ── "그림 틀만" 테마 / Blank-frame theme ──────────────────────────────────────
# 축·제목·텍스트를 모두 제거하고 tick만 남깁니다. PPT로 옮겨 라벨을 다시 그리는
# 논문 제출용 워크플로에 사용. (원본 seurat그래프그리는코드.R의 스타일을 함수화)
theme_blank_frame <- function(tick_len = 8, keep_grid_x = FALSE) {
  th <- ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none",
                   plot.title = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_line(linewidth = 0.5, color = "black"),
                   axis.ticks.y = ggplot2::element_line(linewidth = 0.5, color = "black"),
                   axis.ticks.length.x = ggplot2::unit(tick_len, "pt"),
                   axis.ticks.length.y = ggplot2::unit(tick_len, "pt"))
  if (keep_grid_x) {
    th <- th + ggplot2::theme(
      panel.grid.major.x = ggplot2::element_line(color = "gray92", linewidth = 0.3)
    )
  }
  th
}


# ── 그림 저장 헬퍼 / Figure export helpers ───────────────────────────────────
# EMF (벡터, PPT 편집용) 저장. devEMF 필요.
save_emf <- function(plot, filename, width = 7, height = 7, dir = DIR_FIGURES, draw_grid = FALSE) {
  plot <- plot + theme_blank_frame(keep_grid_x = draw_grid)
  devEMF::emf(file = file.path(dir, filename), width = width, height = height)
  print(plot); grDevices::dev.off()
  invisible(file.path(dir, filename))
}

save_tiff <- function(plot, filename, width = 7, height = 7, dpi = 300, dir = DIR_FIGURES) {
  ggplot2::ggsave(filename, plot = plot, path = dir,
                  width = width, height = height, dpi = dpi, compression = "lzw")
  invisible(file.path(dir, filename))
}

# ── 결과 CSV 저장 헬퍼 ───────────────────────────────────────────────────────
write_result <- function(x, filename, dir = DIR_RESULTS, row.names = FALSE) {
  utils::write.csv(x, file.path(dir, filename), row.names = row.names)
  invisible(file.path(dir, filename))
}
