# =============================================================================
# R/deconvolution.R  ---  면역 deconvolution / bulk 신호의 myeloid 귀속 (Fig 4)
# "bulk SLC7A7 신호가 myeloid 분율로 설명되는가?" — 선행연구 미실시 패널.
# =============================================================================
run_estimate <- function(expr, tmp_dir = tempdir()) {
  stopifnot(requireNamespace("estimate", quietly = TRUE))
  f_in <- file.path(tmp_dir, "expr_input.txt"); f_gct <- file.path(tmp_dir, "expr.gct")
  f_out <- file.path(tmp_dir, "estimate.gct")
  utils::write.table(data.frame(GeneSymbol = rownames(expr), expr, check.names = FALSE),
                     f_in, sep = "\t", quote = FALSE, row.names = FALSE)
  estimate::filterCommonGenes(input.f = f_in, output.f = f_gct, id = "GeneSymbol")
  estimate::estimateScore(f_gct, f_out, platform = "illumina")
  sc <- utils::read.delim(f_out, skip = 2, row.names = 1, check.names = FALSE)[, -1, drop = FALSE]
  out <- as.data.frame(t(sc)); out$sample <- rownames(out); out
}

run_xcell <- function(expr) {
  if (!requireNamespace("xCell", quietly = TRUE)) {
    message("xCell 미설치 → 건너뜀 (remotes::install_github('dviraran/xCell'))"); return(NULL) }
  out <- as.data.frame(t(xCell::xCellAnalysis(expr))); out$sample <- rownames(out); out
}

load_cibersortx <- function(file) {
  d <- utils::read.csv(require_data(file), check.names = FALSE); names(d)[1] <- "sample"; d }

correlate_gene_fractions <- function(expr, fractions, gene = GENE_OF_INTEREST, method = "spearman") {
  if (!gene %in% rownames(expr)) return(NULL)
  v  <- as.numeric(expr[gene, ])
  fr <- fractions[match(colnames(expr), fractions$sample), , drop = FALSE]
  res <- dplyr::bind_rows(lapply(setdiff(names(fr), "sample"), function(cc) {
    y <- suppressWarnings(as.numeric(fr[[cc]]))
    if (all(is.na(y)) || stats::sd(y, na.rm = TRUE) == 0) return(NULL)
    ct <- suppressWarnings(stats::cor.test(v, y, method = method))
    data.frame(cell_type = cc, rho = unname(ct$estimate), p = ct$p.value)
  }))
  res$rho <- as.numeric(res$rho)
  res$p   <- as.numeric(res$p)
  
  res %>%
    dplyr::mutate(FDR = p.adjust(p, "BH")) %>%
    dplyr::arrange(p)
}

# myeloid 분율 보정 후에도 예후 효과가 남는가 (매개 점검)
adjust_for_fraction <- function(expr, clin, fractions, gene = GENE_OF_INTEREST,
                                fraction_col, covars = c("age", "grade", "idh")) {
  v  <- as.numeric(scale(as.numeric(expr[gene, ])))
  fr <- fractions[match(colnames(expr), fractions$sample), fraction_col]
  d  <- data.frame(time = clin$time, event = clin$event, gene = v,
                   myeloid = as.numeric(scale(suppressWarnings(as.numeric(fr)))),
                   clin[, intersect(covars, names(clin)), drop = FALSE])
  keep <- vapply(d, function(x) length(unique(stats::na.omit(x))) > 1, logical(1))
  d <- d[, keep, drop = FALSE]; d <- d[stats::complete.cases(d) & d$time > 0, ]
  rhs <- setdiff(names(d), c("time", "event"))
  s <- summary(survival::coxph(stats::as.formula(paste("survival::Surv(time,event) ~",
       paste(rhs, collapse = "+"))), data = d))
  data.frame(term = rownames(s$conf.int), HR = s$conf.int[, 1], lower = s$conf.int[, 3],
             upper = s$conf.int[, 4], p = s$coefficients[, 5], row.names = NULL)
}
