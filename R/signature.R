# =============================================================================
# R/signature.R  ---  LASSO Cox 기반 유전자 선정 / penalised Cox selection
# -----------------------------------------------------------------------------
# Fig 2: 수송체 gene set → LASSO Cox(10-fold CV) → 선택 유전자 → 단/다변량 Cox.
# ★ SLC7A7 이 LASSO 에서 탈락할 수 있다. 그 경우 candidate-gene 설계
#   ("선행 단일세포 관찰에서 도출된 사전지정 후보")로 전환한다. 양쪽 모두 지원.
# =============================================================================

run_lasso_cox <- function(expr, clin, genes, seed = LASSO_SEED, nfolds = LASSO_NFOLDS,
                          alpha = LASSO_ALPHA, lambda_choice = LASSO_LAMBDA) {
  genes <- intersect(genes, rownames(expr))
  ok <- !is.na(clin$time) & !is.na(clin$event) & clin$time > 0
  if (sum(ok) < 30) stop("생존정보 유효 표본이 부족합니다 (n=", sum(ok), ").", call. = FALSE)
  x <- t(expr[genes, ok, drop = FALSE])
  x <- x[, apply(x, 2, stats::sd, na.rm = TRUE) > 0, drop = FALSE]
  y <- survival::Surv(clin$time[ok], clin$event[ok])

  set.seed(seed)
  cvfit <- glmnet::cv.glmnet(x, y, family = "cox", alpha = alpha, nfolds = nfolds, standardize = TRUE)
  lam <- cvfit[[lambda_choice]]
  co  <- as.matrix(stats::coef(cvfit, s = lam))
  sel <- data.frame(gene = rownames(co)[co[, 1] != 0], coef = co[co[, 1] != 0, 1],
                    stringsAsFactors = FALSE)
  sel <- sel[order(-abs(sel$coef)), ]
  message("LASSO 선택 유전자: ", nrow(sel), " (lambda=", signif(lam, 3), ", n=", sum(ok), ")")
  list(cvfit = cvfit, fit = cvfit$glmnet.fit, lambda = lam, selected = sel, x = x, y = y)
}

plot_lasso <- function(res) {
  op <- graphics::par(mfrow = c(1, 2)); on.exit(graphics::par(op), add = TRUE)
  graphics::plot(res$fit, xvar = "lambda", label = TRUE); graphics::title("LASSO path", line = 2.5)
  graphics::plot(res$cvfit); graphics::title("10-fold CV", line = 2.5)
}

compute_risk_score <- function(expr, selected) {
  g <- intersect(selected$gene, rownames(expr))
  if (!length(g)) return(rep(NA_real_, ncol(expr)))
  as.numeric(crossprod(expr[g, , drop = FALSE], selected$coef[match(g, selected$gene)]))
}

univariate_cox <- function(expr, clin, genes, scale_expr = TRUE) {
  genes <- intersect(genes, rownames(expr))
  res <- lapply(genes, function(g) {
    v <- as.numeric(expr[g, ]); if (scale_expr) v <- as.numeric(scale(v))
    d <- data.frame(time = clin$time, event = clin$event, x = v)
    d <- d[stats::complete.cases(d) & d$time > 0, ]
    if (nrow(d) < 20 || stats::sd(d$x) == 0) return(NULL)
    s <- summary(survival::coxph(survival::Surv(time, event) ~ x, data = d))
    data.frame(gene = g, HR = s$conf.int[1, 1], lower = s$conf.int[1, 3],
               upper = s$conf.int[1, 4], p = s$coefficients[1, 5])
  })
  res <- dplyr::bind_rows(res)
  dplyr::arrange(dplyr::mutate(res, FDR = stats::p.adjust(p, "BH")), p)
}

# IDH·grade·age·MGMT 보정 다변량 Cox (Fig 3 핵심)
multivariable_cox <- function(expr, clin, gene = GENE_OF_INTEREST,
                              covars = c("age", "grade", "idh", "mgmt"), scale_expr = TRUE) {
  if (!gene %in% rownames(expr)) return(NULL)
  v <- as.numeric(expr[gene, ]); if (scale_expr) v <- as.numeric(scale(v))
  d <- data.frame(time = clin$time, event = clin$event, gene = v,
                  clin[, intersect(covars, names(clin)), drop = FALSE])
  keep <- vapply(d, function(x) length(unique(stats::na.omit(x))) > 1, logical(1))
  d <- d[, keep, drop = FALSE]
  d <- d[stats::complete.cases(d) & d$time > 0, ]
  if (nrow(d) < 30) { message("  표본 부족(n=", nrow(d), ")"); return(NULL) }
  rhs <- setdiff(names(d), c("time", "event"))
  fit <- survival::coxph(stats::as.formula(paste("survival::Surv(time, event) ~",
                                                 paste(rhs, collapse = " + "))), data = d)
  s <- summary(fit)
  list(fit = fit, table = data.frame(term = rownames(s$conf.int), HR = s$conf.int[, 1],
       lower = s$conf.int[, 3], upper = s$conf.int[, 4], p = s$coefficients[, 5],
       n = fit$n, row.names = NULL))
}

stratified_cox <- function(expr, clin, gene = GENE_OF_INTEREST, by = "idh") {
  lv <- stats::na.omit(unique(clin[[by]]))
  dplyr::bind_rows(lapply(lv, function(l) {
    idx <- which(clin[[by]] == l)
    r <- multivariable_cox(expr[, idx, drop = FALSE], clin[idx, , drop = FALSE], gene,
                           covars = setdiff(c("age", "grade", "mgmt"), by))
    if (is.null(r)) return(NULL)
    dplyr::mutate(r$table[r$table$term == "gene", ], stratum = paste0(by, "=", l))
  }))
}
