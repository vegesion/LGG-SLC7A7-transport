# =============================================================================
# R/validation.R  ---  예후 검증 / KM · time-dependent ROC · DCA · forest
# Fig 3: 다중 코호트에서 동일 절차 반복 (TCGA LGG/GBM, CGGA 693+325, GEO)
# =============================================================================
`%||%` <- function(a, b) if (is.null(a)) b else a

make_groups <- function(x, method = c("median", "optimal"), time = NULL, event = NULL) {
  method <- match.arg(method)
  cut <- if (method == "optimal" && !is.null(time)) {
    survminer::surv_cutpoint(data.frame(time = time, event = event, x = x),
      time = "time", event = "event", variables = "x")$cutpoint$cutpoint
  } else stats::median(x, na.rm = TRUE)
  factor(ifelse(x > cut, "High", "Low"), levels = c("Low", "High"))
}

km_cohort <- function(co, gene = GENE_OF_INTEREST, method = "median", palette = PALETTE_TWO) {
  d <- data.frame(time = co$clin$time, event = co$clin$event, x = cohort_gene(co, gene))
  d <- d[stats::complete.cases(d) & d$time > 0, ]
  d$group <- make_groups(d$x, method, d$time, d$event)
  fit <- survival::survfit(survival::Surv(time, event) ~ group, data = d)
  survminer::ggsurvplot(fit, data = d, pval = TRUE, conf.int = FALSE, risk.table = TRUE,
    palette = palette, censor.shape = 124, size = 1.1, xlab = "Time (months)",
    ylab = "Overall survival", title = paste0(co$name, " — ", gene),
    ggtheme = ggplot2::theme_classic(base_size = 12))
}

timedep_roc <- function(co, marker, times = ROC_TIMES, label = NULL) {
  d <- data.frame(time = co$clin$time, event = co$clin$event, m = marker)
  d <- d[stats::complete.cases(d) & d$time > 0, ]
  times <- times[times < max(d$time, na.rm = TRUE)]
  if (!length(times)) return(NULL)
  r <- timeROC::timeROC(T = d$time, delta = d$event, marker = d$m, cause = 1, times = times)
  data.frame(cohort = co$name, marker = label %||% "score", time = times, AUC = as.numeric(r$AUC))
}

plot_timedep_roc <- function(co, marker, times = ROC_TIMES) {
  d <- data.frame(time = co$clin$time, event = co$clin$event, m = marker)
  d <- d[stats::complete.cases(d) & d$time > 0, ]
  times <- times[times < max(d$time, na.rm = TRUE)]
  r <- timeROC::timeROC(T = d$time, delta = d$event, marker = d$m, cause = 1, times = times)
  df <- dplyr::bind_rows(lapply(seq_along(times), function(i)
    data.frame(FPR = r$FP[, i], TPR = r$TP[, i],
               t = sprintf("%g mo (AUC=%.2f)", times[i], r$AUC[i]))))
  ggplot2::ggplot(df, ggplot2::aes(FPR, TPR, color = t)) +
    ggplot2::geom_line(linewidth = 1) + ggplot2::geom_abline(linetype = 2, color = "grey60") +
    ggplot2::labs(title = co$name, x = "1 - Specificity", y = "Sensitivity", color = NULL) +
    theme_paper()
}

# Decision Curve Analysis: 유전자 추가가 임상적 유용성을 높이는가
dca_cohort <- function(co, gene = GENE_OF_INTEREST, horizon = 36,
                       covars = c("age", "grade", "idh")) {
  cl <- co$clin
  d <- data.frame(time = cl$time, event = cl$event, gene = cohort_gene(co, gene),
                  cl[, intersect(covars, names(cl)), drop = FALSE])
  keep <- vapply(d, function(x) length(unique(stats::na.omit(x))) > 1, logical(1))
  d <- d[, keep, drop = FALSE]; d <- d[stats::complete.cases(d) & d$time > 0, ]
  if (nrow(d) < 30 || !"gene" %in% names(d)) return(NULL)
  d$outcome <- as.integer(d$time <= horizon & d$event == 1)
  rhs <- setdiff(names(d), c("time", "event", "outcome"))
  m_full <- stats::glm(stats::as.formula(paste("outcome ~", paste(rhs, collapse = "+"))),
                       data = d, family = "binomial")
  m_base <- stats::glm(stats::as.formula(paste("outcome ~",
    paste(setdiff(rhs, "gene"), collapse = "+"))), data = d, family = "binomial")
  th <- seq(0.01, 0.60, by = 0.01)
  nb <- function(p) vapply(th, function(t) {
    tp <- sum(p >= t & d$outcome == 1); fp <- sum(p >= t & d$outcome == 0)
    tp / nrow(d) - fp / nrow(d) * (t / (1 - t)) }, numeric(1))
  prev <- mean(d$outcome)
  dplyr::bind_rows(
    data.frame(threshold = th, net_benefit = nb(stats::predict(m_full, type = "response")),
               model = paste0(gene, " + clinical")),
    data.frame(threshold = th, net_benefit = nb(stats::predict(m_base, type = "response")),
               model = "clinical only"),
    data.frame(threshold = th, net_benefit = prev - (1 - prev) * (th / (1 - th)), model = "treat all"),
    data.frame(threshold = th, net_benefit = 0, model = "treat none")) %>%
    dplyr::mutate(cohort = co$name)
}

plot_dca <- function(dca_df) {
  ggplot2::ggplot(dca_df, ggplot2::aes(threshold, net_benefit, color = model)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::coord_cartesian(ylim = c(-0.05, max(dca_df$net_benefit, na.rm = TRUE) * 1.1)) +
    ggplot2::facet_wrap(~cohort, scales = "free_y") +
    ggplot2::labs(x = "Threshold probability", y = "Net benefit", color = NULL) + theme_paper()
}

forest_multicohort <- function(tab) {
  ggplot2::ggplot(tab, ggplot2::aes(HR, stats::reorder(label, HR))) +
    ggplot2::geom_vline(xintercept = 1, color = "grey60") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lower, xmax = upper), height = 0.2, linewidth = 0.9) +
    ggplot2::geom_point(size = 2.6, color = PALETTE_TWO[1]) +
    ggplot2::scale_x_continuous(trans = "log2") +
    ggplot2::labs(x = "Adjusted HR (per SD)", y = NULL) + theme_paper()
}
