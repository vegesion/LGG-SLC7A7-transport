# =============================================================================
# R/survival.R  ---  TCGA bulk 생존분석 함수 / Cox survival functions
# -----------------------------------------------------------------------------
# 원본 TCGA_allgenes.R / TCGA_GBM.R / TCGA_single.R 의 로직을 함수로 통합.
# LGG와 GBM은 subgroup 정의와 covariate만 다르므로, cohort 설정(cohort spec)을
# 인자로 받아 하나의 스크리닝 엔진(run_cox_screen)으로 처리합니다.
# 통계 절차(Z-score 표준화, 이벤트 임계값, Firth fallback, BH FDR)는 원본과 동일.
# =============================================================================

# ── 표준 생존식 좌변 / survival LHS ──────────────────────────────────────────
.SURV_LHS <- "Surv(paper_Survival..months., paper_Vital.status..1.dead.)"

# ── 안전한 Cox 적합 (separation 시 Firth) / safe Cox with Firth fallback ─────
# 반환: HR, lower95, upper95, p, model 한 행. 실패 시 NULL.
safe_cox <- function(formula, data) {
  warn_msg <- NULL
  fit <- withCallingHandlers(
    tryCatch(survival::coxph(formula, data = data), error = function(e) NULL),
    warning = function(w) { warn_msg <<- conditionMessage(w); invokeRestart("muffleWarning") }
  )

  # 무한대 계수 경고 → Firth penalized Cox (coxphf)
  if (!is.null(warn_msg) && grepl("infinite", warn_msg)) {
    fit <- tryCatch(coxphf::coxphf(formula, data = data), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    return(data.frame(
      HR = exp(fit$coefficients[1]), lower95 = exp(fit$ci.lower[1]),
      upper95 = exp(fit$ci.upper[1]), p = fit$prob[1], model = "firth"
    ))
  }
  if (is.null(fit)) return(NULL)
  s <- summary(fit)
  data.frame(
    HR = s$conf.int[1, "exp(coef)"], lower95 = s$conf.int[1, "lower .95"],
    upper95 = s$conf.int[1, "upper .95"], p = s$coefficients[1, "Pr(>|z|)"], model = "cox"
  )
}

# ── cohort 설정 정의 / cohort specifications ─────────────────────────────────
# 각 subgroup: 필터(filter), covariate 문자열, 최소 이벤트 수.
# exposure 변수는 항상 logCPM_scaled (Z-score 표준화된 발현).
cohort_spec_LGG <- function() list(
  overall = list(filter = NULL,
                 covars  = "age_at_diagnosis + paper_IDH.status + paper_X1p.19q.codeletion + tumor_grade",
                 min_events = COX_MIN_EVENTS_LGG_OVERALL),
  mut     = list(filter = quote(paper_IDH.status == "Mutant"),
                 covars  = "age_at_diagnosis + paper_X1p.19q.codeletion + tumor_grade",
                 min_events = COX_MIN_EVENTS_LGG_MUT),
  wt      = list(filter = quote(paper_IDH.status == "WT"),
                 covars  = "age_at_diagnosis + tumor_grade",
                 min_events = COX_MIN_EVENTS_LGG_WT)
)

cohort_spec_GBM <- function() list(
  overall = list(filter = NULL,
                 covars  = "age_at_diagnosis + paper_MGMT.promoter.status",
                 min_events = COX_MIN_EVENTS_GBM),
  mut     = list(filter = quote(paper_IDH.status == "Mutant"),
                 covars  = "age_at_diagnosis",
                 min_events = COX_MIN_EVENTS_GBM),
  wt      = list(filter = quote(paper_IDH.status == "WT"),
                 covars  = "age_at_diagnosis",
                 min_events = COX_MIN_EVENTS_GBM)
)

# ── 유전자 1개에 대한 subgroup별 Cox / per-gene Cox across subgroups ──────────
analyze_gene_cox <- function(gene, merged_df, spec) {
  df <- dplyr::filter(merged_df, Gene == gene)
  if (isTRUE(COX_SCALE_EXPR)) df$logCPM_scaled <- scale(df$logCPM)[, 1]
  else                         df$logCPM_scaled <- df$logCPM

  fit_subgroup <- function(sg) {
    d <- if (is.null(sg$filter)) df else df[which(eval(sg$filter, envir = df)), ]
    events <- sum(d$paper_Vital.status..1.dead. == 1, na.rm = TRUE)
    if (events < sg$min_events) return(NULL)
    f <- stats::as.formula(paste0(.SURV_LHS, " ~ logCPM_scaled + ", sg$covars))
    safe_cox(f, d)
  }
  res <- lapply(spec, fit_subgroup)

  pick <- function(r, col) if (is.null(r)) NA else r[[col]]
  data.frame(
    gene = gene,
    overall_HR = pick(res$overall, "HR"), overall_lower = pick(res$overall, "lower95"),
    overall_upper = pick(res$overall, "upper95"), overall_p = pick(res$overall, "p"),
    mut_HR = pick(res$mut, "HR"), mut_lower = pick(res$mut, "lower95"),
    mut_upper = pick(res$mut, "upper95"), mut_p = pick(res$mut, "p"),
    wt_HR = pick(res$wt, "HR"), wt_lower = pick(res$wt, "lower95"),
    wt_upper = pick(res$wt, "upper95"), wt_p = pick(res$wt, "p")
  )
}

# ── 유전자 패밀리 전체 스크리닝 / screen a gene family ───────────────────────
# merged_df: expr_long(Gene, Sample, logCPM) 을 clinical에 inner_join 한 것.
# 반환: 유전자별 HR/CI/p + BH FDR 표.
run_cox_screen <- function(merged_df, spec, family = GENE_FAMILY) {
  genes <- unique(merged_df$Gene)
  if (!is.null(family)) genes <- genes[stringr::str_detect(genes, family)]
  res <- pbapply::pblapply(genes, analyze_gene_cox, merged_df = merged_df, spec = spec)
  res <- dplyr::bind_rows(res)
  dplyr::mutate(res,
    FDR_overall = stats::p.adjust(overall_p, "BH"),
    FDR_mutant  = stats::p.adjust(mut_p, "BH"),
    FDR_wt      = stats::p.adjust(wt_p, "BH")
  )
}

# ── expr_long × clinical 병합 헬퍼 / build merged long table ─────────────────
build_merged_df <- function(expr_df, clinical_df) {
  expr_long <- tidyr::pivot_longer(expr_df, cols = -Gene,
                                   names_to = "Sample", values_to = "logCPM")
  clin <- dplyr::rename(as.data.frame(clinical_df), Sample = barcode)
  dplyr::inner_join(expr_long, clin, by = "Sample")
}

# ── Forest plot / forest plot of subgroup HRs ────────────────────────────────
# group_cols: list of list(hr, lower, upper, p, label). 원본 plot_forest 로직 유지.
plot_forest <- function(df, tumor_type, group_cols, log_scale = FALSE) {
  plot_data <- purrr::map_dfr(group_cols, function(g) {
    df %>% dplyr::select(gene, HR = dplyr::all_of(g$hr), lower = dplyr::all_of(g$lower),
                         upper = dplyr::all_of(g$upper), p = dplyr::all_of(g$p)) %>%
      dplyr::mutate(group = g$label, HR_log = HR, CI_lo = lower, CI_hi = upper)
  })
  gene_order <- df %>% dplyr::arrange(overall_HR) %>% dplyr::pull(gene)
  plot_data$gene  <- factor(plot_data$gene,  levels = gene_order)
  plot_data$group <- factor(plot_data$group, levels = purrr::map_chr(group_cols, "label"))

  ggplot2::ggplot(plot_data, ggplot2::aes(x = HR_log, y = gene, color = group)) +
    ggplot2::geom_vline(xintercept = 1, color = "grey60", linewidth = 0.5) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = CI_lo, xmax = CI_hi), height = 0.25,
                            linewidth = 1, position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::geom_point(size = 2.5, position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::scale_color_manual(values = PALETTE_GROUP, name = "Subgroup") +
    ggplot2::scale_x_continuous(name = "Hazard Ratio") +
    ggplot2::labs(y = NULL, title = paste0(tumor_type, " — ", GENE_FAMILY, " family Cox HR")) +
    theme_paper(base_size = 10) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9, face = "italic"),
                   panel.grid.major.x = ggplot2::element_line(color = "gray92", linewidth = 0.3))
}

# subgroup 컬럼 매핑 (LGG=3군, GBM=2군)
forest_groups_lgg <- list(
  list(hr = "overall_HR", lower = "overall_lower", upper = "overall_upper", p = "overall_p", label = "Overall"),
  list(hr = "mut_HR", lower = "mut_lower", upper = "mut_upper", p = "mut_p", label = "Mutant"),
  list(hr = "wt_HR",  lower = "wt_lower",  upper = "wt_upper",  p = "wt_p",  label = "WT")
)
forest_groups_gbm <- list(
  list(hr = "overall_HR", lower = "overall_lower", upper = "overall_upper", p = "overall_p", label = "Overall"),
  list(hr = "wt_HR", lower = "wt_lower", upper = "wt_upper", p = "wt_p", label = "WT")
)

# ── 단일 유전자: 발현 추가 + median 그룹화 / add expression & median split ────
add_gene_expression <- function(gene_symbol, expr_matrix, clinical_df) {
  row_idx <- which(expr_matrix$Gene == gene_symbol)
  if (length(row_idx) == 0) { warning(paste(gene_symbol, "not found")); return(clinical_df) }
  expr_values <- as.numeric(expr_matrix[row_idx, -1])
  if (length(expr_values) != nrow(clinical_df))
    stop(paste(gene_symbol, ": expression length and clinical sample size mismatch"))
  med <- stats::median(expr_values, na.rm = TRUE)
  clinical_df[[paste0(gene_symbol, "_logCPM")]] <- expr_values
  clinical_df[[paste0(gene_symbol, "_group")]]  <-
    factor(ifelse(expr_values >= med, "High", "Low"), levels = c("Low", "High"))
  clinical_df
}

# ── 단일 유전자 상세 검정 / detailed single-gene tests ───────────────────────
# Wilcoxon(발현 vs IDH), chi-sq(group vs IDH), uni/multi/spline Cox, IDH subgroup.
analyze_gene_single <- function(gene, clinical_df) {
  expr_var <- paste0(gene, "_logCPM"); group_var <- paste0(gene, "_group")
  if (!expr_var %in% names(clinical_df) || !group_var %in% names(clinical_df)) {
    message(gene, ": expression/group variable not found"); return(invisible(NULL))
  }
  clinical_df[[group_var]] <- factor(clinical_df[[group_var]], levels = c("Low", "High"))
  f <- function(rhs, data = clinical_df) stats::as.formula(paste0(.SURV_LHS, " ~ ", rhs))

  cat("\n===== Gene:", gene, "=====\n")
  cat("\n[Wilcoxon: expression ~ IDH]\n")
  print(stats::wilcox.test(clinical_df[[expr_var]] ~ clinical_df$paper_IDH.status))
  cat("\n[Chi-square: group vs IDH]\n")
  print(stats::chisq.test(table(clinical_df[[group_var]], clinical_df$paper_IDH.status)))
  cat("\n[Univariable Cox]\n");   print(summary(survival::coxph(f(expr_var), clinical_df)))
  cat("\n[Multivariable Cox]\n"); print(summary(survival::coxph(
    f(paste0(expr_var, " + age_at_diagnosis + paper_IDH.status + paper_X1p.19q.codeletion + tumor_grade")), clinical_df)))
  cat("\n[Non-linear Cox: P-spline]\n"); print(summary(survival::coxph(f(paste0("pspline(", expr_var, ")")), clinical_df)))

  for (grp in c("Mutant", "WT")) {
    sub <- subset(clinical_df, paper_IDH.status == grp)
    cov <- if (grp == "Mutant") "age_at_diagnosis + paper_X1p.19q.codeletion + tumor_grade"
           else                  "age_at_diagnosis + tumor_grade"
    cat("\n--- [IDH-", grp, " subgroup] N =", nrow(sub),
        "Events =", sum(sub$paper_Vital.status..1.dead., na.rm = TRUE), "---\n")
    print(summary(survival::coxph(f(expr_var), sub)))
    print(summary(survival::coxph(f(paste0(expr_var, " + ", cov)), sub)))
  }
  invisible(NULL)
}

# ── Kaplan–Meier plot (High vs Low) ──────────────────────────────────────────
plot_km <- function(gene, clinical_df, group_col = NULL, palette = PALETTE_TWO) {
  group_var <- if (is.null(group_col)) paste0(gene, "_group") else group_col
  df <- data.frame(time = clinical_df$paper_Survival..months.,
                   status = clinical_df$paper_Vital.status..1.dead.,
                   group = clinical_df[[group_var]])
  df <- df[stats::complete.cases(df), ]
  fit <- survival::survfit(survival::Surv(time, status) ~ group, data = df)
  survminer::ggsurvplot(fit, data = df, pval = TRUE, conf.int = FALSE,
                        censor.shape = 124, censor.size = 3, size = 1.1,
                        palette = palette, xlab = "Time (months)",
                        ylab = "Overall survival probability",
                        ggtheme = ggplot2::theme_classic(base_size = 12), title = gene)
}
