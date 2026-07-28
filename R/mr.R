# =============================================================================
# R/mr.R  ---  Mendelian Randomization 함수 (exploratory) / cis-MR + MRcML
# -----------------------------------------------------------------------------
# 원본 "MR통합 찐 최종 코드.R" (가장 완성도 높은 버전)을 정리한 것.
# exposure = eQTLGen cis-eQTL(OpenGWAS eqtl-a-*), outcome = 뇌종양/GBM GWAS.
# ※ 이 분석은 결과의 유의성·일관성이 낮아 논문에는 포함하지 않았습니다(탐색용 보존).
# =============================================================================

# ── biomaRt 연결 (mirror fallback) ───────────────────────────────────────────
connect_mart <- function(mirrors = c("www", "useast", "uswest", "asia")) {
  for (m in mirrors) {
    message(sprintf("biomaRt 연결 시도: %s", m))
    mart <- tryCatch(
      biomaRt::useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = m),
      error = function(e) NULL)
    if (!is.null(mart)) { message("연결 성공: ", m); return(mart) }
  }
  message("모든 mirror 실패"); NULL
}

# ── ENSEMBL id → 유전자 위치 (biomaRt, 실패 시 REST) ─────────────────────────
get_gene_location <- function(ensembl_id, mart = NULL) {
  if (!is.null(mart)) {
    info <- tryCatch(biomaRt::getBM(
      attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position"),
      filters = "ensembl_gene_id", values = ensembl_id, mart = mart),
      error = function(e) NULL)
    if (!is.null(info) && nrow(info) > 0) return(info)
  }
  tryCatch({
    url  <- sprintf("https://rest.ensembl.org/lookup/id/%s?content-type=application/json", ensembl_id)
    resp <- jsonlite::fromJSON(url)
    data.frame(ensembl_gene_id = ensembl_id, chromosome_name = resp$seq_region_name,
               start_position = resp$start, end_position = resp$end)
  }, error = function(e) NULL)
}

# ── cis-SNP 필터 (유전자 ±window) / keep only cis instruments ────────────────
filter_cis_snps <- function(exposure_dat, ensembl_id, mart = NULL, cis_window = MR_CIS_WINDOW) {
  if (!all(c("chr.exposure", "pos.exposure") %in% names(exposure_dat))) return(exposure_dat)
  gi <- get_gene_location(ensembl_id, mart = mart)
  if (is.null(gi) || nrow(gi) == 0) return(exposure_dat)
  snp_chr <- gsub("chr", "", as.character(exposure_dat$chr.exposure))
  mask <- snp_chr == as.character(gi$chromosome_name[1]) &
    exposure_dat$pos.exposure >= (gi$start_position[1] - cis_window) &
    exposure_dat$pos.exposure <= (gi$end_position[1]   + cis_window)
  exposure_dat[mask, ]
}

# ── 유전자 1개 TwoSampleMR + 민감도 분석 / one-gene MR with sensitivity ──────
# 반환 list: summary, harmonised, heterogeneity, steiger, pleiotropy, presso, loo
run_mr_for_gene <- function(gene_id, outcomes, outcome_n, mart = NULL) {
  exposureID <- paste0("eqtl-a-", gene_id)
  exposure_dat <- tryCatch(TwoSampleMR::extract_instruments(
    outcomes = exposureID, p1 = MR_INSTRUMENT_P, clump = TRUE, r2 = MR_CLUMP_R2),
    error = function(e) NULL)
  if (is.null(exposure_dat) || nrow(exposure_dat) == 0) return(NULL)

  exposure_dat <- filter_cis_snps(exposure_dat, gene_id, mart = mart)
  if (nrow(exposure_dat) == 0) return(NULL)

  outcome_dat <- tryCatch(TwoSampleMR::extract_outcome_data(
    snps = exposure_dat$SNP, outcomes = outcomes), error = function(e) NULL)
  if (is.null(outcome_dat) || nrow(outcome_dat) == 0) return(NULL)

  dat_raw <- tryCatch(TwoSampleMR::harmonise_data(exposure_dat, outcome_dat), error = function(e) NULL)
  if (is.null(dat_raw) || nrow(dat_raw) == 0) return(NULL)

  dat <- dat_raw[dat_raw$mr_keep == TRUE, ]
  if (nrow(dat) == 0) return(NULL)

  dat$Fstat <- (dat$beta.exposure / dat$se.exposure)^2  # weak IV 제거
  dat <- dat[dat$Fstat >= MR_FSTAT_MIN, ]
  if (nrow(dat) == 0) return(NULL)
  n_snp <- nrow(dat)

  method_list <- if (n_snp == 1) "mr_wald_ratio"
                 else if (n_snp == 2) "mr_ivw"
                 else c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode")
  res <- tryCatch(TwoSampleMR::mr(dat, method_list = method_list), error = function(e) NULL)
  if (is.null(res) || nrow(res) == 0) return(NULL)
  res$gene <- gene_id; res$nsnp <- n_snp; res$mean_Fstat <- mean(dat$Fstat)

  het <- if (n_snp >= 2) tryCatch(TwoSampleMR::mr_heterogeneity(dat), error = function(e) NULL) else NULL
  if (!is.null(het)) het$gene <- gene_id

  steiger <- tryCatch({
    dat$samplesize.outcome <- outcome_n
    if (all(is.na(dat$samplesize.exposure))) dat$samplesize.exposure <- MR_EQTL_EXPOSURE_N
    st <- TwoSampleMR::steiger_filtering(dat)
    data.frame(gene = gene_id, nsnp = n_snp,
               n_correct_dir = sum(st$steiger_dir == TRUE, na.rm = TRUE),
               n_wrong_dir   = sum(st$steiger_dir == FALSE, na.rm = TRUE),
               n_unknown_dir = sum(is.na(st$steiger_dir)))
  }, error = function(e) NULL)

  plt <- if (n_snp >= 3) tryCatch(TwoSampleMR::mr_pleiotropy_test(dat), error = function(e) NULL) else NULL
  if (!is.null(plt)) plt$gene <- gene_id

  presso_summary <- NULL
  if (n_snp >= 4) {
    presso <- tryCatch(MRPRESSO::mr_presso(
      BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
      SdOutcome = "se.outcome", SdExposure = "se.exposure",
      OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat,
      NbDistribution = 1000, SignifThreshold = 0.05), error = function(e) NULL)
    if (!is.null(presso)) presso_summary <- tryCatch({
      main <- presso$`Main MR results`
      data.frame(gene = gene_id, nsnp = n_snp,
        presso_global_p  = presso$`MR-PRESSO results`$`Global Test`$Pvalue,
        presso_raw_b     = main$`Causal Estimate`[1], presso_raw_p = main$`P-value`[1],
        presso_correct_b = ifelse(nrow(main) > 1, main$`Causal Estimate`[2], NA),
        presso_correct_p = ifelse(nrow(main) > 1, main$`P-value`[2], NA))
    }, error = function(e) NULL)
  }

  loo <- if (n_snp >= 2) tryCatch(TwoSampleMR::mr_leaveoneout(dat), error = function(e) NULL) else NULL
  if (!is.null(loo)) loo$gene <- gene_id

  list(summary = res, harmonised = dat, heterogeneity = het,
       steiger = steiger, pleiotropy = plt, presso_summary = presso_summary, leaveoneout = loo)
}

# ── 유전자 1개 MRcML (BIC/AIC + DP) / constrained ML MR ──────────────────────
run_cml_for_gene <- function(gene_id, outcomes, n_sample, mart = NULL, rho = NULL) {
  exposureID <- paste0("eqtl-a-", gene_id)
  exposure_raw <- tryCatch(TwoSampleMR::extract_instruments(
    outcomes = exposureID, p1 = MR_INSTRUMENT_P, clump = TRUE, r2 = MR_CLUMP_R2),
    error = function(e) NULL)
  if (is.null(exposure_raw) || nrow(exposure_raw) == 0) return(NULL)
  exposure_dat <- filter_cis_snps(exposure_raw, gene_id, mart = mart)
  if (nrow(exposure_dat) == 0) return(NULL)
  outcome_dat <- tryCatch(TwoSampleMR::extract_outcome_data(
    snps = exposure_dat$SNP, outcomes = outcomes), error = function(e) NULL)
  if (is.null(outcome_dat) || nrow(outcome_dat) == 0) return(NULL)
  dat <- tryCatch(TwoSampleMR::harmonise_data(exposure_dat, outcome_dat), error = function(e) NULL)
  if (is.null(dat) || nrow(dat) == 0) return(NULL)
  dat <- dat[dat$mr_keep == TRUE, ]
  if (nrow(dat) == 0) return(NULL)
  dat$Fstat <- (dat$beta.exposure / dat$se.exposure)^2
  dat <- dat[dat$Fstat >= MR_FSTAT_MIN, ]
  if (nrow(dat) < 3) return(NULL)   # MRcML은 최소 3 SNP

  bx <- dat$beta.exposure; bxse <- dat$se.exposure
  by <- dat$beta.outcome;  byse <- dat$se.outcome
  extract_cml <- function(label, r) tryCatch(data.frame(
    gene = gene_id, method = label, nsnp = nrow(dat), mean_Fstat = mean(dat$Fstat),
    theta_BIC = r$MA_BIC_theta, se_BIC = r$MA_BIC_se, pval_BIC = r$MA_BIC_p,
    n_invalid_BIC = length(r$BIC_invalid),
    theta_AIC = r$MA_AIC_theta, se_AIC = r$MA_AIC_se, pval_AIC = r$MA_AIC_p,
    n_invalid_AIC = length(r$AIC_invalid)), error = function(e) NULL)

  rows <- list()
  rows[["mr_cML"]]    <- extract_cml("mr_cML",
    tryCatch(MRcML::mr_cML(bx, by, bxse, byse, n = n_sample, random_start = 100), error = function(e) NULL))
  rows[["mr_cML_DP"]] <- extract_cml("mr_cML_DP",
    tryCatch(MRcML::mr_cML_DP(bx, by, bxse, byse, n = n_sample, random_start = 10,
             random_start_pert = 10, random_seed = 1, num_pert = 200), error = function(e) NULL))
  rows <- rows[!sapply(rows, is.null)]
  if (length(rows) == 0) return(NULL)
  list(summary = do.call(rbind, rows), harmonised = dat)
}

# ── 진행 막대 로그 / simple progress bar message ─────────────────────────────
.mr_progress <- function(i, total, gene_id) {
  pct <- round(i / total * 100)
  message(sprintf("[%d/%d] %s%s %3d%% | %s", i, total,
    strrep("#", round(pct / 5)), strrep("-", 20 - round(pct / 5)), pct, gene_id))
}

# ── 여러 유전자 배치 실행 / batch over a gene list ───────────────────────────
run_mr_batch <- function(ensembl_list, outcomes, outcome_n, mart = NULL) {
  total <- length(ensembl_list); out <- vector("list", total)
  for (i in seq_along(ensembl_list)) {
    .mr_progress(i, total, ensembl_list[i])
    out[[i]] <- run_mr_for_gene(ensembl_list[i], outcomes, outcome_n, mart = mart)
  }
  names(out) <- ensembl_list
  out[!sapply(out, is.null)]
}

run_cml_batch <- function(ensembl_list, outcomes, n_sample, mart = NULL) {
  total <- length(ensembl_list); out <- vector("list", total)
  for (i in seq_along(ensembl_list)) {
    .mr_progress(i, total, ensembl_list[i])
    out[[i]] <- run_cml_for_gene(ensembl_list[i], outcomes, n_sample, mart = mart)
  }
  names(out) <- ensembl_list
  out[!sapply(out, is.null)]
}
