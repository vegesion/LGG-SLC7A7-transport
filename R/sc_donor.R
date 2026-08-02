# =============================================================================
# R/sc_donor.R  ---  donor 수준 단일세포 분석 / donor-level scRNA (main figures)
# -----------------------------------------------------------------------------
# ★ 설계 원칙
#  1) 세포 단위 +/- 이분은 main figure 에서 쓰지 않는다(pseudoreplication + depth 교란).
#  2) gene set 점수는 AddModuleScore 대신 AUCell(순위 기반) → depth 강건.
#  3) donor 수준 회귀에 depth(nFeature) 를 공변량으로.
#  4) ICC(dataset) 가 높으므로 dataset 을 랜덤효과로 넣고, within-dataset 메타분석으로 검증.
# =============================================================================

donor_pseudobulk <- function(obj, genes = NULL, layer = c("data", "counts"),
                             donor_col = "donor_id", min_cells = SC_DONOR_MIN_CELLS) {
  layer <- match.arg(layer)
  m <- Seurat::GetAssayData(obj, assay = "RNA", layer = layer)
  if (!is.null(genes)) m <- m[intersect(genes, rownames(m)), , drop = FALSE]
  d   <- factor(obj[[donor_col]][, 1])
  ind <- Matrix::sparse.model.matrix(~ 0 + d); colnames(ind) <- levels(d)
  n   <- Matrix::colSums(ind); agg <- as.matrix(m %*% ind)
  if (layer == "data") agg <- sweep(agg, 2, n, "/")
  agg[, n >= min_cells, drop = FALSE]
}

# ── UCell: depth 강건 + 메모리 효율적인 gene set 활성도 ──────────────────────
# AUCell 은 유전자×세포 dense 랭킹 행렬을 한 번에 만들어 대형 atlas 에서 RAM 이 폭발한다.
# UCell 은 순위 기반(Mann-Whitney U)으로 동일한 depth 강건성을 유지하면서 청크 처리한다.
# BPCells 온디스크 행렬을 통째로 실체화하지 않도록 세포를 직접 청크로 잘라 넘긴다.
score_auc <- function(obj, gene_sets, seed = AUCELL_SEED, n_cores = 1,
                      chunk_size = 2000, max_rank = 1500) {
  stopifnot(requireNamespace("UCell", quietly = TRUE))
  set.seed(seed)
  cells <- colnames(obj)
  cnt   <- Seurat::GetAssayData(obj, assay = "RNA", layer = "counts")  # lazy handle
  idx   <- split(seq_along(cells), ceiling(seq_along(cells) / chunk_size))
  message("UCell: ", length(cells), " cells / ", length(idx), " chunks")
  
  res <- lapply(seq_along(idx), function(i) {
    if (i %% 20 == 0) message("  chunk ", i, "/", length(idx))
    m <- as(cnt[, idx[[i]], drop = FALSE], "dgCMatrix")   # 청크만 메모리로
    s <- UCell::ScoreSignatures_UCell(m, features = gene_sets,
                                      maxRank = max_rank, ncores = n_cores)
    rm(m); gc(verbose = FALSE)
    s
  })
  out <- do.call(rbind, res)
  colnames(out) <- sub("_UCell$", "", colnames(out))      # downstream 이름 유지
  out[cells, , drop = FALSE]}

# 배치(dataset) 컬럼 자동 탐색
find_batch_col <- function(obj, candidates = SC_BATCH_CANDIDATES) {
  hit <- intersect(candidates, colnames(obj@meta.data))
  if (!length(hit)) NA_character_ else hit[1]
}

# donor 수준 집계표 (outcome/gene/depth/dataset)
donor_table <- function(obj, outcome, gene = GENE_OF_INTEREST, donor_col = "donor_id",
                        depth_col = DEPTH_COVARIATE) {
  bc <- find_batch_col(obj)
  df <- data.frame(
    outcome = if (length(outcome) == 1 && is.character(outcome)) obj[[outcome]][, 1] else outcome,
    gene  = Seurat::FetchData(obj, vars = gene)[, 1],
    depth = obj[[depth_col]][, 1],
    donor = obj[[donor_col]][, 1],
    dataset = if (is.na(bc)) "all" else as.character(obj[[bc]][, 1]), stringsAsFactors = FALSE)
  df %>% dplyr::group_by(donor, dataset) %>%
    dplyr::summarise(outcome = mean(outcome, na.rm = TRUE), gene = mean(gene, na.rm = TRUE),
                     depth = mean(depth, na.rm = TRUE), n = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(n >= SC_DONOR_MIN_CELLS)
}

# ── donor 수준 연속 검정: depth 보정 + dataset 랜덤효과 + within-dataset 메타 ─
donor_continuous_test <- function(obj, outcome, gene = GENE_OF_INTEREST,
                                  adjust_depth = TRUE, mixed = TRUE) {
  agg <- donor_table(obj, outcome, gene)
  form <- if (adjust_depth) outcome ~ gene + depth else outcome ~ gene

  fit_lm <- stats::lm(form, data = agg)

  fit_mm <- NULL
  if (mixed && dplyr::n_distinct(agg$dataset) > 2 && requireNamespace("lme4", quietly = TRUE)) {
    f <- if (adjust_depth) outcome ~ gene + depth + (1 | dataset) else outcome ~ gene + (1 | dataset)

    fit_mm <- try(if (requireNamespace("lmerTest", quietly = TRUE))
                    lmerTest::lmer(f, data = agg) else lme4::lmer(f, data = agg), silent = TRUE)
    if (inherits(fit_mm, "try-error")) fit_mm <- NULL
  }
  list(data = agg, lm = fit_lm, coef = summary(fit_lm)$coefficients,
       mixed = fit_mm, mixed_coef = if (is.null(fit_mm)) NULL else summary(fit_mm)$coefficients,
       meta = meta_within_dataset(agg, adjust_depth),
       spearman = stats::cor.test(agg$gene, agg$outcome, method = "spearman"))}


# within-dataset 추정 → 역분산 가중 메타분석 (배치 지배적일 때의 방어적 근거)
meta_within_dataset <- function(agg, adjust_depth = TRUE, min_n = 4) {
  per <- agg %>% dplyr::group_by(dataset) %>% dplyr::filter(dplyr::n() >= min_n) %>%
    dplyr::group_modify(~{
      f <- if (adjust_depth) stats::lm(outcome ~ gene + depth, data = .x)
           else              stats::lm(outcome ~ gene, data = .x)
      s <- summary(f)$coefficients
      if (!"gene" %in% rownames(s)) return(data.frame())
      data.frame(beta = s["gene", 1], se = s["gene", 2], n = nrow(.x))
    }) %>% dplyr::ungroup() %>% dplyr::filter(is.finite(beta), is.finite(se), se > 0)
  if (!nrow(per)) return(list(per_dataset = per, pooled = NA, se = NA, p = NA))
  w <- 1 / per$se^2; b <- sum(w * per$beta) / sum(w); se <- sqrt(1 / sum(w))
  list(per_dataset = per, pooled = b, se = se, p = 2 * stats::pnorm(-abs(b / se)),
       n_dataset = nrow(per))
}

# ── donor 수준: 발현 vs 세포 아종 비율 회귀 (Fig 6 핵심) ────────────────────
donor_subtype_proportion <- function(obj, subtype_col, gene = GENE_OF_INTEREST,
                                     donor_col = "donor_id", adjust_depth = TRUE) {
  bc <- find_batch_col(obj)
  meta <- data.frame(donor = obj[[donor_col]][, 1],
                     subtype = as.character(obj[[subtype_col]][, 1]),
                     gene = Seurat::FetchData(obj, vars = gene)[, 1],
                     depth = obj[[DEPTH_COVARIATE]][, 1],
                     dataset = if (is.na(bc)) "all" else as.character(obj[[bc]][, 1]),
                     stringsAsFactors = FALSE)
  prop <- meta %>% dplyr::count(donor, dataset, subtype) %>%
    dplyr::group_by(donor) %>% dplyr::mutate(total = sum(n), prop = n / total) %>%
    dplyr::ungroup() %>% dplyr::filter(total >= SC_DONOR_MIN_CELLS)
  gexp <- meta %>% dplyr::group_by(donor) %>%
    dplyr::summarise(gene = mean(gene), depth = mean(depth), .groups = "drop")
  d <- dplyr::left_join(prop, gexp, by = "donor")

  res <- d %>% dplyr::group_split(subtype) %>% lapply(function(x) {
    if (nrow(x) < 10 || stats::sd(x$prop) == 0) return(NULL)
    f <- if (adjust_depth) stats::lm(prop ~ gene + depth, data = x) else stats::lm(prop ~ gene, data = x)
    s <- summary(f)$coefficients
    mm <- if (dplyr::n_distinct(x$dataset) > 2 && requireNamespace("lme4", quietly = TRUE))
      try(lme4::lmer(prop ~ gene + depth + (1 | dataset), data = x), silent = TRUE) else NULL
    data.frame(subtype = x$subtype[1], beta = s["gene", 1], se = s["gene", 2],
               p = s["gene", 4], n_donor = nrow(x),
               beta_mixed = if (is.null(mm) || inherits(mm, "try-error")) NA
                            else lme4::fixef(mm)["gene"])
  }) %>% dplyr::bind_rows() %>%
    dplyr::mutate(FDR = stats::p.adjust(p, "BH")) %>% dplyr::arrange(p)
  list(data = d, result = res)
}

# ── depth 민감도 3종 (supplementary) ────────────────────────────────────────
depth_sensitivity <- function(obj, module_col, gene = GENE_OF_INTEREST,
                              donor_col = "donor_id", depth_col = DEPTH_COVARIATE,
                              n_null = 20, seed = 42) {
  status <- ifelse(Seurat::FetchData(obj, vars = gene)[, 1] > 0, "Positive", "Negative")
  depth  <- obj[[depth_col]][, 1]; donor <- obj[[donor_col]][, 1]; score <- obj[[module_col]][, 1]
  ratio  <- mean(depth[status == "Positive"]) / mean(depth[status == "Negative"])

  glmer_res <- NULL
  if (requireNamespace("lme4", quietly = TRUE)) {
    dd <- data.frame(pos = as.integer(status == "Positive"),
                     depth_s = as.numeric(scale(depth)), donor = donor)
    m <- try(lme4::glmer(pos ~ depth_s + (1 | donor), data = dd, family = "binomial"), silent = TRUE)
    if (!inherits(m, "try-error"))
      glmer_res <- data.frame(OR = exp(lme4::fixef(m)["depth_s"]),
                              p = summary(m)$coefficients["depth_s", 4], row.names = NULL)
  }
  set.seed(seed)
  cnt <- Seurat::GetAssayData(obj, assay = "RNA", layer = "counts")
  samp <- sample(rownames(cnt), min(3000, nrow(cnt)))
  det <- Matrix::rowMeans(as(cnt[samp, ], "dgCMatrix") > 0)
  target <- mean(Seurat::FetchData(obj, vars = gene)[, 1] > 0)
  pool <- names(sort(abs(det - target)))[seq_len(min(n_null, length(det)))]
  dlt <- function(st) { a <- tapply(score, list(donor, st), mean)
    if (!all(c("Positive", "Negative") %in% colnames(a))) return(NA_real_)
    stats::median(a[, "Positive"] - a[, "Negative"], na.rm = TRUE) }
  null_delta <- vapply(pool, function(g)
    dlt(ifelse(Seurat::FetchData(obj, vars = g)[, 1] > 0, "Positive", "Negative")), numeric(1))
  obs <- dlt(status)
  list(depth_ratio = ratio, glmer = glmer_res, observed_delta = obs, null_delta = null_delta,
       empirical_p = mean(abs(stats::na.omit(null_delta)) >= abs(obs)))
}

compare_mg_mdm <- function(obj, gene = GENE_OF_INTEREST, celltype_col = "cell_type",
                           mg = "microglial cell", mdm = c("macrophage", "monocyte")) {
  df <- data.frame(expr = Seurat::FetchData(obj, vars = gene)[, 1],
                   ct = as.character(obj[[celltype_col]][, 1]), donor = obj$donor_id)
  df$lineage <- dplyr::case_when(df$ct == mg ~ "Microglia", df$ct %in% mdm ~ "MDM", TRUE ~ NA_character_)
  agg <- df %>% dplyr::filter(!is.na(lineage)) %>%
    dplyr::group_by(donor, lineage) %>% dplyr::summarise(m = mean(expr), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = lineage, values_from = m) %>%
    dplyr::filter(!is.na(Microglia), !is.na(MDM))
  list(data = agg, test = stats::wilcox.test(agg$Microglia, agg$MDM, paired = TRUE),
       delta_median = stats::median(agg$Microglia - agg$MDM))
}
