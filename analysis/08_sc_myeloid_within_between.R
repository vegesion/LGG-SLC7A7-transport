# =============================================================================
# 08_sc_myeloid_within_between.R            [Fig 6 — 08번 대체본]
# 세포 단위 within-between 분해 혼합모형 (Mundlak / REWB)
# -----------------------------------------------------------------------------
#   y ~ SLC7A7_within + SLC7A7_between + depth_within + depth_between
#         + (1 | donor) + (1 | dataset)
#
#   SLC7A7_between = donor 평균           → "SLC7A7 높은 환자"의 효과 (환자 간)
#   SLC7A7_within  = 세포값 − donor 평균  → "같은 환자 안에서 SLC7A7 높은 세포"의 효과
#   (1|donor)  : pseudoreplication 제거   (1|dataset) : 배치(ICC 0.56) 흡수
#   depth 도 동일하게 분해 → between 효과가 donor 간 depth 차이와 교란되지 않음
#
# 왜 pseudobulk 대신 이 모형인가
#   · pseudobulk 는 donor 안의 세포 이질성을 평균으로 지워버림
#   · 세포 단위 회귀는 pseudoreplication 으로 p 가 붕괴됨
#   · within-between 은 두 층위를 분리해 각각 올바른 자유도로 검정
#   · 추가로 contextual effect (between − within) 를 명시적으로 검정 가능
#
# 입력 : data/processed/{microglia,myeloid}.rds
# 출력 : results/Fig6_withinbetween_modules.csv, ..._genes.csv,
#        figures/Fig6_withinbetween_*.pdf
# =============================================================================

source(here::here("R", "setup.R"))
suppressPackageStartupMessages({
  library(Seurat); library(lme4); library(ggplot2); library(dplyr); library(tidyr)
})
if (requireNamespace("lmerTest", quietly = TRUE)) library(lmerTest)   # p값 (Satterthwaite)

## ── 0. 설정 ─────────────────────────────────────────────────────────────────
CELLS_PER_DONOR <- NA      # donor 당 최대 세포 (계산량 제어; NA 면 전체 사용)
MIN_CELLS_DONOR <- 20       # donor 최소 세포 수
GENE_MIN_PCT    <- 0.10     # 유전자 수준: 최소 검출 세포 비율
RUN_GENE_LEVEL  <- TRUE     # 유전자 수준까지 돌릴지 (오래 걸림)
SEED            <- 42

mg <- readRDS(file.path(DIR_DATA_PROC, "microglia.rds"))
mg <- label_slc_status(mg)                       # (호환용, 이 스크립트에선 미사용)

## ── 1. 세포 단위 설계행렬: within-between 분해 ──────────────────────────────
batch_col <- find_batch_col(mg)
if (is.na(batch_col)) stop("dataset(배치) 컬럼을 찾지 못했습니다.")

meta <- data.frame(
  cell    = colnames(mg),
  g       = Seurat::FetchData(mg, vars = GENE_OF_INTEREST)[, 1],
  depth   = mg[[DEPTH_COVARIATE]][, 1],
  donor   = as.character(mg$donor_id),
  dataset = as.character(mg[[batch_col]][, 1]),
  stringsAsFactors = FALSE
)

# donor 최소 세포 수 필터 + donor 당 서브샘플
set.seed(SEED)
keep_donor <- names(which(table(meta$donor) >= MIN_CELLS_DONOR))
meta <- meta[meta$donor %in% keep_donor, ]
if (!is.na(CELLS_PER_DONOR)) {
  meta <- meta %>% group_by(donor) %>%
    slice_sample(n = CELLS_PER_DONOR) %>% ungroup() %>% as.data.frame()
}
message(sprintf("[설계] 세포 %d, donor %d, dataset %d",
                nrow(meta), dplyr::n_distinct(meta$donor), dplyr::n_distinct(meta$dataset)))

# ★ within-between 분해 (공통 스케일 유지 → 두 계수를 직접 비교 가능)
decompose <- function(x, grp) {
  bmean <- ave(x, grp, FUN = function(v) mean(v, na.rm = TRUE))
  s <- stats::sd(x, na.rm = TRUE); if (!is.finite(s) || s == 0) s <- 1
  list(within = (x - bmean) / s, between = (bmean - mean(x, na.rm = TRUE)) / s)
}

dg <- decompose(meta$g,     meta$donor)
dd <- decompose(meta$depth, meta$donor)
meta$g_within <- dg$within;  meta$g_between <- dg$between
meta$d_within <- dd$within;  meta$d_between <- dd$between
meta$donor <- factor(meta$donor); meta$dataset <- factor(meta$dataset)

# 분해 점검: 두 성분의 분산 비율 (between 이 너무 작으면 환자 간 검정력 부족)
cat(sprintf("\n[분해] var(within) = %.3f, var(between) = %.3f  → between 비중 %.1f%%\n",
            var(meta$g_within), var(meta$g_between),
            100 * var(meta$g_between) / (var(meta$g_within) + var(meta$g_between))))

## ── 2. 적합 함수 ────────────────────────────────────────────────────────────
FORM <- y ~ g_within + g_between + d_within + d_between + (1 | donor) + (1 | dataset)

fit_wb <- function(y, df = meta, label = "") {
  df$y <- as.numeric(y)
  if (stats::sd(df$y, na.rm = TRUE) == 0) return(NULL)
  m <- try(suppressMessages(lmerTest::lmer(FORM, data = df,
             control = lme4::lmerControl(optimizer = "bobyqa",
                                         optCtrl = list(maxfun = 2e5)))), silent = TRUE)
  if (inherits(m, "try-error")) return(NULL)
  cf <- summary(m)$coefficients
  gp <- function(term, col) if (term %in% rownames(cf)) cf[term, col] else NA_real_
  pcol <- if ("Pr(>|t|)" %in% colnames(cf)) "Pr(>|t|)" else NA   # lmerTest 유무
  # contextual effect = between − within (두 층위가 다른가)
  V  <- as.matrix(vcov(m))
  ctx <- gp("g_between", "Estimate") - gp("g_within", "Estimate")
  se_ctx <- sqrt(V["g_between","g_between"] + V["g_within","g_within"] -
                   2*V["g_between","g_within"])
  data.frame(
    response   = label,
    beta_within  = gp("g_within",  "Estimate"), se_within  = gp("g_within",  "Std. Error"),
    p_within   = if (is.na(pcol)) NA_real_ else gp("g_within",  pcol),
    beta_between = gp("g_between", "Estimate"), se_between = gp("g_between", "Std. Error"),
    p_between  = if (is.na(pcol)) NA_real_ else gp("g_between", pcol),
    beta_depth_w = gp("d_within",  "Estimate"), beta_depth_b = gp("d_between", "Estimate"),
    contextual = ctx, se_contextual = se_ctx,
    p_contextual = 2 * stats::pnorm(-abs(ctx / se_ctx)),
    n_cell = nrow(df), n_donor = dplyr::n_distinct(df$donor),
    stringsAsFactors = FALSE)
}

## ── 3. 모듈 수준 (주 분석) ──────────────────────────────────────────────────
tr <- build_transporter_geneset(universe = rownames(mg))
MODULES <- list(
  TRANSPORT   = setdiff(tr, GENE_OF_INTEREST),          # ★ self 제외 (순환논리 차단)
  ARG_ENZYME  = build_arg_enzyme_geneset(universe = rownames(mg)),
  GLYCOLYSIS  = get_msig_genes("HALLMARK_GLYCOLYSIS"),
  OXPHOS      = get_msig_genes("HALLMARK_OXIDATIVE_PHOSPHORYLATION"),
  ROS         = get_msig_genes("HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"),
  HYPOXIA     = get_msig_genes("HALLMARK_HYPOXIA"),
  NFKB        = get_msig_genes("HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
  INFLAMMATORY= get_msig_genes("HALLMARK_INFLAMMATORY_RESPONSE"),
  IFN_GAMMA   = get_msig_genes("HALLMARK_INTERFERON_GAMMA_RESPONSE"),
  IFN_ALPHA   = get_msig_genes("HALLMARK_INTERFERON_ALPHA_RESPONSE"),
  FAO         = get_msig_genes("HALLMARK_FATTY_ACID_METABOLISM"),
  IFN_REACTOME = msig_get("REACTOME_INTERFERON_ALPHA_BETA_SIGNALING", "C2", "CP:REACTOME"),
  PPP_REACTOME = msig_get("REACTOME_PENTOSE_PHOSPHATE_PATHWAY",        "C2", "CP:REACTOME"),
  POLYAMINE    = msig_get("GOBP_POLYAMINE_BIOSYNTHETIC_PROCESS",       "C5", "GO:BP"))
MODULES <- lapply(MODULES, function(g) setdiff(intersect(g, rownames(mg)), GENE_OF_INTEREST))
MODULES <- MODULES[vapply(MODULES, length, 1L) >= 5]

message("[모듈] UCell 점수 계산 (", length(MODULES), "개)")
mg_sub <- subset(mg, cells = meta$cell)
auc <- score_auc(mg_sub, MODULES)                 # UCell (depth 강건, 청크 처리)
auc <- auc[meta$cell, , drop = FALSE]

res_mod <- dplyr::bind_rows(lapply(colnames(auc), function(nm) {
  message("  ", nm); fit_wb(auc[, nm], label = nm)
}))
res_mod <- res_mod %>%
  mutate(fdr_within  = p.adjust(p_within,  "BH"),
         fdr_between = p.adjust(p_between, "BH"),
         fdr_ctx     = p.adjust(p_contextual, "BH"))
write_result(res_mod, "Fig6_withinbetween_modules.csv")
cat("\n===== 모듈 수준 within-between =====\n"); print(as.data.frame(res_mod), digits = 3)

## ── 4. 유전자 수준 (선택) ───────────────────────────────────────────────────
if (RUN_GENE_LEVEL) {
  dat_l <- Seurat::GetAssayData(mg_sub, assay = "RNA", layer = "data")
  dat_l <- as(dat_l[, meta$cell, drop = FALSE], "dgCMatrix")
  pct   <- Matrix::rowMeans(dat_l > 0)
  genes <- setdiff(rownames(dat_l)[pct >= GENE_MIN_PCT], GENE_OF_INTEREST)
  message("[유전자] ", length(genes), "개 × ", nrow(meta), " 세포 — 시간이 걸립니다")

  res_gene <- dplyr::bind_rows(pbapply::pblapply(genes, function(g)
    fit_wb(as.numeric(dat_l[g, ]), label = g)))
  res_gene <- res_gene %>%
    mutate(fdr_within = p.adjust(p_within, "BH"),
           fdr_between = p.adjust(p_between, "BH"))
  write_result(res_gene, "Fig6_withinbetween_genes.csv")
  cat("\n유의 유전자 — within:", sum(res_gene$fdr_within < 0.05, na.rm = TRUE),
      "/ between:", sum(res_gene$fdr_between < 0.05, na.rm = TRUE), "\n")
}

## ── 5. 시각화 ───────────────────────────────────────────────────────────────
# (a) 모듈별 within vs between 계수 비교
pl <- res_mod %>%
  dplyr::select(response, beta_within, se_within, beta_between, se_between) %>%
  pivot_longer(-response, names_to = c(".value", "level"), names_pattern = "(beta|se)_(.*)") %>%
  mutate(level = factor(level, levels = c("within", "between"),
                        labels = c("Within-donor (cell)", "Between-donor (patient)")),
         lo = beta - 1.96*se, hi = beta + 1.96*se,
         y_base = match(response, res_mod$response[order(res_mod$beta_between)]),
         y = y_base + ifelse(level == "Within-donor (cell)", -0.18, 0.18))

pA <- ggplot(pl, aes(beta, y, color = level)) +
  geom_vline(xintercept = 0, color = "grey55") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.7) +
  geom_point(size = 2.2) +
  scale_color_manual(values = c("Within-donor (cell)" = "#7D6608",
                                "Between-donor (patient)" = "#C0392B"), name = NULL) +
  scale_y_continuous(breaks = seq_along(unique(pl$y_base)),
                     labels = res_mod$response[order(res_mod$beta_between)]) +
  labs(x = expression(beta~"(module score per 1 SD "*italic("SLC7A7")*")"), y = NULL,
       title = "Within- vs between-donor effects of SLC7A7") +
  theme_paper(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(DIR_FIGURES, "Fig6_withinbetween_modules.pdf"), pA,
       width = 6.6, height = 0.42*nrow(res_mod) + 2)

# (b) contextual effect (between − within): 두 층위가 다른가
pB <- res_mod %>%
  mutate(lo = contextual - 1.96*se_contextual, hi = contextual + 1.96*se_contextual,
         response = factor(response, levels = response[order(contextual)])) %>%
  ggplot(aes(contextual, response)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey45") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.7) +
  geom_point(aes(shape = fdr_ctx < 0.05), size = 2.5, fill = "white", stroke = 0.8) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), guide = "none") +
  labs(x = "contextual effect (between − within)", y = NULL,
       title = "Do patient-level and cell-level effects differ?",
       caption = "0 에서 벗어나면 두 층위의 연관이 서로 다름 (단순 세포 단위 회귀가 오도될 수 있음)") +
  theme_paper(base_size = 11)
ggsave(file.path(DIR_FIGURES, "Fig6_withinbetween_contextual.pdf"), pB,
       width = 6.2, height = 0.42*nrow(res_mod) + 2)

message("완료: within-between 혼합모형")

# -----------------------------------------------------------------------------
# 해석 가이드
#  · beta_between 유의  → "SLC7A7 이 높은 환자"의 myeloid 가 그 모듈이 높다 (환자 층위)
#                          = 기존 donor 수준 분석과 같은 질문, 더 올바른 자유도
#  · beta_within  유의  → "같은 환자 안에서 SLC7A7 이 높은 세포"가 그 모듈이 높다
#                          = 세포 내재적 공발현. 단, depth 교란에 가장 취약하므로
#                            d_within 계수가 크면 해석 주의
#  · contextual ≠ 0     → 두 층위의 연관이 다름. 세포 단위만 본 기존 문헌이
#                          오도될 수 있다는 근거 (본 연구의 방법론적 기여)
#
# 주의
#  · lmerTest 가 없으면 p 값이 NA 로 나옵니다 → install.packages("lmerTest")
#  · CELLS_PER_DONOR 로 세포를 줄여도 between 층위의 검정력(donor 수)은 그대로입니다.
#    within 층위 검정력만 줄어드니, 최종 실행은 값을 키우거나 NA(전체)로 두세요.
# -----------------------------------------------------------------------------
