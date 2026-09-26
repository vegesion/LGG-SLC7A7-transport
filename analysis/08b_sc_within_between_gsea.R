# =============================================================================
# 08b_sc_within_between_gsea.R          [Fig 6 — 세포 수준 within/between GSEA]
# -----------------------------------------------------------------------------
# 유전자마다 within-between 혼합모형을 적합하고, 그 t 통계량으로 GSEA 를 수행한다.
#
#   gene_ij ~ SLC7A7_within + SLC7A7_between + depth_within + depth_between
#               + (1 | donor) + (1 | dataset)
#
#   ranking(within)  = t(β_within)   → "같은 환자 안에서 SLC7A7 과 공발현하는 축"
#   ranking(between) = t(β_between)  → "SLC7A7 이 높은 환자의 myeloid 가 가진 축"
#
# ★ 두 가지 검정을 항상 함께 보고한다
#     fgsea    : 유전자 독립 가정 → 유의성이 부풀 수 있음
#     cameraPR : inter-gene correlation 보정 (preranked camera)
#   앞선 pseudobulk 분석에서 fgsea 24개 vs camera 0개였던 사례가 있으므로,
#   fgsea 단독 결과는 논문의 근거로 쓰지 않는다.
#
# 입력 : data/processed/microglia.rds
# 출력 : results/Fig6_wb_gene_stats.csv
#        results/Fig6_wb_gsea_{within,between}.csv
#        figures/Fig6_wb_gsea_*.pdf
# =============================================================================

source(here::here("R", "setup.R"))
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(lme4)
  library(fgsea); library(limma); library(msigdbr)
  library(dplyr); library(tidyr); library(ggplot2)
})

## ── 0. 설정 ─────────────────────────────────────────────────────────────────
CELLS_PER_DONOR <- NA      # donor 당 최대 세포 (NA = 전체). within 검정력에만 영향
MIN_CELLS_DONOR <- 20
GENE_MIN_PCT    <- 0.10     # 최소 검출 세포 비율 (0-팽창 완화)
INTER_GENE_COR  <- 0.01     # cameraPR 기본 가정 상관 (관측값이 있으면 아래에서 대체)
SEED            <- 42
N_CORES         <- 1        # >1 이면 parallel::mclapply (Windows 는 1 유지)

set.seed(SEED)
mg <- readRDS(file.path(DIR_DATA_PROC, "microglia.rds"))

## ── 1. 세포 수준 설계: within-between 분해 ──────────────────────────────────
batch_col <- find_batch_col(mg)
if (is.na(batch_col)) stop("dataset(배치) 컬럼을 찾지 못했습니다.")

meta <- data.frame(
  cell    = colnames(mg),
  g       = Seurat::FetchData(mg, vars = GENE_OF_INTEREST)[, 1],
  depth   = mg[[DEPTH_COVARIATE]][, 1],
  donor   = as.character(mg$donor_id),
  dataset = as.character(mg[[batch_col]][, 1]),
  stringsAsFactors = FALSE)

meta <- meta[meta$donor %in% names(which(table(meta$donor) >= MIN_CELLS_DONOR)), ]
if (!is.na(CELLS_PER_DONOR))
  meta <- meta %>% group_by(donor) %>%
    slice_sample(n = min(CELLS_PER_DONOR, dplyr::n())) %>% ungroup() %>% as.data.frame()

decompose <- function(x, grp) {
  bmean <- ave(x, grp, FUN = function(v) mean(v, na.rm = TRUE))
  s <- stats::sd(x, na.rm = TRUE); if (!is.finite(s) || s == 0) s <- 1
  list(within = (x - bmean) / s, between = (bmean - mean(x, na.rm = TRUE)) / s)
}
dg <- decompose(meta$g, meta$donor); dd <- decompose(meta$depth, meta$donor)
meta$g_within <- dg$within; meta$g_between <- dg$between
meta$d_within <- dd$within; meta$d_between <- dd$between
meta$donor <- factor(meta$donor); meta$dataset <- factor(meta$dataset)

N_DONOR <- nlevels(meta$donor)
message(sprintf("[설계] 세포 %s / donor %d / dataset %d",
                format(nrow(meta), big.mark = ","), N_DONOR, nlevels(meta$dataset)))
cat(sprintf("[분해] between 분산 비중 = %.1f%%\n",
            100 * var(meta$g_between) / (var(meta$g_within) + var(meta$g_between))))

## ── 2. 유전자 필터 + 세포×유전자 행렬 (열 추출이 빠르도록 전치) ─────────────
expr <- Seurat::GetAssayData(mg, assay = "RNA", layer = "data")
expr <- as(expr[, meta$cell, drop = FALSE], "dgCMatrix")
pct  <- Matrix::rowMeans(expr > 0)
genes <- setdiff(rownames(expr)[pct >= GENE_MIN_PCT], GENE_OF_INTEREST)
expr_t <- Matrix::t(expr[genes, , drop = FALSE])       # cells × genes (열 추출 O(1))
rm(expr); gc(verbose = FALSE)
message(sprintf("[유전자] %s개 (검출률 ≥ %.0f%%)", format(length(genes), big.mark = ","),
                GENE_MIN_PCT * 100))

## ── 3. 유전자별 LMM — refit() 으로 가속 ─────────────────────────────────────
# lmer 를 매번 처음부터 적합하면 매우 느리다. 모형 구조가 동일하므로
# 한 번만 적합한 뒤 refit(newresp=) 로 반응변수만 갈아끼우면 수십 배 빨라진다.
FORM <- y ~ g_within + g_between + d_within + d_between + (1 | donor) + (1 | dataset)
CTRL <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5),
                          calc.derivs = FALSE)

df0 <- meta; df0$y <- as.numeric(expr_t[, 1])
m0  <- lme4::lmer(FORM, data = df0, control = CTRL)
message("[적합] 기준 모형 완료 → refit() 으로 ", length(genes), "개 유전자 처리")

fit_one <- function(j) {
  y <- as.numeric(expr_t[, j])
  if (stats::sd(y) == 0) return(NULL)
  m <- try(suppressWarnings(suppressMessages(lme4::refit(m0, newresp = y))), silent = TRUE)
  if (inherits(m, "try-error")) return(NULL)
  cf <- summary(m)$coefficients
  if (!all(c("g_within", "g_between") %in% rownames(cf))) return(NULL)
  c(b_w = cf["g_within", "Estimate"],  se_w = cf["g_within", "Std. Error"],
    b_b = cf["g_between", "Estimate"], se_b = cf["g_between", "Std. Error"],
    b_dw = cf["d_within", "Estimate"], b_db = cf["d_between", "Estimate"])
}

t0 <- Sys.time()
res_list <- if (N_CORES > 1 && .Platform$OS.type != "windows") {
  parallel::mclapply(seq_along(genes), fit_one, mc.cores = N_CORES)
} else {
  pbapply::pblapply(seq_along(genes), fit_one)
}
message(sprintf("[적합] 소요 %.1f 분", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

ok  <- !vapply(res_list, is.null, logical(1))
mat <- do.call(rbind, res_list[ok])
stats_df <- data.frame(gene = genes[ok], mat, stringsAsFactors = FALSE) %>%
  mutate(
    t_within  = b_w / se_w,
    t_between = b_b / se_b,
    # within 은 세포 수(수만)라 정규근사, between 은 donor 수에 묶이므로 t 분포 사용
    p_within  = 2 * stats::pnorm(-abs(t_within)),
    p_between = 2 * stats::pt(-abs(t_between), df = max(2, N_DONOR - 3)),
    fdr_within  = p.adjust(p_within,  "BH"),
    fdr_between = p.adjust(p_between, "BH"))
write_result(stats_df, "Fig6_wb_gene_stats.csv")
cat(sprintf("\n[유전자] FDR<0.05 — within %d / between %d (총 %d)\n",
            sum(stats_df$fdr_within < 0.05), sum(stats_df$fdr_between < 0.05),
            nrow(stats_df)))

saveRDS(res_list, file.path(DIR_DATA_PROC, "UCell_gene_scorelist.rds"))
# res_list <- readRDS(file.path(DIR_DATA_PROC, "UCell_gene_scorelist.rds"))
## ── 4. GSEA: fgsea + cameraPR 를 항상 쌍으로 ────────────────────────────────
HALLMARK <- get_hallmark_list()
GOBP <- get_gobp_list()
GOMF <- get_gomf_list()
GOCC <- get_gocc_list()

run_gsea_pair <- function(stat_vec, level, term) {
  stat_vec <- sort(stat_vec[is.finite(stat_vec)], decreasing = TRUE)
  set.seed(SEED)
  fg <- fgsea::fgsea(pathways = term, stats = stat_vec,
                     minSize = 10, maxSize = 500, eps = 0, nPermSimple = 100000)
  idx <- limma::ids2indices(term, names(stat_vec))
  idx <- idx[vapply(idx, length, 1L) >= 10]
  cam <- limma::cameraPR(stat_vec, idx, inter.gene.cor = INTER_GENE_COR,
                         use.ranks = FALSE)
  cam$pathway <- rownames(cam)
  out <- fg %>%
    dplyr::select(pathway, NES, pval, padj, size) %>%
    dplyr::left_join(cam %>% dplyr::select(pathway, Direction,
                                           camera_p = PValue, camera_FDR = FDR),
                     by = "pathway") %>%
    dplyr::mutate(level = level) %>%
    dplyr::arrange(pval)
  out
}

gsea_w <- run_gsea_pair(setNames(stats_df$t_within,  stats_df$gene), "within", HALLMARK)
gsea_w_BP <- run_gsea_pair(setNames(stats_df$t_within,  stats_df$gene), "within", GOBP)
gsea_w_MF <- run_gsea_pair(setNames(stats_df$t_within,  stats_df$gene), "within", GOMF)
gsea_w_CC <- run_gsea_pair(setNames(stats_df$t_within,  stats_df$gene), "within", GOCC)

gsea_b <- run_gsea_pair(setNames(stats_df$t_between, stats_df$gene), "between", HALLMARK) #nPermSimple = 1000000로 해도 불안정, 버리고 pseudobulk를 쓰기로 결정.
write_result(gsea_w, "Fig6_wb_gsea_within.csv")
write_result(gsea_w_BP, "Fig6_wb_gsea_within_BP.csv")
write_result(gsea_w_MF, "Fig6_wb_gsea_within_MF.csv")
write_result(gsea_w_CC, "Fig6_wb_gsea_within_CC.csv")

write_result(gsea_b, "Fig6_wb_gsea_between.csv")

cat("\n===== WITHIN (세포 내 공발현) =====\n")
print(as.data.frame(head(gsea_w, 12)), digits = 3)
cat("\n===== BETWEEN (환자 층위) =====\n")
print(as.data.frame(head(gsea_b, 12)), digits = 3)
cat(sprintf("\n[요약] fgsea padj<0.05 — within %d / between %d\n",
            sum(gsea_w$padj < 0.05, na.rm = TRUE), sum(gsea_b$padj < 0.05, na.rm = TRUE)))
cat(sprintf("       camera FDR<0.05 — within %d / between %d  ← 이쪽을 근거로 사용\n",
            sum(gsea_w$camera_FDR < 0.05, na.rm = TRUE),
            sum(gsea_b$camera_FDR < 0.05, na.rm = TRUE)))

## ── 5. 시각화 ───────────────────────────────────────────────────────────────
both <- dplyr::bind_rows(gsea_w, gsea_b) %>%
  dplyr::mutate(sig_cam = camera_FDR < 0.05,
                lab = sub("HALLMARK_", "", pathway))

# (a) within vs between NES 산점도 — 두 층위가 같은 축을 가리키는가
wide <- both %>% dplyr::select(lab, level, NES, sig_cam) %>%
  tidyr::pivot_wider(names_from = level, values_from = c(NES, sig_cam))
pA <- ggplot(wide, aes(NES_within, NES_between)) +
  geom_hline(yintercept = 0, color = "grey70") + geom_vline(xintercept = 0, color = "grey70") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey60") +
  geom_point(aes(color = sig_cam_within | sig_cam_between), size = 2.2, alpha = .85) +
  ggrepel::geom_text_repel(
    data = subset(wide, sig_cam_within | sig_cam_between | abs(NES_within) > 2 |
                    abs(NES_between) > 2), aes(label = lab), size = 2.7, max.overlaps = 20) +
  scale_color_manual(values = c(`TRUE` = "#D64B4B", `FALSE` = "grey70"),
                     name = "camera FDR<0.05") +
  labs(x = "NES (within-donor, cell level)", y = "NES (between-donor, patient level)",
       title = "Within- vs between-level pathway enrichment") +
  theme_paper(base_size = 11)
ggsave(file.path(DIR_FIGURES, "Fig6_wb_gsea_scatter.pdf"), pA, width = 6.4, height = 5.6)

# (b) 층위별 상위 경로 막대 (camera 로 유의한 것만 강조)
top <- both %>% group_by(level) %>% slice_min(pval, n = 12) %>% ungroup()
pB <- ggplot(top, aes(NES, reorder(lab, NES), fill = sig_cam)) +
  geom_col() + facet_wrap(~level, scales = "free_y") +
  scale_fill_manual(values = c(`TRUE` = "#D64B4B", `FALSE` = "grey78"),
                    name = "camera FDR<0.05") +
  labs(x = "NES", y = NULL, title = "Top pathways by level",
       caption = "회색 = fgsea 로만 유의(상관 보정 후 소멸) → 근거로 사용하지 않음") +
  theme_paper(base_size = 10.5) +
  theme(plot.caption = element_text(hjust = 0, size = 7.5, color = "grey35"))
ggsave(file.path(DIR_FIGURES, "Fig6_wb_gsea_top.pdf"), pB, width = 9.5, height = 5.4)

message("완료: within/between 세포 수준 GSEA")

# -----------------------------------------------------------------------------
# 해석 · 보고 가이드
#  · between-GSEA 는 pseudobulk-GSEA 와 "같은 질문"(환자 층위)에 답한다.
#    다만 depth 를 양 층위에서 보정하고 dataset 랜덤효과를 넣으므로 더 통제적이다.
#    → 본문에는 between 을 싣고, pseudobulk 결과는 supplementary 일치성 확인으로 두면
#      "두 독립적 추정법이 같은 결론"이라는 추가 근거가 된다(하나만 고르지 않아도 됨).
#  · within-GSEA 는 새로운 질문(세포 내재적 공발현)이며 pseudobulk 로는 얻을 수 없다.
#    단 세포 수준이라 depth 교란에 가장 취약하므로 d_within 계수를 함께 점검할 것.
#  · fgsea 와 cameraPR 이 어긋나면 camera 를 따른다. INTER_GENE_COR 는 기본 0.01 이며,
#    Fig4_GSEA_camera_vs_fgsea.csv 에 관측 상관이 있으면 경로별 값으로 바꿔 넣으면 된다.
#
# 런타임 메모
#  · refit() 사용으로 유전자당 수십 ms 수준. 8,000 유전자 × 55,000 세포 기준
#    대략 10~40분(장비에 따라). 느리면 CELLS_PER_DONOR 를 200~300 으로 낮출 것.
#  · CELLS_PER_DONOR 는 within 검정력만 낮추고 between(donor 수) 은 영향받지 않는다.
# -----------------------------------------------------------------------------
