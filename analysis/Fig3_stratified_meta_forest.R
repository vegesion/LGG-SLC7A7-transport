# =============================================================================
# Fig3_stratified_meta_forest.R
# Fig 3 추가 패널 — IDH 층화 forest + 층별 메타분석 (pooled diamond, I²)
# -----------------------------------------------------------------------------
# 입력 : Fig3_stratified_cox_idh.csv  (term, HR, lower, upper, p, n, stratum, cohort)
# 출력 : Fig3D_stratified_meta_forest.pdf / .tiff , Fig3_meta_summary.csv
#
# 통계 : 층화 Cox 의 HR·95%CI 에서 SE 역산 → 역분산 가중 메타분석
#   SE      = (log(upper) - log(lower)) / (2 × 1.96)
#   고정효과 β = Σ(wβ)/Σw ,  w = 1/SE² ,  SE_pooled = √(1/Σw)
#   Cochran's Q = Σw(β - β_pooled)² ,  I² = max(0, (Q-df)/Q) × 100
#   DerSimonian–Laird τ² 로 랜덤효과 병기 (I² > 50 이면 랜덤효과를 대표값으로)
#
# 검증된 기대값:
#   IDH-mutant   k=3  HR 1.465 (1.279–1.678)  p = 3.4e-08  I² = 0%
#   IDH-wildtype k=4  고정 1.175 / 랜덤 1.251 (0.998–1.568)  I² = 77.7%
# =============================================================================
source(here::here("R", "setup.R"))
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2)
})
HAS_PATCHWORK <- requireNamespace("patchwork", quietly = TRUE)

## ── 0. 설정 ─────────────────────────────────────────────────────────────────
IN_FILE      <- file.path(DIR_RESULTS, "Fig3_stratified_cox_idh.csv")   # ← 구글드라이브에서 내려받은 경로
OUT_DIR      <- DIR_RESULTS
GENE         <- "SLC7A7"
COHORT_ORDER <- c("TCGA_LGG", "TCGA_GBM", "CGGA_693", "CGGA_325")
STRATA_ORDER <- c("IDH-mutant", "IDH-wildtype")
PAL          <- c("IDH-mutant" = "black", "IDH-wildtype" = "black")

## ── 1. 로드 + SE 역산 ───────────────────────────────────────────────────────
raw <- read.csv(IN_FILE, stringsAsFactors = FALSE)
raw$cohort  <- gsub("\\\\", "", raw$cohort)      # escape 문자 제거
raw$stratum <- gsub("\\\\", "", raw$stratum)

dat <- raw %>%
  filter(term == "gene") %>%
  mutate(
    stratum = factor(ifelse(grepl("Mutant", stratum, ignore.case = TRUE),
                            "IDH-mutant", "IDH-wildtype"), levels = STRATA_ORDER),
    cohort  = factor(cohort, levels = COHORT_ORDER),
    logHR   = log(HR),
    se      = (log(upper) - log(lower)) / (2 * qnorm(0.975)),
    w       = 1 / se^2
  ) %>%
  filter(!is.na(cohort)) %>%
  arrange(stratum, cohort)

stopifnot(nrow(dat) > 0, all(is.finite(dat$se)))

## ── 2. 층별 메타분석 ────────────────────────────────────────────────────────
meta_one <- function(d) {
  w <- d$w; b <- d$logHR; k <- nrow(d); df <- k - 1
  bf <- sum(w * b) / sum(w); sef <- sqrt(1 / sum(w))
  Q  <- sum(w * (b - bf)^2)
  I2 <- if (Q > 0 && df > 0) max(0, (Q - df) / Q) * 100 else 0
  tau2 <- if (df > 0) max(0, (Q - df) / (sum(w) - sum(w^2) / sum(w))) else 0
  wr <- 1 / (d$se^2 + tau2)
  br <- sum(wr * b) / sum(wr); ser <- sqrt(1 / sum(wr))
  data.frame(
    k = k, n_total = sum(d$n), Q = Q, I2 = I2, tau2 = tau2,
    HR_fixed = exp(bf), lo_fixed = exp(bf - 1.96*sef), hi_fixed = exp(bf + 1.96*sef),
    p_fixed  = 2 * pnorm(-abs(bf / sef)),
    HR_random = exp(br), lo_random = exp(br - 1.96*ser), hi_random = exp(br + 1.96*ser),
    p_random  = 2 * pnorm(-abs(br / ser))
  )
}

meta <- dat %>% group_by(stratum) %>% group_modify(~ meta_one(.x)) %>% ungroup() %>%
  mutate(model = ifelse(I2 > 50, "random", "fixed"),
         HR = ifelse(model == "random", HR_random, HR_fixed),
         lo = ifelse(model == "random", lo_random, lo_fixed),
         hi = ifelse(model == "random", hi_random, hi_fixed),
         p  = ifelse(model == "random", p_random,  p_fixed))

cat("\n===== 메타분석 요약 =====\n"); print(as.data.frame(meta))
write.csv(meta, file.path(OUT_DIR, "Fig3_meta_summary.csv"), row.names = FALSE)

## ── 3. 플롯 행 구성 (층별: 코호트들 → pooled) ───────────────────────────────
blocks <- lapply(levels(dat$stratum), function(st) {
  d <- dat[dat$stratum == st, , drop = FALSE]
  m <- meta[meta$stratum == st, , drop = FALSE]
  rbind(
    data.frame(stratum = st, label = as.character(d$cohort),
               HR = d$HR, lo = d$lower, hi = d$upper, p = d$p, n = d$n,
               weight = d$w / sum(d$w), type = "cohort", stringsAsFactors = FALSE),
    data.frame(stratum = st, label = sprintf("Pooled (%s)", m$model),
               HR = m$HR, lo = m$lo, hi = m$hi, p = m$p, n = m$n_total,
               weight = NA_real_, type = "pooled", stringsAsFactors = FALSE)
  )
})
rows <- do.call(rbind, blocks)
rows$stratum <- factor(rows$stratum, levels = STRATA_ORDER)

# y 좌표: 위(첫 층)에서 아래로, 층 사이 1칸 간격
idx <- seq_len(nrow(rows))
gap <- as.integer(rows$stratum) - 1L
pos <- idx + gap
rows$y <- max(pos) - pos + 1

## ── 4. Pooled 다이아몬드 ────────────────────────────────────────────────────
pooled <- rows[rows$type == "pooled", , drop = FALSE]
dia <- do.call(rbind, lapply(seq_len(nrow(pooled)), function(i) {
  r <- pooled[i, ]
  data.frame(id = paste0("d", i), stratum = r$stratum,
             x = c(r$lo, r$HR, r$hi, r$HR),
             y = c(r$y, r$y + 0.30, r$y, r$y - 0.30))
}))

## ── 5. Panel A — RevMan 5 스타일 forest ─────────────────────────────────────
## RevMan 팔레트 (취향껏 조정)
SQ_FILL  <- "black"   # study 사각형 (RevMan 녹색)
SQ_BORD  <- "black"     # 사각형 테두리
CI_COL   <- "black"     # CI 선
DIA_FILL <- "black"     # pooled 다이아몬드
BAND_G   <- "grey93"    # zebra 회색 행
BAND_W   <- "white"     # zebra 흰색 행

## 5-1. 서브그룹 헤더 행 (§3에서 만들어둔 층 사이 빈 슬롯을 그대로 사용)
hdr <- rows %>%
  group_by(stratum) %>%
  summarise(y = max(y) + 1, .groups = "drop") %>%
  mutate(label = as.character(stratum), type = "header")

top   <- max(hdr$y)
hdr_y <- top + 0.95                      # 컬럼 헤더 행 (pB와 공유)
arr_y <- min(rows$y) - 1.25              # 하단 화살표
ylim  <- c(arr_y - 0.95, hdr_y + 0.75)

## 5-2. zebra 배경 밴드 (위에서부터 번갈아)
band <- data.frame(y = seq_len(top))
band$fill <- ifelse((top - band$y) %% 2 == 1, BAND_G, BAND_W)

## 5-3. y축 라벨 = 헤더 + 들여쓴 study/subtotal (반드시 y 오름차순)
ax <- rbind(
  data.frame(y = rows$y, label = paste0("    ", rows$label), type = rows$type),
  hdr[, c("y", "label", "type")],
  data.frame(y = hdr_y, label = "Cohort", type = "colhead")
)
ax$label[ax$type == "pooled"] <- "    Subtotal (95% CI)"
ax <- ax[order(ax$y), ]
ax_face <- ifelse(ax$type %in% c("header", "colhead", "pooled"), "bold", "plain")

## 5-4. x 범위/눈금
xmin <- min(rows$lo, na.rm = TRUE) * 0.85
xmax <- max(rows$hi, na.rm = TRUE) * 1.15
brk  <- c(0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4)
brk  <- brk[brk >= xmin & brk <= xmax]

pA <- ggplot() +
  # 배경 밴드 (제일 먼저)
  geom_rect(data = band, aes(ymin = y - 0.5, ymax = y + 0.5),
            xmin = -Inf, xmax = Inf, fill = band$fill) +
  # 컬럼 헤더 아래 구분선
  annotate("segment", x = xmin, xend = xmax, y = top + 0.5, yend = top + 0.5,
           color = "grey40", linewidth = 0.4) +
  # 무효과선
  geom_vline(xintercept = 1, color = "black", linewidth = 0.4) +
  # study CI + 가중치 비례 사각형
  geom_errorbarh(data = subset(rows, type == "cohort"),
                 aes(y = y, xmin = lo, xmax = hi),
                 height = 0, linewidth = 0.55, color = CI_COL) +
  geom_point(data = subset(rows, type == "cohort"),
             aes(x = HR, y = y, size = weight),
             shape = 22, fill = SQ_FILL, color = SQ_BORD, stroke = 0.35) +
  # pooled 다이아몬드
  geom_polygon(data = dia, aes(x = x, y = y, group = id),
               fill = DIA_FILL, color = "black", linewidth = 0.3) +
  # 하단 방향 화살표
  annotate("segment", x = 1, xend = xmin * 1.04, y = arr_y, yend = arr_y,
           arrow = grid::arrow(length = grid::unit(0.09, "in"), type = "closed"),
           linewidth = 0.4) +
  annotate("segment", x = 1, xend = xmax * 0.96, y = arr_y, yend = arr_y,
           arrow = grid::arrow(length = grid::unit(0.09, "in"), type = "closed"),
           linewidth = 0.4) +
  scale_size_area(max_size = 6, guide = "none") +
  scale_x_continuous(trans = "log2", breaks = brk, labels = brk,
                     limits = c(xmin, xmax), expand = c(0, 0)) +
  scale_y_continuous(breaks = ax$y, labels = ax$label,
                     limits = ylim, expand = c(0, 0)) +
  labs(x = sprintf("Adjusted HR per 1 SD %s (95%% CI)", GENE), y = NULL) +
  theme_classic(base_size = 11) +
  theme(axis.text        = element_text(color = "black"),
        axis.text.y      = element_text(hjust = 0, face = ax_face, size = 9.2),
        axis.line.y      = element_blank(),
        axis.ticks.y     = element_blank(),
        axis.line.x      = element_line(color = "black", linewidth = 0.4),
        axis.title.x     = element_text(size = 9.5, margin = margin(t = 4)),
        legend.position  = "none",
        plot.margin      = margin(5, 2, 5, 5))

## ── 6. Panel B — 우측 수치 표 (RevMan 컬럼: Weight 추가) ────────────────────
fmt_p <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 1), sprintf("%.3f", p))
tab <- rows %>% mutate(hr_txt = sprintf("%.2f (%.2f–%.2f)", HR, lo, hi),
                       p_txt  = fmt_p(p),
                       n_txt  = as.character(n),
                       w_txt  = ifelse(type == "pooled", "100.0%",
                                       sprintf("%.1f%%", weight * 100)),
                       fw     = ifelse(type == "pooled", "bold", "plain"))
COL_X <- c(hr = 0.00, w = 1.08, p = 1.52, n = 2.06)

pB <- ggplot(tab) +
  geom_rect(data = band, aes(ymin = y - 0.5, ymax = y + 0.5),
            xmin = -Inf, xmax = Inf, fill = band$fill, inherit.aes = FALSE) +
  annotate("segment", x = -0.03, xend = 2.45, y = top + 0.5, yend = top + 0.5,
           color = "grey40", linewidth = 0.4) +
  geom_text(aes(x = COL_X["hr"], y = y, label = hr_txt, fontface = fw), hjust = 0, size = 3.05) +
  geom_text(aes(x = COL_X["w"],  y = y, label = w_txt,  fontface = fw), hjust = 0, size = 3.05,
            color = "grey30") +
  geom_text(aes(x = COL_X["p"],  y = y, label = p_txt,  fontface = fw), hjust = 0, size = 3.05) +
  geom_text(aes(x = COL_X["n"],  y = y, label = n_txt), hjust = 0, size = 3.05, color = "grey35") +
  annotate("text", x = as.numeric(COL_X), y = hdr_y,
           label = c("HR (95% CI)", "Weight", "P", "n"),
           hjust = 0, fontface = "bold", size = 3.1) +
  scale_x_continuous(limits = c(-0.03, 2.45), expand = c(0, 0)) +
  scale_y_continuous(limits = ylim, expand = c(0, 0)) +
  theme_void() + theme(plot.margin = margin(5, 5, 5, 0))

## ── 7. 결합 + 저장 ──────────────────────────────────────────────────────────
cap <- paste(sprintf("%s: I² = %.0f%%, k = %d, n = %d",
                     meta$stratum, meta$I2, meta$k, meta$n_total), collapse = "    |    ")

if (HAS_PATCHWORK) {
  library(patchwork)
  fig <- (pA | pB) + patchwork::plot_layout(widths = c(1, 0.85)) +          # Weight 컬럼 추가분 
    patchwork::plot_annotation(
      caption = cap,
      theme = theme(plot.caption = element_text(hjust = 0, size = 8.5, color = "grey30")))
} else {
  fig <- pA + labs(caption = cap)
}

H <- 0.40 * (nrow(rows) + nrow(hdr) + 3) + 1.6
ggsave(file.path(OUT_DIR, "Fig3D_stratified_meta_forest.pdf"), fig,
       width = 8.2, height = H, useDingbats = FALSE)
ggsave(file.path(OUT_DIR, "Fig3D_stratified_meta_forest.tiff"), fig,
       width = 8.2, height = H, dpi = 300, compression = "lzw")

cat("\n저장: Fig3D_stratified_meta_forest.pdf / .tiff\n")
cat("요약: Fig3_meta_summary.csv\n\n")
print(fig)


library(devEMF)
emf(paste0(file.path(DIR_FIGURES), "/Fig3D_stratified_meta_forest1.emf"), width = 12, height = 5)
print(fig)
dev.off()


pB
pA + theme_blank_frame()
raw_forest <- fig + theme_void()
raw_forest
save
save_emf(pA + theme_blank_frame(), "Fig3D_stratified_meta_forest_1.emf", width = 5, height = 3)
fig + theme_blank_frame()
# -----------------------------------------------------------------------------
# 논문 문장 초안
#  "In IDH-mutant glioma, higher SLC7A7 expression was associated with shorter overall
#   survival in all three cohorts with no detectable between-study heterogeneity
#   (pooled HR 1.47, 95% CI 1.28–1.68, P = 3.4×10⁻⁸, I² = 0%). In contrast, the
#   IDH-wildtype stratum was substantially heterogeneous (I² = 78%), driven by a null
#   association in CGGA-693 (HR 0.96, P = 0.62)."
# -----------------------------------------------------------------------------
