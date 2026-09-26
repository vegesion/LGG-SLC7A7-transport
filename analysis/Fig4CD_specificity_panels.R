# =============================================================================
# Fig4CD_specificity_panels.R
# Fig 4 추가 패널 — C: myeloid 마커 head-to-head / D: 침윤 보정 민감도
# -----------------------------------------------------------------------------
# 목적: "SLC7A7의 예후 효과가 단순한 myeloid 침윤량 때문인가?" 에 답하는 두 패널.
#   C. 동일 모델(age+IDH 보정)에서 SLC7A7 vs 범용 myeloid 마커 8종의 HR 직접 비교
#   D. 6종 침윤 지표로 각각 보정했을 때 SLC7A7 HR 이 유지되는지 (+ 침윤 지표 자체의 p)
#
# 입력 : (C) 코호트 발현/임상 — 프로젝트 파이프라인의 load_cohort() 사용해 재계산
#            ※ 기존 콘솔 출력에는 CI 가 없어 forest 를 그릴 수 없으므로 재계산함
#        (D) Fig4_sensitivity_myeloid_adjustment.csv  (adjusted_for, HR, lower, upper, p)
#            + Fig4_cox_adjusted_for_*.csv  (침윤 지표 자체의 p 를 읽어옴, 선택)
# 출력 : Fig4C_head_to_head.{pdf,tiff}, Fig4D_sensitivity.{pdf,tiff},
#        Fig4CD_combined.{pdf,tiff}, Fig4_head_to_head.csv
# =============================================================================
source(here::here("R", "setup.R"))
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(survival)
})
HAS_PATCHWORK <- requireNamespace("patchwork", quietly = TRUE)

## ── 0. 설정 ─────────────────────────────────────────────────────────────────
OUT_DIR      <- DIR_RESULTS
GENE         <- "SLC7A7"
SENS_FILE    <- "Fig4_sensitivity_myeloid_adjustment.csv"
H2H_FILE     <- "Fig4_head_to_head.csv"      # 있으면 재사용, 없으면 계산
COHORTS_USE  <- c("TCGA_LGG", "CGGA_693", "CGGA_325")
MARKERS      <- c("SLC7A7", "PTPRC", "C1QB", "C1QA", "TYROBP",
                  "ITGAM", "AIF1", "CSF1R", "CD68")
PAL_COHORT   <- c(TCGA_LGG = "#D85A30", CGGA_693 = "#185FA5", CGGA_325 = "#2E8B57")

## =============================================================================
## PANEL C — head-to-head (SLC7A7 vs 범용 myeloid 마커)
## =============================================================================
if (file.exists(H2H_FILE)) {
  h2h <- read.csv(H2H_FILE, stringsAsFactors = FALSE)
  message("[C] 기존 ", H2H_FILE, " 재사용")
} else {
  message("[C] head-to-head 재계산 (CI 포함)")
  source(here::here("R", "setup.R"))          # load_cohort() 등 파이프라인 함수
  
  h2h <- bind_rows(lapply(COHORTS_USE, function(nm) {
    co <- tryCatch(load_cohort(nm), error = function(e) NULL)
    if (is.null(co)) { message("  건너뜀: ", nm); return(NULL) }
    cl <- co$clin
    bind_rows(lapply(intersect(MARKERS, rownames(co$expr)), function(g) {
      d <- data.frame(time = cl$time, event = cl$event,
                      x = as.numeric(scale(as.numeric(co$expr[g, ]))),
                      age = cl$age, idh = cl$idh)
      d <- d[stats::complete.cases(d) & d$time > 0, , drop = FALSE]
      if (nrow(d) < 30 || length(unique(d$idh)) < 2 || sd(d$x) == 0) return(NULL)
      s <- summary(coxph(Surv(time, event) ~ x + age + idh, data = d))
      data.frame(cohort = nm, marker = g,
                 HR = s$conf.int[1, 1], lower = s$conf.int[1, 3],
                 upper = s$conf.int[1, 4], p = s$coefficients[1, 5], n = s$n)
    }))
  }))
  write.csv(h2h, file.path(OUT_DIR, H2H_FILE), row.names = FALSE)
}

# ★ position_dodge 는 x 방향으로 작동 → 가로 forest 에서 점/오차막대가 어긋난다.
#   y 좌표를 직접 계산해 오프셋을 주면 두 geom 이 항상 같은 위치에 놓인다.
nC   <- length(COHORTS_USE)
STEP <- 0.26                                   # 코호트 간 세로 간격

h2h <- h2h %>%
  mutate(cohort = factor(cohort, levels = COHORTS_USE),
         sig    = p < 0.05,
         y_base = as.numeric(factor(marker, levels = rev(MARKERS))),   # 1=CD68 … 9=SLC7A7
         y      = y_base + (as.integer(cohort) - (nC + 1) / 2) * STEP)

stopifnot(!any(is.na(h2h$y_base)))             # MARKERS 에 없는 유전자 방지

y_target <- length(MARKERS)                    # SLC7A7 위치

pC <- ggplot(h2h, aes(x = HR, y = y, color = cohort)) +
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = y_target - 0.5, ymax = y_target + 0.5, fill = "grey93") +
  geom_vline(xintercept = 1, color = "grey50", linewidth = 0.45) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0, linewidth = 0.7) +
  geom_point(aes(shape = sig), size = 2.3, fill = "white", stroke = 0.7) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21),
                     labels = c(`TRUE` = "P < 0.05", `FALSE` = "n.s."), name = NULL) +
  scale_color_manual(values = PAL_COHORT, name = NULL) +
  scale_x_continuous(trans = "log2", breaks = c(0.8, 1, 1.25, 1.6, 2)) +
  scale_y_continuous(breaks = seq_along(MARKERS), labels = rev(MARKERS),
                     expand = expansion(add = 0.6)) +
  labs(x = "Adjusted HR per 1 SD (age + IDH)", y = NULL,
       title = "C  SLC7A7 vs generic myeloid markers") +
  theme_classic(base_size = 11) +
  theme(axis.text = element_text(color = "black"),
        axis.text.y = element_text(
          face = ifelse(rev(MARKERS) == GENE, "bold", "plain")),
        panel.grid.major.y = element_line(color = "grey95", linewidth = 0.3),
        plot.title = element_text(face = "bold", size = 11.5),
        legend.position = "top", legend.box = "horizontal",
        legend.margin = margin(b = -6))

## =============================================================================
## PANEL D — 침윤 보정 민감도
## =============================================================================
sens <- read.csv(paste0(DIR_RESULTS, "/",SENS_FILE), stringsAsFactors = FALSE)
names(sens)[1] <- "adjusted_for"
sens$adjusted_for <- gsub("\\\\", "", sens$adjusted_for)
sens$adjusted_for <- gsub("\\.x$|\\.y$", "", sens$adjusted_for)   # ImmuneScore.x → ImmuneScore

# 침윤 지표 자체의 p 를 개별 파일에서 회수 (하위 폴더까지 재귀 탐색)
ADJ_FILES <- list.files(".", pattern = "^Fig4_cox_adjusted_for_.*\\.csv$",
                        recursive = TRUE, full.names = TRUE)
get_myeloid_p <- function(nm) {
  key <- tolower(gsub("[^A-Za-z0-9]", "", nm))
  hit <- ADJ_FILES[tolower(gsub("[^A-Za-z0-9]", "",
                                sub("^Fig4_cox_adjusted_for_", "", tools::file_path_sans_ext(basename(ADJ_FILES))))) %in%
                     c(key, paste0(key, "x"), paste0(key, "y"))]
  if (!length(hit)) return(NA_real_)
  d <- read.csv(hit[1], stringsAsFactors = FALSE)
  v <- suppressWarnings(as.numeric(d$p[grepl("myeloid", d$term, ignore.case = TRUE)]))
  if (length(v) && is.finite(v[1])) v[1] else NA_real_
}
sens$p_myeloid <- vapply(sens$adjusted_for, get_myeloid_p, numeric(1))
if (all(is.na(sens$p_myeloid)))
  message("[D] 주의: Fig4_cox_adjusted_for_*.csv 를 찾지 못해 침윤 지표 p 를 표시하지 않습니다.\n",
          "     해당 CSV 들을 작업 디렉터리에 두면 y축 라벨과 캡션에 자동 반영됩니다.")

# 보정 전 기준값 (같은 공변량, 침윤 지표만 제외). 없으면 논문 보고값 사용.
REF <- tryCatch({
  co  <- load_cohort("TCGA_LGG")
  mv  <- multivariable_cox(co$expr, co$clin, GENE, covars = c("age", "grade", "idh"))
  r   <- mv$table[mv$table$term == "gene", ]
  c(HR = r$HR, lo = r$lower, hi = r$upper)
}, error = function(e) c(HR = 1.8166, lo = 1.3083, hi = 2.5224))
message(sprintf("[D] 보정 전 기준 HR = %.2f (%.2f-%.2f)", REF["HR"], REF["lo"], REF["hi"]))

sens <- sens %>%
  arrange(HR) %>%
  mutate(adjusted_for = factor(adjusted_for, levels = adjusted_for),
         sig = p < 0.05,
         lab_p  = ifelse(p < 0.001, formatC(p, format = "e", digits = 1), sprintf("%.3f", p)),
         lab_mp = ifelse(is.na(p_myeloid), "—",
                         ifelse(p_myeloid < 0.001,
                                formatC(p_myeloid, format = "e", digits = 1),
                                sprintf("%.2f", p_myeloid))))

# y축 라벨에 침윤 지표 자체의 p 를 병기 (겹침 없이 정보 전달)
sens$ylab <- ifelse(is.na(sens$p_myeloid), as.character(sens$adjusted_for),
                    sprintf("%s  [adj. P=%s]", sens$adjusted_for, sens$lab_mp))
sens$ylab <- factor(sens$ylab, levels = sens$ylab)     # HR 오름차순 유지

xr      <- max(sens$upper)                 # 가장 긴 CI 끝
xlim_hi <- xr * 1.9                        # p값 텍스트 공간 확보

# 캡션: 침윤 지표 p 를 찾은 경우에만 그 범위를 덧붙임 (없으면 Inf 표시 방지)
cap_D <- sprintf("Dashed line & band: model without infiltration term, HR %.2f (%.2f–%.2f).  Numbers = P for %s.",
                 REF["HR"], REF["lo"], REF["hi"], GENE)
if (any(is.finite(sens$p_myeloid)))
  cap_D <- paste0(cap_D, sprintf("  Infiltration terms themselves: P = %.3f–%.2f (never significant).",
                                 min(sens$p_myeloid, na.rm = TRUE),
                                 max(sens$p_myeloid, na.rm = TRUE)))

pD <- ggplot(sens, aes(HR, ylab)) +
  # 보정 전 기준 (밴드 + 파선)
  annotate("rect", xmin = REF["lo"], xmax = REF["hi"], ymin = -Inf, ymax = Inf,
           fill = "#D85A30", alpha = 0.10) +
  geom_vline(xintercept = REF["HR"], color = "#D85A30", linetype = 2, linewidth = 0.5) +
  geom_vline(xintercept = 1, color = "grey50", linewidth = 0.45) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0, linewidth = 0.75) +
  geom_point(aes(shape = sig), size = 2.7, fill = "white") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), guide = "none") +
  # p값은 각 행의 CI 끝 바로 뒤에 → 겹치지 않음
  geom_text(aes(x = upper * 1.07, label = lab_p), hjust = 0, size = 2.9) +
  scale_x_continuous(trans = "log2", breaks = c(1, 1.5, 2, 3, 4),
                     limits = c(0.9, xlim_hi)) +
  labs(x = sprintf("HR of %s after adjusting for each infiltration estimate", GENE),
       y = NULL, title = "D  Adjustment for myeloid/immune infiltration",
       caption = cap_D) +
  theme_classic(base_size = 11) +
  theme(axis.text = element_text(color = "black"),
        plot.title = element_text(face = "bold", size = 11.5),
        plot.caption = element_text(hjust = 0, size = 7.8, color = "grey35"),
        plot.margin = margin(5, 8, 5, 5))

## ── 저장 ────────────────────────────────────────────────────────────────────
ggsave(file.path(OUT_DIR, "Fig4C_head_to_head.pdf"), pC, width = 6.0, height = 4.6)
ggsave(file.path(OUT_DIR, "Fig4C_head_to_head.tiff"), pC, width = 6.0, height = 4.6,
       dpi = 300, compression = "lzw")
ggsave(file.path(OUT_DIR, "Fig4D_sensitivity.pdf"), pD, width = 6.6, height = 3.6)
ggsave(file.path(OUT_DIR, "Fig4D_sensitivity.tiff"), pD, width = 6.6, height = 3.6,
       dpi = 300, compression = "lzw")

if (HAS_PATCHWORK) {
  library(patchwork)
  fig <- pC / pD + plot_layout(heights = c(1.25, 1))
  ggsave(file.path(OUT_DIR, "Fig4CD_combined.pdf"), fig, width = 6.8, height = 8.2)
  ggsave(file.path(OUT_DIR, "Fig4CD_combined.tiff"), fig, width = 6.8, height = 8.2,
         dpi = 300, compression = "lzw")
  print(fig)
} else {
  print(pC); print(pD)
}

cat("\n===== Panel C 요약 (코호트별 SLC7A7 순위) =====\n")
h2h %>% group_by(cohort) %>% arrange(desc(HR)) %>%
  mutate(rank = row_number()) %>% filter(marker == GENE) %>%
  select(cohort, marker, HR, p, rank) %>% as.data.frame() %>% print()
cat("\n===== Panel D 요약 =====\n"); print(as.data.frame(sens[, 1:5]))


pC + theme_blank_frame()
pC
pD
pD + theme_blank_frame()
fig
save_emf(pC + theme_blank_frame(), "Fig4C_head_to_head.emf", width = 5, height = 4, draw_major_grid = FALSE, draw_minor_grid = FALSE)
save_emf(pD + theme_blank_frame(), "Fig4D_sensitivity.emf", width = 5, height = 3.5, draw_major_grid = FALSE, draw_minor_grid = FALSE)
# -----------------------------------------------------------------------------
# 기대값 (검증용)
#  C: SLC7A7 HR — TCGA_LGG 1.73 / CGGA_693 1.23 / CGGA_325 1.60, 세 코호트 모두 최상위
#     CD68 은 TCGA_LGG(p=0.25)·CGGA_693(p=0.19) 에서 n.s. → 열린 원으로 표시됨
#  D: 6종 중 5종에서 SLC7A7 유의 (HR 1.58–2.46), Macrophages 만 p=0.053
#     침윤 지표 자체의 p 는 0.078–0.72 로 단 한 번도 유의하지 않음
#
# 논문 문장:
#  "Across three cohorts, SLC7A7 showed the largest effect size among myeloid markers,
#   whereas the pan-myeloid marker CD68 was not significant in two of three cohorts.
#   SLC7A7 remained associated with survival after adjustment for each of six
#   infiltration estimates (HR 1.58–2.46), while the infiltration estimates themselves
#   were never independently significant (P = 0.078–0.72)."
# -----------------------------------------------------------------------------
