# =============================================================================
# Fig6_UCell_module_panels.R
# Fig 6 — donor 수준 UCell 모듈 회귀: forest + 음성대조 + 귀무분포 + 메타 + null
# -----------------------------------------------------------------------------
# 입력 (모두 UCell 분석 산출물):
#   UC_module_regression.csv  module, predictor(SLC7A7|depth_negctrl), maxRank,
#                             beta_lmer, se_lmer, p_lmer, fdr_lmer,
#                             beta_meta, se_meta, p_meta, fdr_meta, I2,
#                             k_dataset, n_donor, n_pos_dataset, n_gene, source
#   UC_null_distribution.csv  gene, module, beta, beta_meta   (무작위 유전자 200개)
#   UC_empirical_p.csv        module, null_mean, null_sd, beta_lmer, z_vs_null,
#                             p_emp_two, p_emp_BH
# 출력: Fig6A~D 개별 + Fig6_combined.{pdf,tiff} + Fig6_module_verdict.csv
#
# 패널 구성
#   A  모듈 forest + depth 음성대조   → 무엇이 유의하고 무엇이 depth 인공산물인가
#   B  경험적 귀무분포(무작위 200유전자) → 관측 효과가 우연 범위를 벗어나는가
#   C  within-dataset 메타 forest (I²)  → 배치를 넘어 재현되는가
#   D  null 패널(염증 축)               → "유의하지 않음"이 아니라 "효과가 없음"
# =============================================================================
source(here::here("R", "setup.R"))
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })
HAS_RIDGES    <- requireNamespace("ggridges",  quietly = TRUE)
HAS_PATCHWORK <- requireNamespace("patchwork", quietly = TRUE)

## ── 0. 설정 ─────────────────────────────────────────────────────────────────
OUT_DIR <- DIR_RESULTS
GENE    <- "SLC7A7"
MAXRANK <- 1500                      # 민감도용 3000 은 supplementary
FDR_CUT <- 0.05

CATEGORY <- c(
  GLYCOLYSIS = "Metabolic", PPP_REACTOME = "Metabolic", POLYAMINE = "Metabolic",
  ROS = "Metabolic", OXPHOS = "Metabolic", HYPOXIA = "Metabolic",
  ARG_ENZYME = "Metabolic", TRANSPORT = "Metabolic",
  IFN_ALPHA = "Inflammatory", IFN_GAMMA = "Inflammatory", IFN_REACTOME = "Inflammatory",
  NFKB = "Inflammatory", INFLAMMATORY = "Inflammatory",
  PAN_MYELOID = "Myeloid identity", HOMEOSTATIC = "Myeloid identity")
PAL_CAT  <- c(Metabolic = "#C0392B", Inflammatory = "#2471A3", `Myeloid identity` = "#7D6608")
PAL_PRED <- c(SLC7A7 = "#D64B4B", `depth (neg. control)` = "grey55")

## ── 1. 로드 및 판정(verdict) 계산 ───────────────────────────────────────────
reg <- read.csv(file.path(DIR_RESULTS, "/ucell/UC_module_regression.csv"), stringsAsFactors = FALSE)
emp <- read.csv(file.path(DIR_RESULTS, "/ucell/UC_empirical_p.csv"),       stringsAsFactors = FALSE)
nul <- read.csv(file.path(DIR_RESULTS, "/ucell/UC_null_distribution.csv"), stringsAsFactors = FALSE)

base <- reg %>% filter(maxRank == MAXRANK, !grepl("_disj$", module))
tgt  <- base %>% filter(predictor == GENE)
ctl  <- base %>% filter(predictor != GENE) %>%
  dplyr::select(module, beta_ctl = beta_lmer, se_ctl = se_lmer, fdr_ctl = fdr_lmer)

dat <- tgt %>%
  left_join(ctl, by = "module") %>%
  left_join(emp %>% dplyr::select(module, z_vs_null, p_emp_two, p_emp_BH), by = "module") %>%
  mutate(
    category   = factor(CATEGORY[module], levels = names(PAL_CAT)),
    sig_lmer   = fdr_lmer < FDR_CUT,
    sig_meta   = fdr_meta < FDR_CUT,
    ctl_clean  = !(fdr_ctl < FDR_CUT),                 # 음성대조에서 유의하지 않아야 함
    verdict = case_when(
      sig_lmer & !ctl_clean               ~ "depth-confounded",
      sig_lmer &  ctl_clean &  sig_meta   ~ "robust",
      sig_lmer &  ctl_clean & !sig_meta   ~ "specific (meta n.s.)",
      TRUE                                ~ "null"),
    lo = beta_lmer - 1.96 * se_lmer, hi = beta_lmer + 1.96 * se_lmer,
    lo_ctl = beta_ctl - 1.96 * se_ctl, hi_ctl = beta_ctl + 1.96 * se_ctl,
    lo_meta = beta_meta - 1.96 * se_meta, hi_meta = beta_meta + 1.96 * se_meta)

# 음성대조를 통과한 '유의' 모듈의 최대 효과크기 → Panel D 참조 밴드
REF_EFFECT <- max(abs(dat$beta_lmer[dat$sig_lmer & dat$ctl_clean]), na.rm = TRUE)
message(sprintf("[Fig6] n_donor = %d, k_dataset = %d, 참조 효과크기 = %.5f",
                dat$n_donor[1], dat$k_dataset[1], REF_EFFECT))
write.csv(dat %>% dplyr::select(module, category, beta_lmer, se_lmer, fdr_lmer,
                         beta_meta, fdr_meta, I2, n_pos_dataset, k_dataset,
                         beta_ctl, fdr_ctl, z_vs_null, p_emp_two, verdict),
          file.path(OUT_DIR, "Fig6_module_verdict.csv"), row.names = FALSE)
print(dat %>% count(verdict))

ORD <- dat %>% arrange(beta_lmer) %>% pull(module)     # 아래→위 오름차순

## =============================================================================
## Panel A — 모듈 forest + depth 음성대조
## =============================================================================
# ※ position_dodge 는 x 방향이라 가로 forest 에서 어긋남 → y 좌표를 직접 계산
STEP <- 0.20
pA_df <- bind_rows(
  dat %>% transmute(module, predictor = GENE, beta = beta_lmer, lo, hi,
                    sig = sig_lmer, off = +STEP),
  dat %>% transmute(module, predictor = "depth (neg. control)", beta = beta_ctl,
                    lo = lo_ctl, hi = hi_ctl, sig = fdr_ctl < FDR_CUT, off = -STEP)
) %>%
  mutate(y_base = match(module, ORD), y = y_base + off,
         predictor = factor(predictor, levels = names(PAL_PRED)))

pA <- ggplot(pA_df, aes(beta, y, color = predictor)) +
  geom_vline(xintercept = 0, color = "grey55", linewidth = 0.45) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.65) +
  geom_point(aes(shape = sig), size = 2.1, fill = "white", stroke = 0.7) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21),
                     labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "n.s."), name = NULL) +
  scale_color_manual(values = PAL_PRED, name = NULL) +
  scale_y_continuous(breaks = seq_along(ORD), labels = ORD,
                     expand = expansion(add = 0.7)) +
  labs(x = expression(beta~"(donor-level module score per unit "*italic("SLC7A7")*")"),
       y = NULL, title = "A  Module association with depth-matched negative control") +
  theme_classic(base_size = 10.5) +
  theme(axis.text = element_text(color = "black"),
        axis.text.y = element_text(
          color = PAL_CAT[as.character(dat$category[match(ORD, dat$module)])]),
        panel.grid.major.y = element_line(color = "grey96", linewidth = 0.3),
        plot.title = element_text(face = "bold", size = 11),
        legend.position = "top", legend.margin = margin(b = -6))

## =============================================================================
## Panel B — 경험적 귀무분포 (무작위 유전자 200개)
## =============================================================================
nul_p  <- nul %>% filter(module %in% dat$module) %>%
  mutate(module = factor(module, levels = ORD))
obs_p  <- dat %>% mutate(module = factor(module, levels = ORD),
                         lab = ifelse(p_emp_two < 0.001,
                                      formatC(p_emp_two, format = "e", digits = 1),
                                      sprintf("%.3f", p_emp_two)))

if (HAS_RIDGES) {
  library(ggridges)
  pB <- ggplot(nul_p, aes(x = beta, y = module)) +
    ggridges::geom_density_ridges(fill = "grey88", color = "white",
                                  scale = 1.05, rel_min_height = 0.01) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey55", linewidth = 0.4) +
    geom_point(data = obs_p, aes(x = beta_lmer, y = module, color = category),
               size = 2.4, inherit.aes = FALSE) +
    geom_text(data = obs_p, aes(x = beta_lmer, y = module, label = lab),
              vjust = -1.0, size = 2.5, color = "grey25", inherit.aes = FALSE) +
    scale_color_manual(values = PAL_CAT, name = NULL)
} else {
  pB <- ggplot(nul_p, aes(x = beta, y = module)) +
    geom_violin(fill = "grey88", color = NA, scale = "width") +
    geom_vline(xintercept = 0, linetype = 2, color = "grey55") +
    geom_point(data = obs_p, aes(x = beta_lmer, y = module, color = category),
               size = 2.4, inherit.aes = FALSE) +
    scale_color_manual(values = PAL_CAT, name = NULL)
}
pB <- pB +
  labs(x = expression(beta~"from 200 random genes used as predictor"), y = NULL,
       title = "B  Empirical null distribution",
       caption = "Points = observed effect; numbers = two-sided empirical P (uncorrected; none survive BH)") +
  theme_classic(base_size = 10.5) +
  theme(axis.text = element_text(color = "black"),
        plot.title = element_text(face = "bold", size = 11),
        plot.caption = element_text(hjust = 0, size = 7.5, color = "grey35"),
        legend.position = "none")

## =============================================================================
## Panel C — within-dataset 메타분석 forest (배치 강건성)
## =============================================================================
pC_df <- dat %>%
  mutate(module = factor(module, levels = ORD),
         lab_i2 = sprintf("I²=%.0f%%  %s/%s", I2, n_pos_dataset, k_dataset))
xmax_c <- max(pC_df$hi_meta) * 1.05
xr_c   <- diff(range(c(pC_df$lo_meta, pC_df$hi_meta)))

pC <- ggplot(pC_df, aes(beta_meta, module)) +
  geom_vline(xintercept = 0, color = "grey55", linewidth = 0.45) +
  geom_errorbarh(aes(xmin = lo_meta, xmax = hi_meta, color = category),
                 height = 0, linewidth = 0.65) +
  geom_point(aes(shape = sig_meta, color = category), size = 2.1,
             fill = "white", stroke = 0.7) +
  geom_text(aes(x = xmax_c + xr_c * 0.06, label = lab_i2), hjust = 0,
            size = 2.5, color = "grey30") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), guide = "none") +
  scale_color_manual(values = PAL_CAT, guide = "none") +
  coord_cartesian(xlim = c(min(pC_df$lo_meta), xmax_c), clip = "off") +
  labs(x = expression("pooled "*beta*" (inverse-variance, within-dataset)"), y = NULL,
       title = "C  Within-dataset meta-analysis",
       caption = "Right: heterogeneity I² and number of datasets with a positive estimate") +
  theme_classic(base_size = 10.5) +
  theme(axis.text = element_text(color = "black"),
        plot.title = element_text(face = "bold", size = 11),
        plot.caption = element_text(hjust = 0, size = 7.5, color = "grey35"),
        plot.margin = margin(5, 62, 5, 5))

## =============================================================================
## Panel D — null 패널 (염증 축은 "검정력 부족"이 아니라 "효과 없음")
## =============================================================================
pD_df <- dat %>% filter(category == "Inflammatory") %>%
  arrange(beta_lmer) %>%
  mutate(within = lo > -REF_EFFECT & hi < REF_EFFECT,           # CI 가 밴드 안인가
         flag   = ifelse(ctl_clean, "clean null", "depth-confounded"),
         label  = ifelse(ctl_clean, as.character(module), paste0(module, " *")),
         label  = factor(label, levels = label))

# 캡션을 데이터에서 자동 산출 → 데이터가 바뀌어도 문장이 틀리지 않음
n_in  <- sum(pD_df$within & pD_df$ctl_clean)
n_tot <- sum(pD_df$ctl_clean)
cap_D <- sprintf(paste0("Shaded band = largest specific metabolic effect (|β| = %.4f). ",
                        "%d of %d inflammatory modules with a clean negative control have ",
                        "CIs entirely within this band,\nexcluding effects of comparable ",
                        "magnitude. * = also significant for the depth control ",
                        "(not interpretable as a null)."),
                 REF_EFFECT, n_in, n_tot)

pD <- ggplot(pD_df, aes(beta_lmer, label)) +
  annotate("rect", xmin = -REF_EFFECT, xmax = REF_EFFECT,
           ymin = -Inf, ymax = Inf, fill = "#C0392B", alpha = 0.10) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey45", linewidth = 0.45) +
  geom_errorbarh(aes(xmin = lo, xmax = hi, color = flag), height = 0, linewidth = 0.75) +
  geom_point(aes(color = flag), size = 2.6, shape = 21, fill = "white", stroke = 0.8) +
  scale_color_manual(values = c(`clean null` = PAL_CAT[["Inflammatory"]],
                                `depth-confounded` = "grey55"), name = NULL) +
  labs(x = expression(beta~"(95% CI)"), y = NULL,
       title = "D  Inflammatory axes: precise nulls", caption = cap_D) +
  theme_classic(base_size = 10.5) +
  theme(axis.text = element_text(color = "black"),
        plot.title = element_text(face = "bold", size = 11),
        plot.caption = element_text(hjust = 0, size = 7.3, color = "grey35"),
        legend.position = "top", legend.margin = margin(b = -6))

## ── 저장 ────────────────────────────────────────────────────────────────────
sv <- function(p, nm, w, h) {
  ggsave(file.path(OUT_DIR, paste0(nm, ".pdf")),  p, width = w, height = h)
  ggsave(file.path(OUT_DIR, paste0(nm, ".tiff")), p, width = w, height = h,
         dpi = 300, compression = "lzw")
}
sv(pA, "Fig6A_module_forest_negctrl", 6.6, 5.4)
sv(pB, "Fig6B_empirical_null",        6.2, 5.6)
sv(pC, "Fig6C_meta_forest",           6.6, 5.0)
sv(pD, "Fig6D_null_panel",            5.6, 2.8)

if (HAS_PATCHWORK) {
  library(patchwork)
  fig <- (pA | pB) / (pC | pD) + plot_layout(heights = c(1.15, 1))
  sv(fig, "Fig6_combined", 13.0, 10.5)
  print(fig)
} else { print(pA); print(pB); print(pC); print(pD) }

cat("\n===== 모듈 판정 =====\n")
print(dat %>% arrange(desc(beta_lmer)) %>%
        select(module, category, beta_lmer, fdr_lmer, fdr_meta, I2,
               fdr_ctl, p_emp_two, verdict) %>% as.data.frame(), digits = 3)

# -----------------------------------------------------------------------------
# 기대 판정 (maxRank = 1500, n_donor = 67, k = 7 datasets)
#   robust               : PPP_REACTOME, GLYCOLYSIS      ← lmer+meta 유의, 음성대조 통과
#   specific (meta n.s.) : ROS, POLYAMINE                ← 음성대조 통과하나 메타 미달
#   depth-confounded     : OXPHOS, PAN_MYELOID           ← 음성대조에서도 유의 → 주장 불가
#   null                 : NFKB, IFN 계열, INFLAMMATORY, TRANSPORT, ARG_ENZYME, HYPOXIA, HOMEOSTATIC
#
# 논문 문장:
#  "Donor-level module regression with a depth-matched negative control identified
#   glycolysis and the pentose-phosphate pathway as the only modules that were both
#   significant and specific (negative control n.s.) and that survived within-dataset
#   meta-analysis (glycolysis I² = 0%). In contrast, OXPHOS and pan-myeloid modules were
#   also significant for the depth control and were therefore not interpreted.
#   Inflammatory modules showed precise null estimates (e.g. NF-κB β = -0.0007, P = 0.95)."
# -----------------------------------------------------------------------------
