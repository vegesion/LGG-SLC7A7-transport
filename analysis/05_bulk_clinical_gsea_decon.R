# =============================================================================
# 05_bulk_clinical_gsea_decon.R         [Fig 4]
# 임상 상관 · 상관기반 GSEA · 면역 deconvolution(+myeloid 보정 Cox)
# =============================================================================
source(here::here("R", "setup.R"))
suppressPackageStartupMessages({ library(clusterProfiler); library(msigdbr); library(ggpubr);
  library(survival); library(estimate); library(xCell); library(limma); library(fgsea);
  library(dplyr); library(ggplot2); library(devEMF)})

co <- load_cohort(TRAIN_COHORT); expr <- co$expr; clin <- co$clin; v <- cohort_gene(co)

## (a) 임상 상관
cl <- data.frame(expr = v, idh = clin$idh, grade = clin$grade, codel = clin$codel)
for (vr in c("idh", "grade", "codel")) {
  d <- cl[!is.na(cl[[vr]]), ]; if (nrow(d) < 10 || dplyr::n_distinct(d[[vr]]) < 2) next
  p <- ggpubr::ggviolin(d, x = vr, y = "expr", fill = vr, add = "boxplot",
                        add.params = list(fill = "white", width = .12)) +
    ggpubr::stat_compare_means() + ggplot2::labs(x = NULL, y = paste(GENE_OF_INTEREST, "expression")) +
    theme_paper() + ggplot2::theme(legend.position = "none")
  ggplot2::ggsave(file.path(DIR_FIGURES, paste0("Fig4A_violin_", vr, ".pdf")), p, width = 4.2, height = 4.2)
}
write_result(cl, "Fig4_clinical_association.csv")

## (b) 상관 기반 GSEA
ranked <- correlation_ranked_list(expr, gene = GENE_OF_INTEREST, method = "spearman")
write_result(run_gsea_hallmark(ranked), "Fig4_GSEA_hallmark.csv", row.names = TRUE)
write_result(run_gsego(ranked, ont = "BP"), "Fig4_GSEA_GO_BP.csv")
write_result(run_gsego(ranked, ont = "MF"), "Fig4_GSEA_GO_MF.csv")
write_result(run_gsego(ranked, ont = "CC"), "Fig4_GSEA_GO_CC.csv")



## (b-1) CAMERA

## ── Hallmark (scRNA 스크립트와 동일 버전) ──────────────────────────────────
HALLMARK <- local({
  m <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"),
                error = function(e) msigdbr(species = "Homo sapiens", category = "H"))
  split(m$gene_symbol, m$gs_name)
})

## ── 공변량 설정 ─────────────────────────────────────────────────────────────
# PRIMARY: 보정 없음 (기존 상관 기반 fgsea 와 동일한 질문 → 직접 비교 가능)
# SENSITIVITY: age + IDH 보정본도 함께 산출 (리뷰어 대비)
BULK_COVARS <- c("age", "idh")     # 없는 코호트는 자동으로 건너뜀

# =============================================================================
# 핵심 함수 — 한 코호트에서 fgsea 와 CAMERA 를 같은 모형으로 계산
# =============================================================================
#' @param expr  genes × samples, logCPM (TCGA) 또는 log2(x+1) (CGGA)
#' @param clin  data.frame, rownames = colnames(expr)
#' @param covars 보정할 임상변수명 (NULL 이면 무보정)
bulk_gsea_pair <- function(expr, clin = NULL, covars = NULL, cohort = "") {
  
  rownames(clin) <- clin$sample
  stopifnot(GENE_OF_INTEREST %in% rownames(expr))
  x <- scale(as.numeric(expr[GENE_OF_INTEREST, ]))[, 1]
  
  ## 설계행렬
  if (!is.null(covars) && !is.null(clin)) {
    covars <- intersect(covars, colnames(clin))
    cl <- clin[colnames(expr), covars, drop = FALSE]
    ok <- stats::complete.cases(cl)
    if (sum(ok) < 50) { message("  ", cohort, ": 공변량 결측 과다 → 무보정으로 진행"); covars <- NULL } #무보정 : SLC7A7만 공변량으로 넣음
  } else covars <- NULL
  
  if (is.null(covars) || !length(covars)) {
    keep <- rep(TRUE, ncol(expr))
    design <- model.matrix(~ x)
  } else {
    keep <- ok
    design <- model.matrix(as.formula(paste("~ x +", paste(covars, collapse = " + "))),
                           data = cbind(x = x, cl)[keep, , drop = FALSE])
  }
  colnames(design) <- make.names(colnames(design))
  COEF <- "x"
  
  ## self-gene 제외 (순환논리 차단)
  E <- expr[setdiff(rownames(expr), GENE_OF_INTEREST), keep, drop = FALSE]
  E <- E[matrixStats::rowVars(as.matrix(E)) > 0, , drop = FALSE]
  
  ## limma
  fit <- eBayes(lmFit(E, design))
  tt  <- topTable(fit, coef = COEF, number = Inf, sort.by = "none")
  st  <- setNames(tt$t, rownames(tt))
  
  ## (1) fgsea — 관례적 방법
  set.seed(GSEA_SEED)
  fg <- fgsea(HALLMARK, st, minSize = GSEA_MIN_SIZE, maxSize = GSEA_MAX_SIZE, nproc = 1)
  
  ## (2) CAMERA — 유전자간 상관을 데이터에서 추정
  idx <- ids2indices(HALLMARK, rownames(E))
  cm  <- camera(E, idx, design, contrast = COEF, inter.gene.cor = NA, sort = FALSE)
  
  ## (3) FRY — 회전검정 (보조)
  fr  <- fry(E, idx, design, contrast = COEF, sort = FALSE)
  
  pw <- intersect(fg$pathway, rownames(cm))
  tibble(
    cohort          = cohort,
    n               = sum(keep),
    adjusted        = if (is.null(covars) || !length(covars)) "none" else paste(covars, collapse = "+"),
    pathway         = pw,
    setSize         = fg$size[match(pw, fg$pathway)],
    fgsea_NES       = fg$NES[match(pw, fg$pathway)],
    fgsea_padj      = fg$padj[match(pw, fg$pathway)],
    camera_dir      = cm[pw, "Direction"],
    camera_p        = cm[pw, "PValue"],
    camera_FDR      = cm[pw, "FDR"],
    inter_gene_cor  = cm[pw, "Correlation"],
    fry_p           = fr[pw, "PValue"],
    fry_FDR         = fr[pw, "FDR"],
    leadingEdge     = fg$leadingEdge,
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# 실행 — 4개 코호트 × (무보정 · 보정)
# =============================================================================
coh <- load_all_cohorts()      # R/cohorts.R  (list: $expr, $clin)

res <- lapply(names(coh), function(nm) {
  message("[camera] ", nm)
  a <- bulk_gsea_pair(coh[[nm]]$expr, coh[[nm]]$clin, covars = NULL,        cohort = nm)
  b <- bulk_gsea_pair(coh[[nm]]$expr, coh[[nm]]$clin, covars = BULK_COVARS, cohort = nm)
  rbind(a, b)
}) %>% bind_rows() %>% arrange(camera_FDR)

res_csv <- res
res_csv$leadingEdge <- as.character(res_csv$leadingEdge)
write.csv(res_csv, file.path(DIR_RESULTS, "Fig4_GSEA_camera_vs_fgsea.csv"), row.names = FALSE)
saveRDS(res,       file.path(DIR_DATA_PROC, "pseudobulk_GSEA_result_bulk.rds"))

# res <- readRDS(file.path(DIR_DATA_PROC, "pseudobulk_GSEA_result_bulk.rds"))
## ── 요약: 훈련 코호트(무보정)에서 두 방법 비교 ──────────────────────────────
main <- res %>% filter(cohort == TRAIN_COHORT, adjusted == "none")
message("\n=== ", TRAIN_COHORT, " : fgsea vs CAMERA ===")
message("  fgsea padj<0.05 : ", sum(main$fgsea_padj  < 0.05, na.rm = TRUE), " / ", nrow(main))
message("  CAMERA FDR<0.05 : ", sum(main$camera_FDR  < 0.05, na.rm = TRUE))
message("  두 방법 모두    : ", sum(main$fgsea_padj < 0.05 & main$camera_FDR < 0.05, na.rm = TRUE))
message("  inter-gene cor  : median ", round(median(main$inter_gene_cor, na.rm = TRUE), 3),
        "  (max ", round(max(main$inter_gene_cor, na.rm = TRUE), 3), ")")
print(main %>% arrange(camera_FDR) %>%
        select(pathway, setSize, fgsea_NES, fgsea_padj, camera_dir, camera_p, camera_FDR, inter_gene_cor) %>%
        head(15), digits = 3)

## ── 그림: 두 방법 일치도 ────────────────────────────────────────────────────
p <- ggplot(main, aes(fgsea_NES, -log10(camera_FDR))) +
  geom_point(aes(colour = fgsea_padj < 0.05 & camera_FDR < 0.05,
                 size = inter_gene_cor), alpha = .8) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey40") +
  geom_vline(xintercept = 0, colour = "grey70") +
  scale_colour_manual(values = c(`FALSE` = "grey65", `TRUE` = PALETTE_TWO[1]),
                      name = "both significant") +
  scale_size_continuous(name = "inter-gene r", range = c(3, 12)) +
  labs(x = "fgsea NES", y = expression(-log[10]~"CAMERA FDR"),
       title = paste0(TRAIN_COHORT, " — Hallmark: fgsea vs correlation-adjusted CAMERA")) +
  theme_bw(base_size = 10)

pdf(file.path(DIR_FIGURES, "Fig4_GSEA_camera_vs_fgsea.pdf"), width = 6.5, height = 5)
p
dev.off()
save_emf(p, "Fig4_GSEA_camera_vs_fgsea.emf")
message("\n※ 논문 표기 원칙")
message("  · CAMERA 를 primary 로 선언하고 fgsea 를 관례적 참조로 병기")
message("  · inter_gene_cor 열을 반드시 함께 보고 (CAMERA 기본 프리셋은 0.01)")
message("  · scRNA 표와 동일한 열 구성 → 방법 갈아타기가 아니라 일관 보고로 읽힘")




## (c) Deconvolution + myeloid 보정
frac <- tryCatch(run_estimate(expr), error = function(e) { message("ESTIMATE 실패: ", conditionMessage(e)); NULL })
frac$sample <- gsub("\\.", "-", frac$sample)
xc <- run_xcell(expr)

# frac <- read.csv(file.path(DIR_RESULTS,"Fig4_deconvolution_scores.csv"))
# cors <- read.csv(file.path(DIR_RESULTS,"Fig4_gene_vs_fractions.csv"))

if (!is.null(xc)) frac <- if (is.null(frac)) xc else dplyr::left_join(frac, xc, by = "sample")
if (!is.null(frac)) {
  write_result(frac, "Fig4_deconvolution_scores.csv")
  cors <- correlate_gene_fractions(expr, frac); write_result(cors, "Fig4_gene_vs_fractions.csv")
  print(utils::head(cors, 15))
  mye_cols <- c("Macrophages", "Macrophages M1", "Macrophages M2", "Monocytes",
                "ImmuneScore.x", "MicroenvironmentScore")     # 광범위 지표도 포함
  sens <- dplyr::bind_rows(lapply(intersect(mye_cols, names(frac)), function(cc) {
    a <- adjust_for_fraction(expr, clin, frac, fraction_col = cc)
    g <- a[a$term == "gene", ]
    data.frame(adjusted_for = cc, HR = g$HR, lower = g$lower, upper = g$upper, p = g$p)
  }))
  print(sens)
  write_result(sens, "Fig4_sensitivity_myeloid_adjustment.csv")
  for (mye in mye_cols){
    adj <- adjust_for_fraction(expr, clin, frac, fraction_col = mye)
    print(adj); write_result(adj, paste0("Fig4_cox_adjusted_for_", mye, ".csv"))
  }
  top <- utils::head(cors[order(-abs(cors$rho)), ], 20)
  ggplot2::ggsave(file.path(DIR_FIGURES, "Fig4C_deconvolution_corr.pdf"),
    ggplot2::ggplot(top, ggplot2::aes(rho, stats::reorder(cell_type, rho), fill = FDR < 0.05)) +
      ggplot2::geom_col() +
      ggplot2::scale_fill_manual(values = c(`TRUE` = "#D64B4B", `FALSE` = "grey75")) +
      ggplot2::labs(x = "Spearman rho", y = NULL) + theme_paper(), width = 6, height = 6)
  
  p <- ggplot2::ggplot(top, ggplot2::aes(rho, stats::reorder(cell_type, rho), fill = FDR < 0.05)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#FC8D59", `FALSE` = "grey75")) +
    ggplot2::labs(x = "Spearman rho", y = NULL) + theme_paper()
  save_emf(p + theme_blank_frame(), "Fig4C_deconvolution_corr.emf", width = 2)

}
message("완료: 임상/GSEA/deconvolution")


# 시각화(visualization) ------------------------------------------------------

res_adj <- res %>%
  filter(adjusted != "none")

write.csv(res_adj, "gsea_mixed.csv")
dim(res_adj)

res_lgg <- res_adj %>%
  filter(cohort == "TCGA_LGG") %>%
  filter(camera_FDR <= 0.05)

res_gbm <- res_adj %>%
  filter(cohort == "TCGA_GBM") %>%
  filter(camera_FDR <= 0.05)

res_325 <- res_adj %>%
  filter(cohort == "CGGA_325") %>%
  filter(camera_FDR <= 0.05)

res_693 <- res_adj %>%
  filter(cohort == "CGGA_693") %>%
  filter(camera_FDR <= 0.05)



suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr)
  library(forcats); library(ggplot2); library(scales); library(tibble)
})


# =============================================================================
# 0.  설정 — 여기만 수정하면 됩니다
# =============================================================================

## 본인 객체 이름으로 교체. 이름(왼쪽)이 그림의 x축 라벨 겸 순서가 됩니다.
df_list <- list(
  "TCGA-LGG" = res_lgg,
  "TCGA-GBM" = res_gbm,
  "CGGA-693" = res_693,
  "CGGA-325" = res_325
)

FDR_CUT     <- 0.05   # 테두리를 진하게 칠할 CAMERA FDR 기준
MIN_COHORT  <- 2      # 최소 몇 개 코호트에서 유의해야 그림에 넣을지 (1 = 합집합 전부)
SIZE_CAP    <- 8      # -log10(FDR) 상한. 소수 극단값이 점 크기를 독식하는 것 방지
OUT_DIR     <- "figures"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 0.  MSigDB Hallmark systematic name (M-number) 룩업
#     msigdbr 10.x 는 category -> collection 으로 인자명이 바뀌어서 둘 다 대응
# =============================================================================
get_hallmark_ids <- function(species = "Homo sapiens") {
  args <- names(formals(msigdbr::msigdbr))
  df <- if ("collection" %in% args) {
    msigdbr::msigdbr(species = species, collection = "H")
  } else {
    msigdbr::msigdbr(species = species, category = "H")
  }
  df %>%
    dplyr::distinct(gs_name, gs_id) %>%
    tibble::deframe()          # named vector: HALLMARK_XXX -> "M5890"
}

HALLMARK_ID <- get_hallmark_ids()
# 한 번 뽑아서 고정하고 싶으면:  dput(HALLMARK_ID)  결과를 스크립트에 붙여넣으세요

# =============================================================================
# 1.  4개 df 결합 + 정리
# =============================================================================

NUM_COLS <- c("setSize", "fgsea_NES", "fgsea_padj",
              "camera_p", "camera_FDR", "inter_gene_cor",
              "fry_p", "fry_FDR")

dat <- bind_rows(
  lapply(df_list, function(d) {
    d %>%
      # leadingEdge는 list-column이라 downstream에서 걸림. cohort는 .id로 다시 붙음
      select(-any_of(c("leadingEdge", "stringsAsFactors", "cohort"))) %>%
      mutate(across(any_of(NUM_COLS), as.numeric))
  }),
  .id = "cohort"
) %>%
  mutate(cohort = factor(cohort, levels = names(df_list)))

stopifnot(all(c("pathway", "fgsea_NES", "camera_p", "camera_FDR") %in% names(dat)))

## 안전장치: 혹시 Down이 섞여 있으면 알려줌 (현재 데이터는 110/110 Up)
if ("camera_dir" %in% names(dat)) {
  n_down <- sum(dat$camera_dir == "Down", na.rm = TRUE)
  message(sprintf("camera_dir: Up %d / Down %d", sum(dat$camera_dir == "Up", na.rm = TRUE), n_down))
  if (n_down > 0)
    warning("Down이 존재합니다. 아래 scale_fill_gradientn을 발산형으로 바꾸는 것을 고려하세요.")
}


# =============================================================================
# 2.  라벨 정리 (HALLMARK_ 제거 + 약어 복원)
# =============================================================================

pretty_label <- function(x, with_id = TRUE, id_sep = " ") {
  base <- x %>%
    str_remove("^HALLMARK_") %>%
    str_replace_all("_", " ") %>%
    str_to_sentence()
  
  fixes <- c(
    "Tnfa signaling via nfkb"           = "TNF\u03b1 signaling via NF-\u03baB",
    "Interferon gamma response"         = "IFN-\u03b3 response",
    "Interferon alpha response"         = "IFN-\u03b1 response",
    "Il6 jak stat3 signaling"           = "IL6\u2013JAK\u2013STAT3 signaling",
    "Il2 stat5 signaling"               = "IL2\u2013STAT5 signaling",
    "Epithelial mesenchymal transition" = "Epithelial\u2013mesenchymal transition",
    "Kras signaling up"                 = "KRAS signaling (up)",
    "Kras signaling dn"                 = "KRAS signaling (down)",
    "P53 pathway"                       = "p53 pathway",
    "Tgf beta signaling"                = "TGF-\u03b2 signaling",
    "Pi3k akt mtor signaling"           = "PI3K\u2013AKT\u2013mTOR signaling",
    "Mtorc1 signaling"                  = "mTORC1 signaling",
    "Wnt beta catenin signaling"        = "WNT\u2013\u03b2-catenin signaling",
    "Reactive oxygen species pathway"   = "Reactive oxygen species",
    "Myc targets v1"                    = "MYC targets V1",
    "Myc targets v2"                    = "MYC targets V2",
    "E2f targets"                       = "E2F targets",
    "G2m checkpoint"                    = "G2/M checkpoint",
    "Oxidative phosphorylation"         = "Oxidative phosphorylation",
    "Dna repair"                        = "DNA repair",
    "Uv response up"                    = "UV response (up)",
    "Uv response dn"                    = "UV response (down)",
    "Unfolded protein response"         = "Unfolded protein response",
    "Il6 jak stat3"                     = "IL6\u2013JAK\u2013STAT3"
  )
  lab <- unname(ifelse(base %in% names(fixes), fixes[base], base))
  if (!with_id) return(lab)
  
  id <- unname(HALLMARK_ID[x])          # 매칭 안 되면 NA
  ifelse(is.na(id), lab, paste0(lab, id_sep, "(", id, ")"))
}


# =============================================================================
# 3.  기능적 카테고리 (facet 행 그룹)
#     — 면역 블록이 빽빽하고 대사 블록이 성긴 대비가 곧 메시지입니다
# =============================================================================

CAT_MAP <- list(
  "Immune / inflammatory" = c(
    "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",   "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING",   "HALLMARK_IL2_STAT5_SIGNALING",
    "HALLMARK_COMPLEMENT",                "HALLMARK_ALLOGRAFT_REJECTION"),
  
  "Stromal / hypoxic" = c(
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_ANGIOGENESIS",
    "HALLMARK_HYPOXIA",        "HALLMARK_COAGULATION",
    "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_APICAL_JUNCTION",
    "HALLMARK_APICAL_SURFACE",  "HALLMARK_MYOGENESIS"),
  
  "Signaling / stress" = c(
    "HALLMARK_KRAS_SIGNALING_UP", "HALLMARK_P53_PATHWAY", "HALLMARK_APOPTOSIS",
    "HALLMARK_PI3K_AKT_MTOR_SIGNALING", "HALLMARK_MTORC1_SIGNALING",
    "HALLMARK_NOTCH_SIGNALING",  "HALLMARK_UV_RESPONSE_UP",
    "HALLMARK_DNA_REPAIR",       "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
    "HALLMARK_PROTEIN_SECRETION","HALLMARK_ESTROGEN_RESPONSE_EARLY",
    "HALLMARK_ESTROGEN_RESPONSE_LATE", "HALLMARK_ANDROGEN_RESPONSE"),
  
  "Metabolic / proliferative" = c(
    "HALLMARK_XENOBIOTIC_METABOLISM", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
    "HALLMARK_ADIPOGENESIS", "HALLMARK_PEROXISOME",
    "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_GLYCOLYSIS",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_HEME_METABOLISM",
    "HALLMARK_BILE_ACID_METABOLISM", "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
    "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2",
    "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_MITOTIC_SPINDLE")
)
CAT_LEVELS <- c(names(CAT_MAP), "Other")
cat_lookup <- stack(setNames(CAT_MAP, names(CAT_MAP)))  # values = pathway, ind = category

assign_cat <- function(p) {
  out <- as.character(cat_lookup$ind)[match(p, as.character(cat_lookup$values))]
  factor(ifelse(is.na(out), "Other", out), levels = CAT_LEVELS)
}


# =============================================================================
# 4.  행 순서 결정: 카테고리 → 재현 코호트 수 → NES 중앙값
# =============================================================================

plot_dat <- dat %>%
  mutate(category = assign_cat(pathway),
         label    = pretty_label(pathway),
         pass_fdr = camera_FDR <= FDR_CUT,
         size_val = pmin(-log10(camera_FDR), SIZE_CAP)) %>%
  group_by(pathway) %>%
  mutate(n_cohort = n_distinct(cohort),
         med_nes  = median(fgsea_NES, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(n_cohort >= MIN_COHORT)

## y축은 아래에서 위로 쌓이므로 오름차순 정렬 → 재현성 높은 것이 각 facet 상단
lev <- plot_dat %>%
  distinct(label, category, n_cohort, med_nes) %>%
  arrange(category, n_cohort, med_nes) %>%
  pull(label)

plot_dat <- plot_dat %>% mutate(label = factor(label, levels = lev))

message(sprintf("표시 pathway %d개 / 점 %d개 (MIN_COHORT = %d)",
                nlevels(plot_dat$label), nrow(plot_dat), MIN_COHORT))


# =============================================================================
# 5.  격자 완성 — 모든 (pathway x cohort) 칸을 강제로 생성
#
#     Hallmark 50개는 4개 코호트 전부에서 검정되었으므로 "검정 안 됨"인 칸은
#     존재하지 않습니다. 따라서 점이 없는 칸 = camera_p > 0.05, 예외 없음.
#     complete()로 만들면 res와 이름을 대조할 필요가 없어 누락이 원천 차단됩니다.
# =============================================================================

grid_dat <- plot_dat %>%
  select(cohort, pathway, fgsea_NES, camera_p, camera_FDR, size_val, pass_fdr) %>%
  complete(pathway, cohort) %>%                    # 빠진 조합을 NA 행으로 생성
  mutate(category = assign_cat(pathway),
         label    = factor(pretty_label(pathway), levels = lev),
         is_sig   = !is.na(camera_p))

dot_dat <- filter(grid_dat,  is_sig)   # 점 (camera_p <= 0.05)
bg_dat  <- filter(grid_dat, !is_sig)   # X  (camera_p >  0.05)

## 검산: 칸 수 = 점 + X 가 반드시 성립해야 합니다
n_cell <- nlevels(droplevels(grid_dat$label)) * nlevels(grid_dat$cohort)
message(sprintf("격자 %d칸 = 점 %d + X %d  %s",
                n_cell, nrow(dot_dat), nrow(bg_dat),
                if (nrow(grid_dat) == n_cell) "[OK]" else "[!! 불일치]"))
stopifnot(nrow(grid_dat) == n_cell)
# =============================================================================
# 6.  Fig 4A — dotplot
# =============================================================================

theme_pub <- function(base_size = 9) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major   = element_line(linewidth = 0.25, colour = "grey90"),
      panel.grid.minor   = element_blank(),
      panel.border       = element_rect(linewidth = 0.4, colour = "grey40"),
      panel.spacing.y    = unit(2.5, "pt"),
      strip.background   = element_rect(fill = "grey94", colour = "grey40", linewidth = 0.4),
      # strip.text.y.left  = element_text(angle = 0, hjust = 0, face = "bold", size = base_size - 0.2),
      strip.text.y.left  = element_blank(),
      axis.text.x        = element_text(face = "bold", size = base_size),
      axis.text.y        = element_text(size = base_size - 0.5, colour = "grey15"),
      axis.title         = element_blank(),
      axis.ticks         = element_line(linewidth = 0.3),
      legend.key         = element_blank(),
      legend.background  = element_blank(),
      legend.title       = element_text(size = base_size - 0.5, face = "bold"),
      legend.text        = element_text(size = base_size - 1),
      legend.box.spacing = unit(4, "pt"),
      plot.title         = element_text(face = "bold", size = base_size + 2),
      plot.subtitle      = element_text(size = base_size - 0.5, colour = "grey30")
    )
}

p4a <- ggplot(plot_dat, aes(x = cohort, y = label))

if (!is.null(bg_dat) && nrow(bg_dat) > 0) {
  p4a <- p4a + geom_point(data = bg_dat, shape = 4, size = 3,
                          colour = "grey82", stroke = 0.45)
}

p4a <- p4a +
  geom_point(aes(size = size_val, fill = fgsea_NES, colour = pass_fdr),
             shape = 21, stroke = 0) +
  
  facet_grid(rows = vars(category), scales = "free_y", space = "free_y", switch = "y", labeller = labeller(category = function(x) rep("", length(x)))) +
  
  ## 전부 Up이므로 순차형 스케일. Down이 생기면 여기를 발산형으로 교체
  scale_fill_gradientn(
    colours = c("#FEF0D9", "#FDCC8A", "#FC8D59", "#E34A33", "#8C2D04"),
    limits  = c(1.3, 3.8), oob = scales::squish,
    name    = "fgsea NES",
    guide   = guide_colourbar(barwidth = unit(0.4, "cm"), barheight = unit(3, "cm"),
                              frame.colour = "grey70", ticks.colour = "grey40", order = 1)) +
  
  scale_size_continuous(
    range  = c(1.6, 12),
    limits = c(0, SIZE_CAP),
    breaks = c(1.3, 3, 5, 8),
    labels = c("1.3", "3", "5", "\u22658"),
    name   = expression(-log[10]~italic(FDR)[CAMERA]),
    guide  = guide_legend(override.aes = list(fill = "grey75", colour = "grey25"), order = 2)) +
  
  scale_colour_manual(
    values = c("TRUE" = "grey12", "FALSE" = "grey72"),
    labels = c("TRUE" = "< 0.05", "FALSE" = "> 0.05"),
    name   = "CAMERA FDR",
    guide  = guide_legend(override.aes = list(size = 4, fill = "grey85"), order = 3)) +
  
  scale_x_discrete(position = "top", expand = expansion(add = 0.65)) +
  scale_y_discrete(expand = expansion(add = 0.65)) +
  theme_pub()


print(p4a)
print(p4a+theme_blank_frame(keep_grid_x = TRUE) + theme(
  panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
  panel.grid.minor = element_line(colour = "grey92", linewidth = 0.2)
))

save_emf(p4a, "fig4.emf", height = 9.5, width = 3.5, legend_p = "right", draw_minor_grid = TRUE, draw_major_grid = TRUE)
# =============================================================================
# 7.  저장 — 높이는 행 수에 맞춰 자동
# =============================================================================

n_rows <- nlevels(plot_dat$label)
n_cats <- n_distinct(plot_dat$category)
fig_h  <- 1.5 + 0.185 * n_rows + 0.12 * n_cats   # inch
fig_w  <- 6.7                                    # ~170 mm, 2단 저널 full width

ggsave(file.path(OUT_DIR, "Fig4A_camera_dotplot.pdf"), p4a,
       width = fig_w, height = fig_h, device = cairo_pdf)   # cairo: γ, κ, – 글리프 보존

ggsave(file.path(OUT_DIR, "Fig4A_camera_dotplot.tiff"), p4a,
       width = fig_w, height = fig_h, dpi = 300, compression = "lzw")

message(sprintf("saved: %s (%.1f x %.1f in)", OUT_DIR, fig_w, fig_h))


# =============================================================================
# 8.  Fig 4B (선택) — fgsea가 잡고 CAMERA가 놓친 것
#     res(전체 adjusted 행)가 있어야 실행됩니다.
#     "OXPHOS는 상관구조 보정 후 살아남지 못한다"를 한 패널로 보여주는 용도.
# =============================================================================

if (exists("res")) {
  
  cmp <- res %>%
    select(-any_of(c("leadingEdge", "stringsAsFactors"))) %>%
    mutate(across(any_of(NUM_COLS), as.numeric),
           cohort = factor(gsub("_", "-", as.character(cohort)),
                           levels = levels(plot_dat$cohort)),
           status = case_when(
             fgsea_padj <= 0.05 & camera_p <= 0.05 ~ "Both",
             fgsea_padj <= 0.05 & camera_p >  0.05 ~ "fgsea only",
             TRUE                                  ~ "Neither"),
           status = factor(status, c("Both", "fgsea only", "Neither")),
           label  = pretty_label(pathway)) %>%
    filter(!is.na(cohort))
  
  HILITE <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_MYC_TARGETS_V1",
              "HALLMARK_MTORC1_SIGNALING", "HALLMARK_E2F_TARGETS",
              "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_GLYCOLYSIS")
  
  p4b <- ggplot(cmp, aes(-log10(fgsea_padj), -log10(camera_p))) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted",
                colour = "grey65", linewidth = 0.35) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               colour = "grey45", linewidth = 0.35) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed",
               colour = "grey45", linewidth = 0.35) +
    geom_point(aes(fill = inter_gene_cor, shape = status),
               size = 2.1, colour = "grey25", stroke = 0.35, alpha = 0.9) +
    facet_wrap(~ cohort, nrow = 1) +
    scale_shape_manual(values = c("Both" = 21, "fgsea only" = 24, "Neither" = 22),
                       name = NULL) +
    scale_fill_viridis_c(option = "mako", direction = -1,
                         name = "Inter-gene\ncorrelation") +
    labs(x = expression(-log[10]~italic(P)[adj]~"(fgsea)"),
         y = expression(-log[10]~italic(P)~"(CAMERA)")) +
    theme_pub() +
    theme(axis.title = element_text(size = 8.5))
  
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p4b <- p4b + ggrepel::geom_text_repel(
      data = filter(cmp, pathway %in% HILITE, fgsea_padj <= 0.05),
      aes(label = label), size = 2.2, colour = "grey15",
      min.segment.length = 0, segment.size = 0.25,
      box.padding = 0.3, max.overlaps = 30)
  } else {
    message("ggrepel 미설치 — Fig4B 라벨 생략")
  }
  
  print(p4b)
  ggsave(file.path(OUT_DIR, "Fig4B_fgsea_vs_camera.pdf"), p4b,
         width = 7.2, height = 2.6, device = cairo_pdf)
}


# =============================================================================
# 9.  캡션에 쓸 수치 자동 산출
# =============================================================================

summ <- dat %>%
  group_by(cohort) %>%
  summarise(n_p05   = sum(camera_p   <= 0.05, na.rm = TRUE),
            n_fdr05 = sum(camera_FDR <= FDR_CUT, na.rm = TRUE),
            .groups = "drop")

core_p   <- dat %>% count(pathway) %>% filter(n == n_distinct(dat$cohort)) %>% pull(pathway)
core_fdr <- dat %>% filter(camera_FDR <= FDR_CUT) %>% count(pathway) %>%
  filter(n == n_distinct(dat$cohort)) %>% pull(pathway)

cat("\n--- caption numbers ---\n")
print(as.data.frame(summ))
cat(sprintf("\ncamera_p <= 0.05 in all %d cohorts : %d sets\n",
            n_distinct(dat$cohort), length(core_p)))
cat(sprintf("camera_FDR <= %.2f in all %d cohorts: %d sets\n",
            FDR_CUT, n_distinct(dat$cohort), length(core_fdr)))
cat("\ncore (FDR, all cohorts):\n"); cat(paste0("  ", core_fdr, collapse = "\n"), "\n")

if ("camera_dir" %in% names(dat))
  cat(sprintf("\ndirection: %d/%d Up\n",
              sum(dat$camera_dir == "Up", na.rm = TRUE), nrow(dat)))


