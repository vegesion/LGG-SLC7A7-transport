# =============================================================================
# 05_bulk_clinical_gsea_decon.R         [Fig 4]
# 임상 상관 · 상관기반 GSEA · 면역 deconvolution(+myeloid 보정 Cox)
# =============================================================================
source(here::here("R", "setup.R"))
suppressPackageStartupMessages({ library(clusterProfiler); library(msigdbr); library(ggpubr);
  library(survival); library(estimate); library(xCell); library(limma); library(fgsea);
  library(dplyr); library(ggplot2)})

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
}
message("완료: 임상/GSEA/deconvolution")
