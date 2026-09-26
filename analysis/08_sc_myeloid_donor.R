# =============================================================================
# 08_sc_myeloid_donor.R                  [Fig 6]
# 아종 UMAP · donor 수준 발현 vs 아종 비율 회귀(핵심) · MG vs MDM · pseudobulk DEG/GSEA
# =============================================================================

# =============================================================================
# 08_add_microglia_camera.R
#   → 08_sc_myeloid_donor.R 뒤에 추가 (또는 이 블록으로 pseudobulk 파트를 대체)
# -----------------------------------------------------------------------------
# 설계: microglia subset → donor pseudobulk → depth(nCount)+dataset 보정
#       → limma 유전자 수준 DEG (★ Fig6 주 결과)
#       → CAMERA / FRY / fgsea 병기 (bulk 표와 동일 형식)
#
# 원칙 (Methods 에 그대로 쓸 수 있게 코드에 반영)
#   1) 분석 단위 = donor (세포 단위 검정 금지)
#   2) 발현 = 연속값 (이분화 금지)
#   3) 공변량 = nCount_RNA(depth) + dataset
#   4) self-gene(SLC7A7) 은 모든 랭킹·집합에서 제외
#   5) gene set 검정은 CAMERA 를 primary (fgsea 는 관례적 병기)
#
# 출력
#   results/Fig6_microglia_DEG.csv               ★ Fig6 주 결과 (유전자 수준)
#   results/Fig6_microglia_GSEA_camera_vs_fgsea.csv  ★ 논문 Table (scRNA)
#   results/Supp_microglia_ICC_depth.csv         Methods 용 숫자 1줄
#   figures/Fig6_microglia_volcano.pdf           ★ Fig6 그림
# =============================================================================
source(here::here("R", "setup.R"))
suppressPackageStartupMessages({
  library(Matrix); library(edgeR); library(limma); library(fgsea); library(msigdbr)
  library(lme4); library(dplyr); library(ggplot2); library(ggrepel)
})


mye <- readRDS(file.path(DIR_DATA_PROC, "myeloid.rds"))
mg  <- readRDS(file.path(DIR_DATA_PROC, "microglia.rds"))
sub_col <- intersect(c("annotation_level_4","celltype_original","annotation_level_3"),
                     colnames(mye@meta.data))[1]

ggplot2::ggsave(file.path(DIR_FIGURES, "Fig6A_umap_myeloid_subtype.pdf"),
                DimPlot(mye, group.by = sub_col, label = TRUE, raster = TRUE) + theme_paper(), width = 8, height = 6)

## donor 수준 발현 vs 아종 비율 (depth 보정 + dataset 혼합모형)
pr <- donor_subtype_proportion(mye, sub_col)
print(pr$result); write_result(pr$result, "Fig6_donor_subtype_proportion_regression.csv")
write_result(pr$data, "Fig6_donor_subtype_proportions.csv")
top_sub <- utils::head(pr$result$subtype[order(pr$result$p)], 4)
ggplot2::ggsave(file.path(DIR_FIGURES, "Fig6B_proportion_regression.pdf"),
                ggplot2::ggplot(dplyr::filter(pr$data, subtype %in% top_sub), ggplot2::aes(gene, prop)) +
                  ggplot2::geom_point(size = 1.8, alpha = .8) + ggplot2::geom_smooth(method = "lm") +
                  ggpubr::stat_cor(method = "spearman", size = 3) +
                  ggplot2::facet_wrap(~subtype, scales = "free_y") +
                  ggplot2::labs(x = paste(GENE_OF_INTEREST, "(donor mean)"), y = "Subtype proportion") +
                  theme_paper(), width = 7.5, height = 6)

## microglia vs MDM
cmp <- compare_mg_mdm(mye); print(cmp$test); write_result(cmp$data, "Fig6_microglia_vs_MDM.csv")
ggplot2::ggsave(file.path(DIR_FIGURES, "Fig6C_mg_vs_mdm.pdf"),
                ggpubr::ggpaired(tidyr::pivot_longer(cmp$data, c(Microglia, MDM),
                                                     names_to = "lineage", values_to = "expr"),
                                 x = "lineage", y = "expr", id = "donor", line.color = "grey80") +
                  ggpubr::stat_compare_means(paired = TRUE) +
                  ggplot2::labs(x = NULL, y = paste(GENE_OF_INTEREST, "(donor mean)")) + theme_paper(),
                width = 4, height = 4.5)


# graph -------------------------------------------------------------------

p <- ggpubr::ggpaired(tidyr::pivot_longer(cmp$data, c(Microglia, MDM),
                                          names_to = "lineage", values_to = "expr"),
                      x = "lineage", y = "expr", id = "donor", line.color = "grey80") +
  ggpubr::stat_compare_means(paired = TRUE, geom = "blank") +
  ggplot2::labs(x = NULL, y = paste(GENE_OF_INTEREST, "(donor mean)")) + theme_paper()


save_emf(p + theme_blank_frame(), "Fig6_microglia_vs_MDM.emf")
# ── 파라미터 ────────────────────────────────────────────────────────────────
MG_LABELS <- c("TAM-MG aging sig", "TAM-MG pro-infl I",
               "TAM-MG pro-infl II", "TAM-MG prolif")
LABEL_COL        <- "annotation_level_4"
DEPTH_COL        <- "nCount_RNA"      # ★ nFeature 아님 (UMI 총량)
MIN_DONOR_CELLS  <- 50                # microglia 는 donor당 세포수가 적음
MIN_DATASET_DONOR<- 3

HALLMARK <- local({
  m <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"),
                error = function(e) msigdbr(species = "Homo sapiens", category = "H"))
  split(m$gene_symbol, m$gs_name)
})

GO_BP <- local({
  m <- tryCatch(
    msigdbr(
      species = "Homo sapiens",
      collection = "C5",
      subcollection = "GO:BP"
    ),
    error = function(e) {
      msigdbr(
        species = "Homo sapiens",
        category = "C5",
        subcategory = "GO:BP"
      )
    }
  )
  split(m$gene_symbol, m$gs_name)
})

# =============================================================================
# §1. microglia subset → donor pseudobulk
# =============================================================================
mye <- readRDS(file.path(DIR_DATA_PROC, "myeloid.rds"))
stopifnot(LABEL_COL %in% colnames(mye@meta.data), DEPTH_COL %in% colnames(mye@meta.data))

mg <- mye[, mye@meta.data[[LABEL_COL]] %in% MG_LABELS]
message(sprintf("[1] microglia cells = %d", ncol(mg)))

md <- mg@meta.data
donor_col <- if ("donor_id" %in% colnames(md)) "donor_id" else "donor"

cnt <- SeuratObject::LayerData(mg, assay = "RNA", layer = "counts")
dat <- SeuratObject::LayerData(mg, assay = "RNA", layer = "data")   # LogNormalize

donor_f <- factor(as.character(md[[donor_col]]))
ok      <- donor_f %in% names(which(table(donor_f) >= MIN_DONOR_CELLS))
cnt <- cnt[, ok, drop = FALSE]; dat <- dat[, ok, drop = FALSE]
md  <- md[ok, , drop = FALSE]; donor_f <- droplevels(donor_f[ok])

# 희소 지시행렬로 집계 (AggregateExpression 버전 버그 회피)
ind <- Matrix::sparse.model.matrix(~ 0 + donor_f); colnames(ind) <- levels(donor_f)
ncell <- Matrix::colSums(ind)

pb <- as.matrix(cnt %*% ind)                              # pseudobulk counts
dm <- sweep(as.matrix(dat %*% ind), 2, ncell, "/")        # donor 평균 정규화 발현

donor_meta <- md %>%
  mutate(donor = as.character(.data[[donor_col]])) %>%
  group_by(donor) %>%
  summarise(depth   = mean(.data[[DEPTH_COL]]),
            dataset = names(which.max(table(author))),
            n_cell  = dplyr::n(), .groups = "drop") %>% as.data.frame()
rownames(donor_meta) <- donor_meta$donor
donor_meta <- donor_meta[colnames(pb), ]

# donor 1명뿐인 dataset 제외 (더미가 자유도 소진)
k2 <- donor_meta$dataset %in% names(which(table(donor_meta$dataset) >= MIN_DATASET_DONOR))
pb <- pb[, k2, drop = FALSE]; dm <- dm[, k2, drop = FALSE]
donor_meta <- donor_meta[k2, ]; donor_meta$dataset <- factor(donor_meta$dataset)

message(sprintf("    donors = %d | datasets = %d | median cells/donor = %.0f",
                nrow(donor_meta), nlevels(donor_meta$dataset), median(donor_meta$n_cell)))

# =============================================================================
# 08_sc_myeloid_donor.R  —  §1b ~ §5b 교체 블록
# -----------------------------------------------------------------------------
# 변경 핵심 : 예측변수 x = donor 평균 LogNormalize  →  pseudobulk TMM log-CPM
# 전제      : 기존 §1 을 실행해서 pb, dm, donor_meta 가 메모리에 있어야 함
#             (기존 §2 §3 §5 는 이 블록으로 통째 대체)
#
# 이 블록이 고치는 것
#   1) filterByExpr 이 design(=x 포함) 에 의존하던 순환  → null design 으로 분리
#   2) 예측변수 QC (0-count donor, library size 상관) 강제 출력
#   3) donor 세포수 편차를 sample weight 로 반영
#   4) depth-only "음성대조" → 발현량·검출률 매칭 대조 유전자 기반 empirical p
#   5) 헤더에 약속돼 있었지만 코드에 없던 volcano / Fig6_microglia_DEG.csv 생성
# =============================================================================

USE_SAMPLE_WEIGHTS <- TRUE   # donor당 세포수 편차가 크면 TRUE 권장
PRIOR_COUNT        <- 1      # GOI count 가 낮으면 2~3 까지 올려도 됨
N_CTRL             <- 200    # §5b 매칭 대조 유전자 개수

# =============================================================================
# §1b. 필터링 · TMM  ( x 와 무관한 null design 사용 → 순환 제거 )
# =============================================================================
stopifnot(GENE_OF_INTEREST %in% rownames(pb),
          identical(colnames(pb), rownames(donor_meta)))

design0 <- model.matrix(~ depth + dataset, data = donor_meta)   # x 미포함
colnames(design0) <- make.names(colnames(design0))

dge  <- DGEList(counts = pb)
keep <- filterByExpr(dge, design = design0, min.count = PB_FILTER_MIN_COUNT)
stopifnot(GENE_OF_INTEREST %in% names(keep))           # 이름 없으면 append 돼서 조용히 깨짐
message(sprintf("[1b] %s 가 filterByExpr 를 자력 통과: %s",
                GENE_OF_INTEREST, keep[[GENE_OF_INTEREST]]))
keep[GENE_OF_INTEREST] <- TRUE                          # 예측변수는 강제 보존

dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
message(sprintf("     genes = %d | donors = %d", nrow(dge), ncol(dge)))

# =============================================================================
# §1c. 예측변수 x = TMM log-CPM
# =============================================================================
lcpm  <- edgeR::cpm(dge, log = TRUE, prior.count = PRIOR_COUNT)
x_raw <- lcpm[GENE_OF_INTEREST, ]

# self-gene 을 library size 에서 뺀 버전 (compositional sanity check)
cnt_goi <- as.numeric(pb[GENE_OF_INTEREST, colnames(dge)])
eff_lib <- dge$samples$lib.size * dge$samples$norm.factors
x_excl  <- log2((cnt_goi + PRIOR_COUNT) / ((eff_lib - cnt_goi) / 1e6))

# =============================================================================
# §1d. 예측변수 QC  ★ 이 출력 안 보고 넘어가면 안 됨
# =============================================================================
sr <- function(a, b) stats::cor(a, b, method = "spearman")

qc <- data.frame(
  donor                = colnames(dge),
  count_goi            = cnt_goi,
  lib_size             = dge$samples$lib.size,
  norm_factor          = dge$samples$norm.factors,
  n_cell               = donor_meta$n_cell,
  depth                = donor_meta$depth,
  dataset              = donor_meta$dataset,
  x_logcpm             = as.numeric(x_raw),
  x_logcpm_excl_self   = x_excl,
  x_donor_mean_lognorm = as.numeric(dm[GENE_OF_INTEREST, colnames(dge)])
)
write.csv(qc, file.path(DIR_RESULTS, "Supp_predictor_QC.csv"), row.names = FALSE)

message(sprintf(
  "[1d] %s counts/donor  min=%d  median=%.0f  zero=%d/%d
     rho(x, lib_size)=%.2f   rho(x, depth)=%.2f   rho(x, n_cell)=%.2f
     rho(x, x_excl_self)=%.3f   rho(x, 기존 donor-mean x)=%.3f",
  GENE_OF_INTEREST, min(qc$count_goi), stats::median(qc$count_goi),
  sum(qc$count_goi == 0), nrow(qc),
  sr(qc$x_logcpm, qc$lib_size), sr(qc$x_logcpm, qc$depth), sr(qc$x_logcpm, qc$n_cell),
  sr(qc$x_logcpm, qc$x_logcpm_excl_self), sr(qc$x_logcpm, qc$x_donor_mean_lognorm)))

if (any(qc$count_goi == 0))
  warning("count=0 donor 존재 → 해당 donor 의 x 는 사실상 -log2(lib.size). ",
          "MIN_DONOR_CELLS 를 올리거나 그 donor 를 제외할 것.")
if (abs(sr(qc$x_logcpm, qc$lib_size)) > 0.5)
  warning("x 가 library size 와 강하게 상관 → depth 보정 후에도 잔여 교란 가능.")
if (sr(qc$x_logcpm, qc$x_donor_mean_lognorm) < 0.8)
  warning("두 예측변수 정의가 크게 불일치 → 어느 쪽을 쓰느냐로 결론이 바뀔 수 있음. ",
          "논문에는 둘 다 보고할 것.")

donor_meta$x <- as.numeric(scale(x_raw))   # 단위: SLC7A7 log-CPM 1 SD

# =============================================================================
# §2. Methods 용 숫자 — ICC(dataset) 와 depth 상관   (새 x 기준으로 재계산)
# =============================================================================
icc_df  <- data.frame(donor = colnames(dge), expr = as.numeric(x_raw),
                      depth = donor_meta$depth, dataset = donor_meta$dataset)
fit_icc <- lme4::lmer(expr ~ depth + (1 | dataset), data = icc_df, REML = TRUE)
vc  <- as.data.frame(lme4::VarCorr(fit_icc))
ICC <- vc$vcov[vc$grp == "dataset"] / sum(vc$vcov)

icc_out <- data.frame(
  n_donor = nrow(icc_df),
  n_dataset = length(unique(icc_df$dataset)),
  ICC_dataset = ICC,
  rho_expr_depth = sr(icc_df$expr, icc_df$depth),
  rho_expr_ncell = sr(icc_df$expr, donor_meta$n_cell),
  depth_var_by_dataset_R2 = summary(lm(depth ~ dataset, icc_df))$r.squared
)
write.csv(icc_out, file.path(DIR_RESULTS, "Supp_microglia_ICC_depth.csv"), row.names = FALSE)
print(icc_out, digits = 3)

if (icc_out$n_dataset < 5)
  message("  ! dataset 수 < 5 → ICC 는 점추정 의미 없음. 서술적으로만 보고할 것.")
if (icc_out$depth_var_by_dataset_R2 > 0.8)
  warning("depth 가 dataset 더미로 거의 설명됨(R2>0.8) → 둘이 준공선. ",
          "depth 를 dataset 내 centering 한 버전으로 대체 검토.")

# =============================================================================
# §3. limma 유전자 수준 DEG  ★ Fig6 주 결과
# =============================================================================
design <- model.matrix(~ x + depth + dataset, data = donor_meta)
colnames(design) <- make.names(colnames(design))
COEF <- "x"

stopifnot(COEF %in% colnames(design))
if (qr(design)$rank < ncol(design))
  stop("design 이 rank-deficient. MIN_DATASET_DONOR 를 올리거나 depth 를 빼야 함.")
message(sprintf("[3] residual df = %d  (donor %d - 파라미터 %d)",
                ncol(dge) - ncol(design), ncol(dge), ncol(design)))
if (ncol(dge) - ncol(design) < 8)
  warning("잔차 df < 8 → eBayes 로도 버티기 어려움. dataset 더미 수를 줄일 것.")

v <- if (USE_SAMPLE_WEIGHTS) {
  limma::voomWithQualityWeights(dge, design, plot = FALSE)
} else {
  limma::voom(dge, design, plot = FALSE)
}
fit <- eBayes(lmFit(v, design))

deg <- topTable(fit, coef = COEF, number = Inf, sort.by = "P")
deg$gene <- rownames(deg)

# ★ self-gene 제외 후 FDR 재계산 (순환논리 차단)
deg <- deg[deg$gene != GENE_OF_INTEREST, ]
deg$adj.P.Val <- p.adjust(deg$P.Value, "BH")
deg <- deg[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
write.csv(deg, file.path(DIR_RESULTS, "Fig6_microglia_DEG.csv"), row.names = FALSE)
write.csv(deg, file.path(DIR_RESULTS, "Fig6_pseudobulk_DEG_limma.csv"), row.names = FALSE)

message(sprintf("    significant (FDR<0.05): %d  (up %d / down %d)",
                sum(deg$adj.P.Val < 0.05),
                sum(deg$adj.P.Val < 0.05 & deg$logFC > 0),
                sum(deg$adj.P.Val < 0.05 & deg$logFC < 0)))

## Fig6 volcano  (헤더에 약속돼 있었으나 누락돼 있던 것)
p_vol <- ggplot2::ggplot(deg, ggplot2::aes(logFC, -log10(P.Value))) +
  ggplot2::geom_point(ggplot2::aes(color = adj.P.Val < 0.05), size = .8, alpha = .6) +
  ggplot2::scale_color_manual(values = c(`FALSE` = "grey75", `TRUE` = "#B2182B"),
                              guide = "none") +
  ggrepel::geom_text_repel(data = utils::head(deg, 25),
                           ggplot2::aes(label = gene), size = 2.6, max.overlaps = 20) +
  ggplot2::labs(x = sprintf("log2FC per 1 SD of %s (pseudobulk log-CPM)", GENE_OF_INTEREST),
                y = expression(-log[10]~italic(P))) +
  theme_paper()
ggplot2::ggsave(file.path(DIR_FIGURES, "Fig6_microglia_volcano.pdf"),
                p_vol, width = 6, height = 5)

## 관심 유전자 표
FOCUS <- c("NFKB1","NFKBIA","RELA",
           "GLRX","G6PD","TALDO1","TXN","GPX4","CAT",
           "LPL","CD9","ABCA1","TREM2","SPP1","APOE","GPNMB",
           "SMS","ODC1","SRM","AZIN1",
           "SLC3A2","SLC7A1","SLC7A5","SLC7A6","SLC7A11","ASS1","ASL")
foc <- deg[deg$gene %in% FOCUS, ] %>% dplyr::arrange(P.Value)
print(foc, digits = 3)
write.csv(foc, file.path(DIR_RESULTS, "Fig6_microglia_DEG_focus.csv"), row.names = FALSE)

# =============================================================================
# §4. Gene set: CAMERA (primary) + FRY + fgsea
#     ※ inter.gene.cor 는 추정(NA) 과 limma 권장 고정값(0.01) 을 병기
# =============================================================================
E   <- v$E[rownames(v) != GENE_OF_INTEREST, , drop = FALSE]
W   <- v$weights[rownames(v) != GENE_OF_INTEREST, , drop = FALSE]
idx <- ids2indices(HALLMARK, rownames(E))

cm  <- camera(E, idx, design, contrast = COEF, weights = W,
              inter.gene.cor = NA,   sort = FALSE)
cm2 <- camera(E, idx, design, contrast = COEF, weights = W,
              inter.gene.cor = 0.01, sort = FALSE)
fr  <- fry(E, idx, design, contrast = COEF, weights = W, sort = FALSE)

st <- setNames(deg$t, deg$gene)
set.seed(GSEA_SEED)
fg <- fgsea(HALLMARK, st, minSize = GSEA_MIN_SIZE, maxSize = GSEA_MAX_SIZE, nproc = 1)
fg_bp <- fgsea(GO_BP, st, minSize = GSEA_MIN_SIZE, maxSize = GSEA_MAX_SIZE, nproc = 1)

fg_csv <- fg; fg_csv$leadingEdge <- vapply(fg$leadingEdge, paste, character(1), collapse = ";")
fg_bp_csv <- fg_bp; fg_bp_csv$leadingEdge <- vapply(fg_bp$leadingEdge, paste, character(1), collapse = ";")

saveRDS(fg, file.path(DIR_DATA_PROC, "pseudobulk_GSEA_result.rds"))
saveRDS(fg_bp, file.path(DIR_DATA_PROC, "pseudobulk_GSEA_BP_result.rds"))

write_result(fg_csv, "Fig6_pseudobulk_GSEA.csv")
write_result(fg_bp_csv, "Fig6_pseudobulk_GSEA_BP.csv")

pw <- intersect(fg$pathway, rownames(cm))
gs <- data.frame(
  dataset_type      = "scRNA_microglia_donor",
  n                 = ncol(E),
  pathway           = pw,
  setSize           = fg$size[match(pw, fg$pathway)],
  fgsea_NES         = fg$NES[match(pw, fg$pathway)],
  fgsea_padj        = fg$padj[match(pw, fg$pathway)],
  camera_dir        = cm[pw, "Direction"],
  camera_p          = cm[pw, "PValue"],
  camera_FDR        = cm[pw, "FDR"],
  inter_gene_cor    = cm[pw, "Correlation"],
  camera_FDR_fixed  = cm2[pw, "FDR"],          # inter.gene.cor = 0.01
  fry_p             = fr[pw, "PValue"],
  fry_FDR           = fr[pw, "FDR"]
) %>% dplyr::arrange(camera_p)
write.csv(gs, file.path(DIR_RESULTS, "Fig6_microglia_GSEA_camera_vs_fgsea.csv"),
          row.names = FALSE)
print(utils::head(gs, 12), digits = 3)

# =============================================================================
# §5. depth 음성대조 (기존 §5 유지 — 단, 이건 "특이성 증거" 가 아님)
#     depth 는 이미 주모형의 공변량이라 x 효과는 depth 에 직교화돼 있음.
#     부호가 반대라는 건 x 와 depth 가 상관 있다는 뜻일 뿐. 참고용으로만.
# =============================================================================
dsn_d <- model.matrix(~ scale(depth) + dataset, data = donor_meta)
colnames(dsn_d) <- make.names(colnames(dsn_d))
v_d   <- limma::voom(dge, dsn_d, plot = FALSE)
fit_d <- eBayes(lmFit(v_d, dsn_d))
tt_d  <- topTable(fit_d, coef = 2, number = Inf, sort.by = "none")
set.seed(GSEA_SEED)
fg_d  <- fgsea(HALLMARK,
               setNames(tt_d$t, rownames(tt_d))[setdiff(rownames(tt_d), GENE_OF_INTEREST)],
               minSize = GSEA_MIN_SIZE, maxSize = GSEA_MAX_SIZE, nproc = 1)

d2 <- merge(gs[, c("pathway", "fgsea_NES")],
            data.frame(pathway = fg_d$pathway, depth_only_NES = fg_d$NES), by = "pathway")
d2$opposite_sign <- sign(d2$fgsea_NES) != sign(d2$depth_only_NES)
write.csv(d2, file.path(DIR_RESULTS, "Supp_microglia_depth_negative_control.csv"),
          row.names = FALSE)

# =============================================================================
# §5b. 매칭 대조 유전자 음성대조  ★ 이게 진짜 특이성 검정
#   발현량(AveExpr) · 검출률이 GOI 와 가장 비슷한 유전자 N개를 각각 예측변수로
#   똑같이 돌려서, 관측된 pathway 점수가 그 null 분포에서 얼마나 극단적인지 평가.
#   "발현 높은 아무 유전자나 넣어도 같은 pathway 나오는 거 아니냐" 를 정면으로 막음.
# =============================================================================
ave <- rowMeans(lcpm)
det <- rowMeans(pb[rownames(dge), colnames(dge), drop = FALSE] > 0)
d2g <- sqrt(((ave - ave[GENE_OF_INTEREST]) / stats::sd(ave))^2 +
              ((det - det[GENE_OF_INTEREST]) / stats::sd(det))^2)

ctrl_genes <- setdiff(names(sort(d2g))[seq_len(min(N_CTRL + 1, length(d2g)))],
                      GENE_OF_INTEREST)
ctrl_genes <- ctrl_genes[apply(lcpm[ctrl_genes, , drop = FALSE], 1, stats::sd) > 0]
message(sprintf("[5b] 매칭 대조 유전자 %d개 (AveExpr %.2f±?, 검출률 %.2f 기준)",
                length(ctrl_genes), ave[GENE_OF_INTEREST], det[GENE_OF_INTEREST]))

idx_all <- ids2indices(HALLMARK, rownames(v))

set_z <- function(tv, drop_gene = NULL) {         # 집합 수준 competitive z-score
  tv[names(tv) %in% c(GENE_OF_INTEREST, drop_gene)] <- NA_real_
  m <- mean(tv, na.rm = TRUE); s <- stats::sd(tv, na.rm = TRUE)
  vapply(idx_all, function(i) {
    ti <- tv[i]; ti <- ti[!is.na(ti)]
    if (length(ti) < 5) NA_real_ else (mean(ti) - m) / (s / sqrt(length(ti)))
  }, numeric(1))
}

obs_z  <- set_z(fit$t[, COEF])
null_z <- vapply(ctrl_genes, function(g) {
  dsn <- design
  dsn[, COEF] <- as.numeric(scale(lcpm[g, ]))
  set_z(eBayes(lmFit(v, dsn))$t[, COEF], drop_gene = g)   # voom weight 는 재사용(근사)
}, numeric(length(idx_all)))

emp_p <- (rowSums(abs(null_z) >= abs(obs_z), na.rm = TRUE) + 1) / (ncol(null_z) + 1)

nc <- data.frame(
  pathway       = names(obs_z),
  obs_z         = as.numeric(obs_z),
  null_median_z = apply(null_z, 1, stats::median, na.rm = TRUE),
  null_q95_absz = apply(abs(null_z), 1, stats::quantile, .95, na.rm = TRUE),
  emp_p         = as.numeric(emp_p),
  emp_FDR       = p.adjust(emp_p, "BH"),
  n_control     = ncol(null_z)
) %>% dplyr::arrange(emp_p)
write.csv(nc, file.path(DIR_RESULTS, "Supp_matched_control_gene_null.csv"), row.names = FALSE)
print(utils::head(nc, 12), digits = 3)

message(sprintf("  ▶ 매칭 대조 대비 empirical FDR<0.05 pathway: %d / %d",
                sum(nc$emp_FDR < 0.05, na.rm = TRUE), nrow(nc)))

message("\ndone.")
