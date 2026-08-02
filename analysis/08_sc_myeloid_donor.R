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

suppressPackageStartupMessages({
  library(Matrix); library(edgeR); library(limma); library(fgsea); library(msigdbr)
  library(lme4); library(dplyr); library(ggplot2); library(ggrepel)
})

# ── 파라미터 ────────────────────────────────────────────────────────────────
MG_LABELS <- c("TAM-MG aging sig", "TAM-MG pro-infl I",
               "TAM-MG pro-infl II", "TAM-MG prolif")
LABEL_COL        <- "annotation_level_4"
DEPTH_COL        <- "nCount_RNA"      # ★ nFeature 아님 (UMI 총량)
MIN_DONOR_CELLS  <- 50                # microglia 는 donor당 세포수가 적음
MIN_DATASET_DONOR<- 2

HALLMARK <- local({
  m <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"),
                error = function(e) msigdbr(species = "Homo sapiens", category = "H"))
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
# §2. Methods 용 숫자 — ICC(dataset) 와 depth 상관  ★ 각 1문장씩만 보고
# =============================================================================
x_raw <- as.numeric(dm[GENE_OF_INTEREST, ])
icc_df <- data.frame(donor = colnames(dm), expr = x_raw,
                     depth = donor_meta$depth, dataset = donor_meta$dataset)
fit_icc <- lmer(expr ~ depth + (1 | dataset), data = icc_df, REML = TRUE)
vc  <- as.data.frame(VarCorr(fit_icc))
ICC <- vc$vcov[vc$grp == "dataset"] / sum(vc$vcov)

icc_out <- data.frame(
  n_donor      = nrow(icc_df),
  n_dataset    = nlevels(icc_df$dataset),
  ICC_dataset  = ICC,
  rho_expr_depth = cor(icc_df$expr, icc_df$depth, method = "spearman"),
  rho_expr_ncell = cor(icc_df$expr, donor_meta$n_cell, method = "spearman"),
  depth_var_by_dataset_R2 = summary(lm(depth ~ dataset, icc_df))$r.squared
)
write.csv(icc_out, file.path(DIR_RESULTS, "Supp_microglia_ICC_depth.csv"), row.names = FALSE)
print(icc_out, digits = 3)
message("  ▶ Methods 문장: \"Dataset accounted for ", round(100*ICC), "% of donor-level variance ",
        "in SLC7A7 (ICC = ", round(ICC, 2), "); dataset and nCount_RNA were included as covariates.\"")

# =============================================================================
# §3. limma 유전자 수준 DEG  ★ Fig6 주 결과
# =============================================================================
x      <- scale(x_raw)[, 1]
design <- model.matrix(~ x + donor_meta$depth + donor_meta$dataset)
colnames(design) <- make.names(colnames(design))
COEF   <- "x"

dge <- DGEList(pb)
keep <- filterByExpr(dge, design = design, min.count = PB_FILTER_MIN_COUNT)
keep[GENE_OF_INTEREST] <- TRUE                      # 예측변수는 남겨 두고 랭킹에서만 제거
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
message(sprintf("[3] genes after filterByExpr: %d", nrow(dge)))

v   <- voom(dge, design, plot = FALSE)
fit <- eBayes(lmFit(v, design))
deg <- topTable(fit, coef = COEF, number = Inf, sort.by = "P")
deg$gene <- rownames(deg)

# ★ self-gene 제외 후 FDR 재계산 (순환논리 차단)
deg <- deg[deg$gene != GENE_OF_INTEREST, ]
deg$adj.P.Val <- p.adjust(deg$P.Value, "BH")
deg <- deg[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
write.csv(deg, file.path(DIR_RESULTS, "Fig6_pseudobulk_DEG_limma.csv"), row.names = FALSE)

message(sprintf("    significant (FDR<0.05): %d  (up %d / down %d)",
                sum(deg$adj.P.Val < 0.05),
                sum(deg$adj.P.Val < 0.05 & deg$logFC > 0),
                sum(deg$adj.P.Val < 0.05 & deg$logFC < 0)))

## 관심 유전자 표 (논문 본문용)
FOCUS <- c("NFKB1","NFKBIA","RELA",                       # NF-kB 축 (Rotoli 2018)
           "GLRX","G6PD","TALDO1","TXN","GPX4","CAT",     # PPP / NADPH / redox
           "LPL","CD9","ABCA1","TREM2","SPP1","APOE","GPNMB",  # 지질 · LAM 감별
           "SMS","ODC1","SRM","AZIN1",                    # 폴리아민 (아르기닌 축)
           "SLC3A2","SLC7A1","SLC7A5","SLC7A6","SLC7A11","ASS1","ASL")
foc <- deg[deg$gene %in% FOCUS, ] %>% arrange(P.Value)
print(foc, digits = 3)
write.csv(foc, file.path(DIR_RESULTS, "Fig6_microglia_DEG_focus.csv"), row.names = FALSE)


# =============================================================================
# §4. Gene set: CAMERA (primary) + FRY + fgsea — bulk 표와 동일 형식
# =============================================================================
E   <- v$E[rownames(v) != GENE_OF_INTEREST, , drop = FALSE]
W   <- v$weights[rownames(v) != GENE_OF_INTEREST, , drop = FALSE]
idx <- ids2indices(HALLMARK, rownames(E))

cm <- camera(E, idx, design, contrast = COEF, weights = W, inter.gene.cor = NA, sort = FALSE)
fr <- fry(E, idx, design, contrast = COEF, weights = W, sort = FALSE)

st <- setNames(deg$t, deg$gene)
set.seed(GSEA_SEED)
fg <- fgsea(HALLMARK, st, minSize = GSEA_MIN_SIZE, maxSize = GSEA_MAX_SIZE, nproc = 1)

fg_csv <- fg; fg_csv$leadingEdge <- as.character(fg_csv$leadingEdge)
saveRDS(fg, file.path(DIR_DATA_PROC, "pseudobulk_GSEA_result.rds"))
write_result(fg_csv, "Fig6_pseudobulk_GSEA.csv")

pw <- intersect(fg$pathway, rownames(cm))
gs <- data.frame(
  dataset_type   = "scRNA_microglia_donor",
  n              = ncol(E),
  pathway        = pw,
  setSize        = fg$size[match(pw, fg$pathway)],
  fgsea_NES      = fg$NES[match(pw, fg$pathway)],
  fgsea_padj     = fg$padj[match(pw, fg$pathway)],
  camera_dir     = cm[pw, "Direction"],
  camera_p       = cm[pw, "PValue"],
  camera_FDR     = cm[pw, "FDR"],
  inter_gene_cor = cm[pw, "Correlation"],
  fry_p          = fr[pw, "PValue"],
  fry_FDR        = fr[pw, "FDR"]
) %>% arrange(camera_p)
write.csv(gs, file.path(DIR_RESULTS, "Fig6_microglia_GSEA_camera_vs_fgsea.csv"), row.names = FALSE)
print(head(gs, 12), digits = 3)

message(sprintf("\n  fgsea padj<0.05 : %d / %d", sum(gs$fgsea_padj < 0.05, na.rm = TRUE), nrow(gs)))
message(sprintf("  CAMERA FDR<0.05 : %d", sum(gs$camera_FDR < 0.05, na.rm = TRUE)))
message(sprintf("  inter-gene cor  : median %.3f (max %.3f)",
                median(gs$inter_gene_cor, na.rm = TRUE), max(gs$inter_gene_cor, na.rm = TRUE)))

# =============================================================================
# §5. depth 음성대조 — Results 에 한 문장 (D2)
# =============================================================================
dsn_d <- model.matrix(~ scale(donor_meta$depth) + donor_meta$dataset)
colnames(dsn_d) <- make.names(colnames(dsn_d))
v_d   <- voom(dge, dsn_d, plot = FALSE)
fit_d <- eBayes(lmFit(v_d, dsn_d))
tt_d  <- topTable(fit_d, coef = 2, number = Inf, sort.by = "none")
set.seed(GSEA_SEED)
fg_d  <- fgsea(HALLMARK, setNames(tt_d$t, rownames(tt_d))[setdiff(rownames(tt_d), GENE_OF_INTEREST)],
               minSize = GSEA_MIN_SIZE, maxSize = GSEA_MAX_SIZE, nproc = 1)

d2 <- merge(gs[, c("pathway", "fgsea_NES")],
            data.frame(pathway = fg_d$pathway, depth_only_NES = fg_d$NES), by = "pathway")
d2$opposite_sign <- sign(d2$fgsea_NES) != sign(d2$depth_only_NES)
write.csv(d2, file.path(DIR_RESULTS, "Supp_microglia_depth_negative_control.csv"), row.names = FALSE)
print(head(d2[order(-abs(d2$fgsea_NES)), ], 8), digits = 3)
message(sprintf("\n  ▶ Results 문장: depth 단독 예측변수 사용 시 상위 pathway 중 %d/%d 가 반대 부호",
                sum(head(d2[order(-abs(d2$fgsea_NES)), ], 8)$opposite_sign), 8))

message("\ndone.")
message(" Fig6 주 그림 : figures/Fig6_microglia_volcano.pdf")
message(" 논문 Table   : results/Fig6_microglia_GSEA_camera_vs_fgsea.csv (bulk 표와 같은 열 구성)")

