# =============================================================================
# 08_sc_myeloid_donor.R                  [Fig 6]
# 아종 UMAP · donor 수준 발현 vs 아종 비율 회귀(핵심) · MG vs MDM · pseudobulk DEG/GSEA
# =============================================================================
source(here::here("R", "setup.R"))
suppressPackageStartupMessages({
  library(Seurat); library(edgeR); library(limma); library(fgsea); library(ggpubr) })

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

## donor pseudobulk DEG (paired limma) + Hallmark GSEA
mg <- label_slc_status(mg)
pb <- run_pseudobulk_limma(mg); write_result(pb$result, "Fig6_pseudobulk_DEG.csv")
saveRDS(pb, file.path(DIR_DATA_PROC, "pb_fit.rds"))
gs <- run_fgsea_pseudobulk(pb$fit, coef = pb$coef); write_result(gs$table, "Fig6_pseudobulk_GSEA.csv")
print(utils::head(gs$table, 15))
ggplot2::ggsave(file.path(DIR_FIGURES, "Fig6D_pseudobulk_GSEA_NES.pdf"),
                plot_nes_bar(gs$table, top_n = 12), width = 7, height = 5)
message("완료: Fig6")


# pseudobulk_연속형donor수준 ---------------------------------------------------
library(edgeR)

mg <- readRDS(file.path(DIR_DATA_PROC, "microglia.rds"))
pb <- donor_pseudobulk(mg, layer = "counts")          # genes × donors (전체 세포 합산)

# donor 수준 공변량
md <- data.frame(donor = mg$donor_id, 
                 g = FetchData(mg, vars = GENE_OF_INTEREST)[,1],
                 depth = mg[[DEPTH_COVARIATE]][,1],
                 ds = as.character(mg[[find_batch_col(mg)]][,1]),
                 mt = PercentageFeatureSet(mg, pattern = "^MT-")) %>%
  group_by(donor) %>% summarise(across(c(g, depth, mt), mean), ds = first(ds))
md$g <- as.numeric(scale(md$g))
md <- md[match(colnames(pb), md$donor), ]

dge <- DGEList(pb); keep <- filterByExpr(dge); dge <- calcNormFactors(dge[keep, ])
design <- model.matrix(~ g + depth + mt + ds, data = md)   # ★ depth·MT·batch 모두 보정

v   <- voom(dge, design); fit <- eBayes(lmFit(v, design))
res <- topTable(fit, coef = "g", number = Inf)
res <- res[rownames(res) != GENE_OF_INTEREST, ]            # self 제외
head(res, 20); sum(res$adj.P.Val < 0.05)

# gs <- run_fgsea_pseudobulk(res$fit, coef = pb$coef); write_result(gs$table, "Fig6_pseudobulk_GSEA.csv")

# 이상치/레버리지 점검: 상위 유전자 몇 개를 donor 산점도로
for (g in c("CPEB4","COTL1","QKI")) {
  y <- log2(edgeR::cpm(dge)[g, ] + 1)
  plot(md$g, y, main = g, xlab = "SLC7A7 (donor z)", ylab = "log2 CPM")
  abline(lm(y ~ md$g), col = "red")
}


ranks <- setNames(res$t, rownames(res)); ranks <- sort(ranks, decreasing = TRUE)
set.seed(GSEA_SEED)
gs <- fgsea::fgsea(get_hallmark_list(), ranks, minSize = 10, maxSize = 500, eps = 0)
write_result(gs$table, "Fig6_pseudobulk_GSEA_modified.csv")
ggplot2::ggsave(file.path(DIR_FIGURES, "Fig6D_pseudobulk_GSEA_NES_modified.pdf"),
                plot_nes_bar(gs, top_n = 10), width = 7, height = 5)
head(gs[order(gs$pval), c("pathway","NES","padj")], 15)


#딱 하나만 확인하고 넘어가세요. OXPHOS가 진짜인지 미토콘드리아 유전자 잔재인지 봐야 해요:
le <- gs$leadingEdge[[which(gs$pathway == "HALLMARK_OXIDATIVE_PHOSPHORYLATION")]]
head(le, 30)


#하나 다시 해볼 만한 게 생겼어요. FAO + OXPHOS 항진은 lipid-associated macrophage(LAM/TREM2⁺) 표현형의 교과서적 특징이에요. 예전에 efferocytosis 모듈(TREM2·APOE·LPL·ABCA1)을 봤을 때 null이었는데, 그건 AddModuleScore + 세포 단위로 본 거였잖아요. 지금 파이프라인으로 다시 보면 결과가 달라질 수 있어요:
lam <- intersect(c("TREM2","APOE","LPL","ABCA1","GPNMB","CD9","SPP1","FABP5",
                   "LIPA","CTSB","CTSD","NR1H3","PLIN2"), rownames(dge))
# 위 연속형 limma fit 에서 바로 검정 (같은 design, self 제외 불필요)
limma::topTable(fit, coef = "g", number = Inf)[lam, c("logFC","adj.P.Val")]
# 또는 gene set 수준으로
fgsea::fgsea(list(LAM = lam), ranks, minSize = 5)[, c("pathway","NES","padj")]



#지질을 들여오고(LPL↑) 내보내지 않으며(ABCA1↓) 태우는(FAO·OXPHOS↑) 방향이에요. 세 갈래가 대사적으로 일관돼요. 즉 지질을 축적·수출하는 LAM이 아니라, 유입해서 산화 연료로 쓰는 상태입니다.
#그래서 서술은 이렇게 가는 게 정확해요: "an oxidative, lipid-utilizing myeloid state that is distinct from the canonical TREM2⁺ lipid-associated and SPP1⁺ TAM programs." 알려진 두 축 어디에도 안 들어가는 새로운 축이라는 건 논문에 유리한 주장이에요 — 다만 "새롭다"고 쓰려면 그걸 보여줘야 하니, 직교성 패널을 하나 만드세요:

axes <- list(
  LAM_TREM2  = c("TREM2","APOE","GPNMB","CD9","LPL","SPP1","FABP5","CTSD"),
  SPP1_TAM   = c("SPP1","MARCO","VCAN","FN1","CD44"),
  Homeostatic= c("P2RY12","TMEM119","CX3CR1","SALL1","SELPLG"),
  IFN_TAM    = c("ISG15","IFIT1","IFIT3","MX1","STAT1","IRF7"),
  Inflam_TAM = c("IL1B","TNF","CXCL8","NFKBIA","CCL3"),
  OxPhos_FAO = c(le[1:30], "HADHB","CPT1A","ACADM")   # 이번에 나온 축
)
fgsea::fgsea(lapply(axes, intersect, names(ranks)), ranks, minSize = 4)[, c("pathway","NES","padj")]

