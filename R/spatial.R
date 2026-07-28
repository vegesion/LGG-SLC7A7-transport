# =============================================================================
# R/spatial.R  ---  공간전사체 / spatial transcriptomics (Fig 8, 유일한 차별점)
# SLC7A7 공간 분포 · niche 주석 · 저산소/괴사 중첩 · spot 수준 SLC7A7 vs ASS1
# =============================================================================
load_spatial_sample <- function(dir, slice = basename(dir)) {
  obj <- Seurat::Load10X_Spatial(data.dir = dir, slice = slice)
  obj <- Seurat::SCTransform(obj, assay = "Spatial", verbose = FALSE)
  obj$sample <- slice; obj
}

load_spatial_all <- function(root = DIR_SPATIAL) {
  dirs <- list.dirs(root, recursive = FALSE)
  if (!length(dirs)) stop("공간 데이터 폴더가 없습니다: ", root, call. = FALSE)
  stats::setNames(lapply(dirs, load_spatial_sample), basename(dirs))
}

add_spatial_modules <- function(obj, sets = NULL) {
  if (is.null(sets)) sets <- list(
    hypoxia   = get_msig_genes("HALLMARK_HYPOXIA"),
    myeloid   = c("AIF1","CD68","CSF1R","ITGAM","TYROBP","C1QA","C1QB","P2RY12","TMEM119"),
    necrosis  = c("VEGFA","CA9","ADM","NDRG1","SLC2A1","LDHA","PGK1"),
    transport = build_transporter_geneset(universe = rownames(obj)))
  for (nm in names(sets)) {
    g <- intersect(sets[[nm]], rownames(obj)); if (length(g) < 3) next
    obj <- Seurat::AddModuleScore(obj, features = list(g), name = paste0(nm, "_"), seed = 42)
    names(obj@meta.data)[names(obj@meta.data) == paste0(nm, "_1")] <- paste0(nm, "_score")
  }
  obj
}

spatial_spot_correlation <- function(obj, gene = GENE_OF_INTEREST,
                                     partners = c("ASS1", "ASL", "SLC3A2"),
                                     scores = c("hypoxia_score","myeloid_score","necrosis_score")) {
  vars <- c(gene, intersect(partners, rownames(obj)))
  d <- Seurat::FetchData(obj, vars = vars)
  for (s in intersect(scores, names(obj@meta.data))) d[[s]] <- obj[[s]][, 1]
  x <- d[[gene]]
  dplyr::bind_rows(lapply(setdiff(names(d), gene), function(v) {
    ct <- suppressWarnings(stats::cor.test(x, d[[v]], method = "spearman"))
    data.frame(sample = obj$sample[1], variable = v, rho = unname(ct$estimate), p = ct$p.value)
  })) %>% dplyr::mutate(FDR = stats::p.adjust(p, "BH"))
}

plot_spatial_gene <- function(obj, gene = GENE_OF_INTEREST, blank = FALSE) {
  p <- Seurat::SpatialFeaturePlot(obj, features = gene) +
    ggplot2::labs(title = paste(obj$sample[1], gene))
  if (blank) p <- p + theme_blank_frame(); p
}

spatial_niche_test <- function(obj, gene = GENE_OF_INTEREST, score = "hypoxia_score") {
  v <- Seurat::FetchData(obj, vars = gene)[, 1]; s <- obj[[score]][, 1]
  grp <- factor(ifelse(s >= stats::quantile(s, 0.75, na.rm = TRUE), "niche-high", "rest"),
                levels = c("rest", "niche-high"))
  list(test = stats::wilcox.test(v ~ grp),
       medians = tapply(v, grp, stats::median, na.rm = TRUE),
       data = data.frame(expr = v, group = grp, sample = obj$sample[1]))
}
