# =============================================================================
# R/gene_sets.R  ---  수송체 gene set curation / transport-centric gene set
# -----------------------------------------------------------------------------
# 선행연구(ASL)는 아르기닌 "대사효소"를 썼음. 본 연구는 "수송 축"으로 재정의한다.
# 근거: 막수송체는 GO 'metabolic process' 계열에서 체계적으로 제외되므로,
#       대사 gene set 만 쓰면 y+L/CAT 계열이 구조적으로 누락된다.
# =============================================================================

build_transporter_geneset <- function(manual = TRANSPORTER_MANUAL,
                                      msig_sets = MSIGDB_TRANSPORT_SETS,
                                      universe = NULL) {
  msig <- dplyr::bind_rows(
    msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP"),
    msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
  )
  from_msig <- msig %>% dplyr::filter(gs_name %in% msig_sets)
  genes <- sort(unique(c(from_msig$gene_symbol, manual)))
  if (!is.null(universe)) genes <- intersect(genes, universe)

  prov <- data.frame(gene = genes,
                     manual = genes %in% manual,
                     msigdb = genes %in% from_msig$gene_symbol, stringsAsFactors = FALSE)
  attr(genes, "provenance") <- prov
  message("transporter gene set: ", length(genes), " genes (manual ",
          sum(prov$manual), " / MSigDB ", sum(prov$msigdb), ")")
  genes
}

build_arg_enzyme_geneset <- function(manual = ARG_ENZYME_SET, universe = NULL) {
  g <- sort(unique(manual)); if (!is.null(universe)) g <- intersect(g, universe); g
}

get_msig_genes <- function(gs_name, category = "H", subcategory = NULL) {
  m <- msigdbr::msigdbr(species = "Homo sapiens", category = category, subcategory = subcategory)
  unique(m$gene_symbol[m$gs_name == gs_name])
}
