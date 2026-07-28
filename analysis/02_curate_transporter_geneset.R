# =============================================================================
# 02_curate_transporter_geneset.R      [Fig 2 준비]
# 수송 축 gene set curation (MSigDB 3종 + 수동 y+L/CAT) — 선행연구(대사효소)와 차별화
# 출력: results/transporter_geneset.csv
# =============================================================================
source(here::here("R", "setup.R"))
suppressPackageStartupMessages({ library(msigdbr) })

expr <- load_expr_matrix(require_data(basename(PATH_LGG_NORM)))
tr  <- build_transporter_geneset(universe = expr$Gene)
enz <- build_arg_enzyme_geneset(universe = expr$Gene)

prov <- attr(tr, "provenance")
write_result(prov, "transporter_geneset.csv")
write_result(data.frame(gene = enz), "arg_enzyme_geneset.csv")
cat("\n수송 축:", length(tr), "| 대사효소 축:", length(enz),
    "| 겹침:", length(intersect(tr, enz)), "\n")
stopifnot(GENE_OF_INTEREST %in% tr)
message("완료: gene set curation")
# Methods: "Unlike prior work focusing on arginine metabolic enzymes, we curated a
#  transport-centric gene set, as membrane transporters are systematically excluded
#  from GO 'metabolic process' terms."
