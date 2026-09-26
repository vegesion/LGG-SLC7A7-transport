
library(harmony)
microglia <- RunHarmony(
  object = microglia,
  group.by.vars = "donor_id"
)

microglia <- RunUMAP(
  microglia,
  reduction = "harmony",
  dims = 1:30,
  reduction.name = "umap.harmony"
)
microglia@reductions
names(microglia@meta.data)

cnt <- GetAssayData(
  microglia,
  assay="RNA",
  layer="counts"
)

cds <- new_cell_data_set(
  cnt,
  cell_metadata = microglia@meta.data,
  gene_metadata =
    data.frame(
      gene_short_name = rownames(cnt),
      row.names = rownames(cnt)
    )
)

reducedDims(cds)$HARMONY <-
  Embeddings(
    microglia,
    "harmony"
  )

reducedDims(cds)


cds <- cluster_cells(
  cds,
  reduction_method = "HARMONY"
)

cds <- learn_graph(
  cds,
  use_partition = FALSE,
  close_loop = FALSE
)




args(RunHarmony)
getAnywhere(RunHarmony.Seurat)
args(RunHarmony)
harmony::RunHarmony
packageVersion("harmony")
packageVersion("Seurat")
find("RunHarmony")
