cnt <- SeuratObject::LayerData(mg, assay = "RNA", layer = "counts")
cnt
v <- cnt@matrix[1:200]
print(v[1,1])
all(v == round(v))                      # FALSE 면 확정
summary(v); range(Matrix::colSums(cnt))  # colSums 이 2배 안이면 확정
SeuratObject::Layers(mye)                # counts 말고 다른 레이어 있는지


cnt["SLC7A7", 1:20]
as.matrix(cnt["SLC7A7", 1:20])
Seurat::Assays(microglia)
colnames(microglia[[]])
list(
  dimensions = dim(microglia),
  assays = Assays(microglia),
  default_assay = DefaultAssay(microglia),
  metadata = colnames(microglia[[]]),
  reductions = Reductions(microglia)
)
