


# feature끼리 그리기 -----------------------------------------------------------



features <- c(
  "IFN_ALPHA",
  "NFKB",
  "TRANSPORT",
  "ARG_ENZYME",
  "PPP_REACTOME",
  "GLYCOLYSIS"
)

## 모든 pathway에서 공통 color scale 계산
lims <- range(as.matrix(colData(cds)[, features]), na.rm = TRUE)

plist <- lapply(features, function(f){
  
  plot_cells(
    cds,
    color_cells_by = f,
    cell_size = 0.6,
    label_cell_groups = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    label_roots = FALSE
  ) +
    ggtitle(f) +
    scale_color_viridis_c(
      option = "C",
      limits = lims,
      oob = scales::squish
    ) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5)
    )
  
})

p <- wrap_plots(plist, ncol = 3, guides = "collect")



p



# 유전자와 함께 그리기 -------------------------------------------------------------

library(patchwork)

p1 <- plot_cells(
  cds,
  genes = "SLC7A7",
  cell_size = 0.6,
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  label_roots = FALSE
) +
  ggtitle("SLC7A7")

p2 <- plot_cells(
  cds,
  color_cells_by = "IFN_ALPHA",
  cell_size = 0.6,
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  label_roots = FALSE
) +
  ggtitle("IFN_ALPHA")

p3 <- plot_cells(
  cds,
  color_cells_by = "NFKB",
  cell_size = 0.6,
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  label_roots = FALSE
) +
  ggtitle("NFKB")

p4 <- plot_cells(
  cds,
  color_cells_by = "HYPOXIA",
  cell_size = 0.6,
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  label_roots = FALSE
) +
  ggtitle("HYPOXIA")

(p1 | p2) / (p3 | p4)



# pseudotime 여러개 붙여 그리기 ---------------------------------------------------


df <- data.frame(
  pseudotime = pseudotime(cds),
  SLC7A7 = as.numeric(exprs(cds)["SLC7A7", ]),
  IFN_ALPHA = colData(cds)$IFN_ALPHA
)

df_long <- pivot_longer(
  df,
  cols = c(SLC7A7, IFN_ALPHA),
  names_to = "Feature",
  values_to = "Value"
)

ggplot(df_long, aes(pseudotime, Value)) +
  geom_point(size = 0.3, alpha = 0.2) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  facet_wrap(~Feature, scales = "free_y") +
  theme_classic()







df <- data.frame(
  pseudotime = pseudotime(cds),
  score = colData(cds)$NFKB
)

ggplot(df, aes(pseudotime, score)) +
  geom_point(size = 0.3, alpha = 0.3) +
  geom_smooth(method = "loess", se = TRUE, color = "red") +
  theme_classic() +
  labs(
    x = "Pseudotime",
    y = "IFN_ALPHA UCell score"
  )
