# =============================================================================
# Hallmark GSEA integration plot
# pseudobulk / scRNA within / spatial
#
# Color  : NES
# Size   : -log10(camera_FDR)
# Border : camera FDR significance
# Grey   : camera FDR >= 0.05
#
# Designed for publication figure
# =============================================================================
source(here::here("R", "setup.R"))
library(tidyverse)
library(scales)

# -----------------------------------------------------------------------------
# 0. Input files
# -----------------------------------------------------------------------------
FILE_TCGA <- file.path(DIR_RESULTS, "Fig4_GSEA_camera_vs_fgsea.csv")
FILE_PB <- file.path(DIR_RESULTS, "Fig6_microglia_GSEA_camera_vs_fgsea.csv")
FILE_WITHIN <- file.path(DIR_RESULTS, "Fig6_wb_gsea_within.csv")
FILE_SPATIAL <- file.path(DIR_RESULTS, "Fig8_spatial_gsea_HALL.csv")

# output
OUT_DIR <- DIR_RESULTS
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)


# -----------------------------------------------------------------------------
# 1. Read data
# -----------------------------------------------------------------------------
tcga <- read.csv(FILE_TCGA, check.names = FALSE)
pb <- read.csv(FILE_PB, check.names = FALSE)
within <- read.csv(FILE_WITHIN, check.names = FALSE)
spatial <- read.csv(FILE_SPATIAL, check.names = FALSE)


# -----------------------------------------------------------------------------
# 2. Standardize columns
# -----------------------------------------------------------------------------

# bulk TCGA_LGG
tcga2 <- tcga %>%
  dplyr::filter(cohort == "TCGA_LGG") %>%
  transmute(
    pathway,
    NES = fgsea_NES,
    camera_p,
    camera_FDR,
    Direction = camera_dir,
    analysis = "Bulk"
  )


# pseudobulk
pb2 <- pb %>%
  transmute(
    pathway,
    NES = fgsea_NES,
    camera_p,
    camera_FDR,
    Direction = camera_dir,
    analysis = "Pseudobulk"
  )

# scRNA within
within2 <- within %>%
  transmute(
    pathway,
    NES,
    camera_p,
    camera_FDR,
    Direction,
    analysis = "scRNA within"
  )

# spatial
spatial2 <- spatial %>%
  transmute(
    pathway,
    NES,
    camera_p,
    camera_FDR,
    Direction,
    analysis = "Spatial"
  )


# -----------------------------------------------------------------------------
# 3. Combine
# -----------------------------------------------------------------------------

gsea <- bind_rows(
  tcga2,
  pb2,
  within2,
  spatial2
) %>%
  mutate(
    analysis = factor(
      analysis,
      levels = c(
        "Bulk",
        "Pseudobulk",
        "scRNA within",
        "Spatial"
      )
    )
  )


# -----------------------------------------------------------------------------
# 4. Clean pathway names
# -----------------------------------------------------------------------------

gsea <- gsea %>%
  mutate(
    pathway_clean = pathway %>%
      str_remove("^HALLMARK_") %>%
      str_replace_all("_", " ") %>%
      str_to_title()
  )


# -----------------------------------------------------------------------------
# 5. Camera-FDR significance
# -----------------------------------------------------------------------------

gsea <- gsea %>%
  mutate(
    camera_sig = camera_FDR < 0.05,
    
    # cap extreme values so one very significant pathway
    # doesn't dominate the size scale
    neglog10_FDR = -log10(pmax(camera_FDR, 1e-300)),
    neglog10_FDR_plot = pmin(neglog10_FDR, 15)
  )


# -----------------------------------------------------------------------------
# 6. Select pathways
#
# Use UNION of camera-FDR significant pathways.
# This is much better than selecting top pathways independently in each panel,
# because the same pathway can then be compared across all three analyses.
# -----------------------------------------------------------------------------

sig_pathways <- gsea %>%
  filter(camera_FDR < 0.05) %>%
  distinct(pathway) %>%
  pull(pathway)


# -----------------------------------------------------------------------------
# 7. Rank pathways
#
# Primary:
#   strongest overall NES
#
# Secondary:
#   number of analyses in which pathway is camera-FDR significant
#
# This places pathways with reproducible signals near the top.
# -----------------------------------------------------------------------------

pathway_order <- gsea %>%
  filter(pathway %in% sig_pathways) %>%
  group_by(pathway, pathway_clean) %>%
  summarise(
    n_sig = sum(camera_FDR < 0.05, na.rm = TRUE),
    mean_abs_NES = mean(abs(NES), na.rm = TRUE),
    max_abs_NES = max(abs(NES), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(
    desc(n_sig),
    desc(mean_abs_NES),
    desc(max_abs_NES)
  ) %>%
  pull(pathway)


# -----------------------------------------------------------------------------
# 8. Add non-significant rows for all selected pathways
# -----------------------------------------------------------------------------

plot_df <- gsea %>%
  filter(pathway %in% sig_pathways) %>%
  mutate(
    pathway = factor(
      pathway,
      levels = rev(pathway_order)
    )
  )


# -----------------------------------------------------------------------------
# 9. Plot
# -----------------------------------------------------------------------------

p <- ggplot(
  plot_df,
  aes(
    x = analysis,
    y = pathway
  )
) +
  
  # non-significant background
  geom_point(
    aes(
      size = neglog10_FDR_plot
    ),
    shape = 21,
    fill = "grey92",
    colour = "grey75",
    stroke = 0.25
  ) +
  
  # significant points
  geom_point(
    data = plot_df %>%
      filter(camera_FDR < 0.05),
    aes(
      size = neglog10_FDR_plot,
      fill = NES
    ),
    shape = 21,
    colour = "black",
    stroke = 0.25,
    alpha = 0.95
  ) +
  
  # NES color
  scale_fill_gradient2(
    low = "#3B6FB6",
    mid = "white",
    high = "#C43D3D",
    midpoint = 0,
    name = "NES",
    limits = c(
      floor(min(plot_df$NES, na.rm = TRUE)),
      ceiling(max(plot_df$NES, na.rm = TRUE))
    )
  ) +
  
  # FDR -> point size
  scale_size_continuous(
    name = expression(-log[10]("camera FDR")),
    range = c(2, 12),
    breaks = c(1, 2, 3, 5, 10, 15),
    limits = c(0, 15)
  ) +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  theme_minimal(
    base_size = 11
  ) +
  
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    
    panel.grid.major.y = element_line(
      colour = "grey92",
      linewidth = 0.35
    ),
    
    axis.text.x = element_text(
      size = 11,
      face = "bold",
      colour = "black"
    ),
    
    axis.text.y = element_text(
      size = 9,
      colour = "black"
    ),
    
    axis.ticks = element_blank(),
    
    legend.position = "right",
    
    legend.title = element_text(
      size = 9,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 8
    ),
    
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0
    ),
    
    plot.subtitle = element_text(
      size = 10,
      colour = "grey35",
      hjust = 0
    ),
    
    plot.margin = margin(
      10, 15, 10, 10
    )
  )


# -----------------------------------------------------------------------------
# 10. Add title/subtitle
# -----------------------------------------------------------------------------
# 
p <- p +
  labs(
    title = "Hallmark pathway activity across analytical levels",
    subtitle = "Color indicates normalized enrichment score (NES); point size indicates camera FDR"
  )


# -----------------------------------------------------------------------------
# 11. Print
# -----------------------------------------------------------------------------

print(p)


# -----------------------------------------------------------------------------
# 12. Save
# -----------------------------------------------------------------------------

ggsave(
  file.path(
    OUT_DIR,
    "Hallmark_GSEA_pseudobulk_within_spatial.pdf"
  ),
  p,
  width = 9,
  height = max(6.5, 0.28 * length(pathway_order) + 2),
  units = "in",
  device = cairo_pdf
)


p_aixs <- p + theme_blank_frame(lp = "right", keep_grid_x = TRUE)+
  theme( panel.grid.major = element_line(colour = "grey92", linewidth = 0.3), panel.grid.minor = element_line(colour = "grey92", linewidth = 0.2) )

save_emf(p, "Fig8_Hallmark_GSEA_pseudobulk_within_spatial.emf", height = 9.5, width = 3.5, legend_p = "right", draw_minor_grid = TRUE, draw_major_grid = TRUE)




# -----------------------------------------------------------------------------
# 13. Save integrated table
# -----------------------------------------------------------------------------

write.csv(
  gsea %>%
    filter(pathway %in% sig_pathways) %>%
    arrange(pathway, analysis),
  file.path(
    OUT_DIR,
    "Hallmark_GSEA_integrated_cameraFDR.csv"
  ),
  row.names = FALSE
)


message(
  "\nSaved to: ",
  OUT_DIR,
  "\nSelected pathways: ",
  length(pathway_order)
)