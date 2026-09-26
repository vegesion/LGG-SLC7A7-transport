# ============================================================
# Fig 6A (final) — Hallmark GSEA, NES rank order + camera robustness
#   input : Fig6_microglia_GSEA_camera_vs_fgsea.csv
#   fgsea = 효과 크기 / camera = robustness. fry 미사용.
# ============================================================
source(here::here("R", "setup.R"))
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(stringr); library(patchwork)
})

FDR    <- 0.05
COL_UP <- "#B2182B"
COL_DN <- "#2166AC"

tidy_name <- function(x) {
  x <- str_remove(x, "^HALLMARK_")
  str_to_sentence(str_replace_all(x, "_", " "))
}

theme_fig <- function(base = 9) {
  theme_classic(base_size = base) +
    theme(axis.text       = element_text(colour = "black"),
          axis.line       = element_line(linewidth = 0.3),
          axis.ticks      = element_line(linewidth = 0.3),
          plot.title      = element_text(face = "bold", size = base + 1),
          plot.subtitle   = element_text(size = base - 1, colour = "grey35"),
          legend.key.size = unit(3, "mm"),
          legend.title    = element_text(size = base - 1),
          legend.text     = element_text(size = base - 2))
}

# ---- 데이터 ------------------------------------------------
d <- fread(file.path(DIR_RESULTS, "Fig6_microglia_GSEA_camera_vs_fgsea.csv"))

d <- d[fgsea_padj < FDR][order(fgsea_NES)]          # 26 sets, NES 오름차순
d[, dir    := factor(fifelse(fgsea_NES > 0, "Up", "Down"), c("Up", "Down"))]
d[, robust := camera_FDR < FDR]
d[, label  := factor(tidy_name(pathway), levels = tidy_name(pathway))]
d[, cam_nlp := -log10(pmax(camera_FDR, .Machine$double.xmin))]

# ---- (좌) NES lollipop --------------------------------------
p1 <- ggplot(d, aes(x = fgsea_NES, y = label)) +
  geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = fgsea_NES, yend = label, colour = dir),
               linewidth = 0.5) +
  geom_point(aes(size = setSize, colour = dir,
                 fill = fifelse(robust, as.character(dir), NA_character_)),
             shape = 21, stroke = 0.7) +
  scale_colour_manual(values = c(Up = COL_UP, Down = COL_DN), guide = "none") +
  scale_fill_manual(values = c(Up = COL_UP, Down = COL_DN),
                    na.value = "white", guide = "none") +
  scale_size_continuous(range = c(1.6, 4.5), name = "Set size",
                        breaks = c(50, 100, 200)) +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0.10))) +
  labs(x = "NES (fgsea)", y = NULL,
       title    = "Hallmark pathway enrichment",
       subtitle = sprintf("FDR < %.2f (n = %d donors) · filled = camera FDR < %.2f",
                          FDR, unique(d$n), FDR)) +
  theme_fig()

# ---- (우) camera FDR 막대 -----------------------------------
p2 <- ggplot(d, aes(x = cam_nlp, y = label)) +
  geom_col(aes(fill = robust), width = 0.62) +
  geom_vline(xintercept = -log10(FDR), linetype = "22",
             colour = "grey35", linewidth = 0.35) +
  scale_fill_manual(values = c(`TRUE` = "grey25", `FALSE` = "grey80"),
                    guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = expression(-log[10]~"camera FDR"), y = NULL) +
  theme_fig() +
  theme(axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y  = element_blank())

pA <- p1 + p2 + plot_layout(widths = c(1, 0.34))


dir.create("figures", showWarnings = FALSE)
ggsave("figures/Fig6A_hallmark.pdf", pA, width = 7.2, height = 5.4,
       units = "in", device = cairo_pdf)
ggsave("figures/Fig6A_hallmark.png", pA, width = 7.2, height = 5.4,
       units = "in", dpi = 400)
p1+theme_blank_frame()
save_emf(p1+theme_blank_frame(),"Fig6A_hallmark_1.emf", width = 1.5, height = 9)

# ---- source data (supplementary table) ----------------------
fwrite(d[order(-fgsea_NES),
         .(pathway, setSize, NES = fgsea_NES, fgsea_FDR = fgsea_padj,
           camera_direction = camera_dir, camera_FDR, inter_gene_cor)],
       "figures/Fig6A_source_data.csv")
