# =============================================================================
# Figure 2 — Correlation plots and scatterplots
# -----------------------------------------------------------------------------
# Description:
#   Three correlation panels (maternal, cordonal, maternal vs cordonal diagonal)
#   assembled with a shared colorbar, plus a faceted scatterplot of key
#   maternal vs cordonal analyte pairs.
#
#   IMPORTANT — data privacy:
#     The raw pollutant dataset (pll_gr) contains patient-level data and
#     cannot be shared publicly. This script exports the correlation matrices
#     and p-value matrices to Figure2_correlation_input.xlsx so that the
#     figures can be reproduced from those summary inputs.
#     To reproduce the figures from scratch, load the correlation matrices
#     from that file instead of recomputing from raw data.
#
# Input (raw, not publicly available):
#   - pll_gr: pollutant dataset (produced by 1_0_miRNA_normalization.R)
#
# Input (public, for reproduction):
#   - Figure2_correlation_input.xlsx   (exported by this script on first run)
#
# Output:
#   - Figure2_correlation_input.xlsx   (correlation + p-value matrices)
#   - Figure2_correlations.tiff
#
# Author: Ilaria Cosentini
# Date:   June 2026 (last update)
# =============================================================================


# --- Libraries ---------------------------------------------------------------
library(readxl)
library(writexl)
library(corrplot)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(cowplot)


# =============================================================================
# SECTION 1: CORRELATION MATRICES
# -----------------------------------------------------------------------------
# If you have access to the original data, run block 1a to compute and export.
# Otherwise, run block 1b to load from the exported file.
# =============================================================================

# ---- 1a: Compute from raw data (requires pll_gr) ----------------------------
# Uncomment this block if running with original data.

# pll_gr_M <- pll_gr[, c(1:9, 12)]
# colnames(pll_gr_M) <- gsub("^M_", "", colnames(pll_gr_M))
# data_crr_M  <- cor(as.matrix(pll_gr_M), method = "pearson")
# testRes_M   <- corrplot::cor.mtest(pll_gr_M, conf.level = 0.95)
#
# pll_gr_C <- na.omit(pll_gr[, 15:22])
# colnames(pll_gr_C) <- gsub("^C_", "", colnames(pll_gr_C))
# data_crr_C  <- cor(as.matrix(pll_gr_C), method = "pearson")
# testRes_C   <- corrplot::cor.mtest(pll_gr_C, conf.level = 0.95)
#
# # Maternal vs Cordonal: diagonal of cross-correlation matrix
# data_crr_MC <- cor(as.matrix(pll_gr_M),
#                    as.matrix(pll_gr_C[row.names(pll_gr_C) %in% row.names(pll_gr_M), ]),
#                    method = "pearson")
# testRes_MC  <- corrplot::cor.mtest(
#                  cbind(pll_gr_M, pll_gr_C[row.names(pll_gr_C) %in% row.names(pll_gr_M), ]),
#                  conf.level = 0.95)
#
# M_single_col     <- t(matrix(diag(data_crr_MC), nrow = 1,
#                              dimnames = list(NULL, colnames(data_crr_MC))))
# row.names(M_single_col) <- gsub("^M_","",row.names(data_crr_MC_sb))
#colnames(M_single_col) <- ""
#testRes_MC$p <- t(matrix(diag(testRes_MC$p), nrow = 1, dimnames = list(NULL, names(testRes_MC$p))))
#row.names(testRes_MC$p) <- gsub("^M_","",row.names(testRes_MC$p))
#M_single_col[c(2:4),1] <- NA
#
# # --- Export correlation matrices (public input) ---
# writexl::write_xlsx(
#   list(
#     corr_maternal     = as.data.frame(data_crr_M),
#     pval_maternal     = as.data.frame(testRes_M$p),
#     corr_cordonal     = as.data.frame(data_crr_C),
#     pval_cordonal     = as.data.frame(testRes_C$p),
#     corr_MvC_diag     = as.data.frame(data_crr_MC),
#     pval_MvC_diag     = as.data.frame(testRes_MC$p),
#     scatter_input     = pll_gr[, c(1, 15, 5, 19, 8, 9, 12, 20:22)]
#   ),
#   "Figure2_correlation_input.xlsx"
# )
# message("Correlation matrices exported to Figure2_correlation_input.xlsx")




# ---- 1b: Load from exported file (public reproduction) ----------------------
corr_M      <- as.matrix(read_excel("Figure2_correlation_input.xlsx", sheet = "corr_maternal"))
pval_M      <- as.matrix(read_excel("Figure2_correlation_input.xlsx", sheet = "pval_maternal"))
corr_C      <- as.matrix(read_excel("Figure2_correlation_input.xlsx", sheet = "corr_cordonal"))
pval_C      <- as.matrix(read_excel("Figure2_correlation_input.xlsx", sheet = "pval_cordonal"))
M_single_col <- as.matrix(read_excel("Figure2_correlation_input.xlsx", sheet = "corr_MvC_diag"))
p_single_col <- as.matrix(read_excel("Figure2_correlation_input.xlsx", sheet = "pval_MvC_diag"))
df_scatter  <- as.data.frame(read_excel("Figure2_correlation_input.xlsx", sheet = "scatter_input"))

# Restore rownames (first column is the analyte name)
rownames(corr_M)       <- colnames(corr_M)
rownames(pval_M)       <- colnames(pval_M)
rownames(corr_C)       <- colnames(corr_C)
rownames(pval_C)       <- colnames(pval_C)
rownames(M_single_col) <- rownames(corr_M)[seq_len(nrow(M_single_col))]
colnames(M_single_col) <- ""

cor_palette <- colorRampPalette(c("blue", "white", "red"))(200)


# =============================================================================
# SECTION 2: CORRPLOT PANELS (saved as temp PNG for assembly)
# =============================================================================

make_corrplot_png <- function(corr_mat, pval_mat, title,
                              coef_col = "white", width = 16) {
  tmp <- tempfile(fileext = ".png")
  png(tmp, res = 300, units = "cm", width = width, height = 16)
  corrplot::corrplot(
    corr_mat, p.mat = pval_mat,
    method       = "color",
    insig        = "blank",
    addCoef.col  = coef_col,
    number.cex   = 1,
    type         = "upper",
    tl.cex       = 1.5,
    cl.cex       = 1,
    bg           = FALSE,
    title        = title,
    mar          = c(0, 0, 2, 0),
    cex.main     = 2,
    cl.pos       = "n",
    tl.srt       = 45,
    tl.col       = "black",
    col          = cor_palette
  )
  dev.off()
  tmp
}

tmp1 <- make_corrplot_png(corr_M,      pval_M, "MATERNAL SIDE")
tmp2 <- make_corrplot_png(corr_C,      pval_C, "CORD SIDE")
tmp3 <- make_corrplot_png(M_single_col, p_single_col, "MATERNALvsCORD SIDE", width = 16)


# =============================================================================
# SECTION 3: COLORBAR PANEL
# =============================================================================

tmp4 <- tempfile(fileext = ".png")
png(tmp4, res = 300, units = "cm", width = 4, height = 16)

p_legend <- ggplot(data.frame(x = 1, y = seq(-1, 1, length.out = 200)),
                   aes(x = x, y = y, fill = y)) +
  geom_tile() +
  scale_fill_gradientn(colours = c("blue", "white", "red"), limits = c(-1, 1)) +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.5),
                     labels = seq(-1, 1, by = 0.5),
                     position = "right") +
  theme_void() +
  theme(
    legend.position     = "none",
    axis.text.y.right   = element_text(size = 14, margin = margin(l = 5)),
    axis.ticks.y.right  = element_line(),
    plot.margin         = margin(0, 10, 20, 30)
  )
print(p_legend)
dev.off()


# =============================================================================
# SECTION 4: SCATTERPLOT
# =============================================================================

df_scatter$id <- row.names(df_scatter)
df_scatter    <- na.omit(df_scatter)

df_long_sc <- df_scatter %>%
  pivot_longer(cols = -id,
               names_to  = c("compartment", "analyte"),
               names_sep = "_",
               values_to = "value")

df_plot <- df_long_sc %>%
  pivot_wider(names_from = compartment, values_from = value) %>%
  rename(maternal = M, fetal = C)

remove_outliers <- function(x) {
  x > quantile(x, 0.25, na.rm = TRUE) - 1.5 * IQR(x, na.rm = TRUE) &
  x < quantile(x, 0.75, na.rm = TRUE) + 1.5 * IQR(x, na.rm = TRUE)
}

df_clean <- df_plot %>%
  group_by(analyte) %>%
  filter(remove_outliers(maternal), remove_outliers(fetal)) %>%
  ungroup()

df_clean$analyte <- factor(df_clean$analyte,
                           levels = c("Se", "DDE", "PCB138", "PCB153", "PCB180"))

SCATTER <- ggplot(df_clean, aes(maternal, fetal, color = analyte)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  stat_regline_equation(label.x.npc = "left", label.y.npc = 1,    size = 5) +
  stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~~~~")),
           label.x.npc = "left", label.y.npc = 0.95, size = 5) +
  facet_wrap(~ analyte, nrow = 1, scales = "free") +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Maternal concentration", y = "Cord blood concentration") +
  theme_bw(base_size = 12) +
  theme(
    strip.text    = element_text(size = 14, face = "bold"),
    panel.border  = element_rect(color = "black", fill = NA),
    legend.position = "none",
    panel.grid    = element_blank()
  )


# =============================================================================
# SECTION 5: ASSEMBLE AND SAVE
# =============================================================================

img1 <- ggdraw() + draw_image(tmp1) + draw_label("A", x = 0.02, y = 0.97, fontface = "bold", size = 20)
img2 <- ggdraw() + draw_image(tmp2) + draw_label("B", x = 0.02, y = 0.97, fontface = "bold", size = 20)
img3 <- ggdraw() + draw_image(tmp3) + draw_label("C", x = 0.02, y = 0.97, fontface = "bold", size = 20)
img4 <- ggdraw() + draw_image(tmp4)

SCATTER_D <- ggdraw() +
  draw_plot(SCATTER) +
  draw_label("D", x = 0.02, y = 0.97, fontface = "bold", size = 20)

riga1 <- plot_grid(img1, img2, img3, img4, nrow = 1, rel_widths = c(1, 1, 1, 0.17))
riga2 <- plot_grid(SCATTER_D, nrow = 1)

figura_finale <- plot_grid(riga1, riga2, ncol = 1, rel_heights = c(1, 1))

tiff("Figure2_correlations.tiff", res = 300, units = "cm", width = 50, height = 30)
print(figura_finale)
dev.off()

# Clean up temp files
file.remove(tmp1, tmp2, tmp3, tmp4)

message("Figure saved to Figure2_correlations.tiff")
