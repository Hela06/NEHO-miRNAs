# =============================================================================
# Figure 4 — Venn diagrams
# -----------------------------------------------------------------------------
# Description:
#   1. Derives p-values from BWQS credible intervals for volcano-style plots
#   2. Builds Euler/Venn diagrams comparing miRNAs significant in:
#      - BWQS mixture analyses (EEs, OCs, TOT)
#      - Individual pollutant linear regressions (FDR < 0.05)
#   Diagrams are produced separately for maternal and cordonal, then combined.
#
# Input:
#   - output_bWQS/FINAL_MODELS_significant_TOT.xlsx     (output of 2_1_BWQS.R)
#   - output_bWQS/FINAL_MODELS_significant_divided.xlsx (output of 2_1_BWQS.R)
#   - output_bWQS/FINAL_MODELS_nonsignificant.xlsx      (output of 2_1_BWQS.R)
#   - regression_MM.xlsx                                (output of 2_0)
#   - regression_CC.xlsx                                (output of 2_0)
#
# Output:
#   - Figure4_Venn.tiff
#
# Author: Ilaria Cosentini
# Date:   June 2026 (last update)
# =============================================================================


# --- Libraries ---------------------------------------------------------------
library(readxl)
library(dplyr)
library(data.table)
library(eulerr)
library(RColorBrewer)
library(ggpubr)
library(cowplot)


# --- Load data ---------------------------------------------------------------
bwqs_TOT  <- as.data.frame(readxl::read_excel("./output_bWQS/FINAL_MODELS_significant_TOT.xlsx"))
bwqs_div  <- as.data.frame(readxl::read_excel("./output_bWQS/FINAL_MODELS_significant_divided.xlsx"))
bwqs_insg <- as.data.frame(readxl::read_excel("./output_bWQS/FINAL_MODELS_nonsignificant.xlsx"))
df_MM     <- as.data.frame(readxl::read_excel("./regression_MM.xlsx"))
df_CC     <- as.data.frame(readxl::read_excel("./regression_CC.xlsx"))


# =============================================================================
# SECTION 1: COMPUTE P-VALUES FROM BWQS CREDIBLE INTERVALS
# =============================================================================

compute_p_value <- function(lower, upper) {
  theta_hat <- (upper + lower) / 2
  SE        <- (upper - theta_hat) / 1.96
  Z         <- theta_hat / SE
  2 * (1 - pnorm(abs(Z)))
}

add_labels <- function(df) {
  df %>%
    mutate(
      p_value   = mapply(compute_p_value, `2.5%`, `95%`),
      MET_LIP   = ifelse(grepl("_CD$", miRNA),
                         paste0("C_WQS_MIX_", ANALYSIS),
                         paste0("M_WQS_MIX_", ANALYSIS)),
      miRNA_fin = gsub("_CD$", "", miRNA)
    )
}

# Keep only beta1 rows and add derived columns
bwqs_TOT$ANALYSIS  <- "TOT"
bwqs_TOT  <- add_labels(subset(bwqs_TOT,  mix_name == "beta1"))
bwqs_div  <- add_labels(subset(bwqs_div,  mix_name == "beta1"))
bwqs_insg <- add_labels(subset(bwqs_insg, mix_name == "beta1"))


# =============================================================================
# SECTION 2: ASSEMBLE df_bwqs_beta1
# =============================================================================

df_bwqs_beta1 <- bind_rows(bwqs_TOT, bwqs_div, bwqs_insg)

df_bwqs_beta1 <- df_bwqs_beta1 %>%
  mutate(
    ASSOCIATION = case_when(
      mean > 0 & p_value < 0.05 ~ "POSITIVE",
      mean < 0 & p_value < 0.05 ~ "NEGATIVE",
      TRUE                       ~ "NO ASSOCIATION"
    ),
    SIDE       = ifelse(grepl("^M_", MET_LIP), "MATERNAL", "CORD"),
    miRNA_annt = ifelse(p_value < 0.05, miRNA_fin, NA)
  )

df_bwqs_beta1$SIDE <- factor(df_bwqs_beta1$SIDE, levels = c("MATERNAL", "CORD"))


# =============================================================================
# SECTION 3: EULER / VENN DIAGRAMS
# =============================================================================

# Helper to extract significant miRNAs from BWQS
sig_bwqs <- function(side, analysis) {
  unique(df_bwqs_beta1$miRNA_fin[
    df_bwqs_beta1$p_value < 0.05 &
    df_bwqs_beta1$SIDE    == side &
    df_bwqs_beta1$ANALYSIS == analysis
  ])
}

# Helper to extract significant miRNAs from linear regression
sig_lm <- function(df, pollutant_pattern) {
  ids <- df$miRNA[df$p.adjust.FDR < 0.05 & grepl(pollutant_pattern, df$pollutant)]
  unique(gsub("_CD$", "", ids))
}

# Color palettes
colors_EEs  <- brewer.pal(4, "Pastel1")
colors_PCBs <- c("#CCEBC5", "#DECBE4", "#D9D9D9", "#FFED6F")
colors_WQS  <- c(colors_PCBs[1], "#B3CDE3", colors_EEs[1])


# --- Maternal EEs ------------------------------------------------------------
list_M_EEs <- list(
  EEs = sig_bwqs("MATERNAL", "EEs"),
  Cu  = sig_lm(df_MM, "M_Cu$"),
  Se  = sig_lm(df_MM, "M_Se$"),
  Zn  = sig_lm(df_MM, "M_Zn$")
)
a_p <- plot(euler(list_M_EEs),
            quantities = list(type = "counts", cex = 2),
            labels = list(cex = 3), fills = colors_EEs, edges = FALSE)

# --- Cordonal EEs ------------------------------------------------------------
list_C_EEs <- list(
  EEs = sig_bwqs("CORD", "EEs"),
  Cu  = sig_lm(df_CC, "C_Cu$"),
  Se  = sig_lm(df_CC, "C_Se$"),
  Zn  = sig_lm(df_CC, "C_Zn$")
)
b_p <- plot(euler(list_C_EEs),
            quantities = list(type = "counts", cex = 2),
            labels = list(cex = 3), fills = colors_EEs, edges = FALSE)

# --- Maternal POPs -----------------------------------------------------------
list_M_PCBs <- list(
  POPs = sig_bwqs("MATERNAL", "OCs"),
  HCB  = sig_lm(df_MM, "M_HCB$"),
  DDE  = sig_lm(df_MM, "M_DDE$"),
  PCBs = sig_lm(df_MM, "M_PCB")
)
c_p <- plot(euler(list_M_PCBs),
            quantities = list(type = "counts", cex = 2),
            labels = list(cex = 3), fills = colors_PCBs, edges = FALSE)

# --- Cordonal POPs -----------------------------------------------------------
list_C_PCBs <- list(
  POPs = sig_bwqs("CORD", "OCs"),
  HCB  = sig_lm(df_CC, "C_HCB$"),
  DDE  = sig_lm(df_CC, "C_DDE$"),
  PCBs = sig_lm(df_CC, "C_PCB")
)
d_p <- plot(euler(list_C_PCBs),
            quantities = list(type = "counts", cex = 2),
            labels = list(cex = 3), fills = colors_PCBs, edges = FALSE)

# --- Maternal WQS mixture overlap --------------------------------------------
list_M_WQS <- list(
  POPs = sig_bwqs("MATERNAL", "OCs"),
  TOT  = sig_bwqs("MATERNAL", "TOT"),
  EEs  = sig_bwqs("MATERNAL", "EEs")
)
e_p <- plot(euler(list_M_WQS),
            quantities = list(type = "counts", cex = 2),
            labels = list(cex = 3), fills = colors_WQS, edges = FALSE)

# --- Cordonal WQS mixture overlap --------------------------------------------
list_C_WQS <- list(
  POPs = sig_bwqs("CORD", "OCs"),
  TOT  = sig_bwqs("CORD", "TOT"),
  EEs  = sig_bwqs("CORD", "EEs")
)
f_p <- plot(euler(list_C_WQS),
            quantities = list(type = "counts", cex = 2),
            labels = list(cex = 3), fills = colors_WQS, edges = FALSE)


# =============================================================================
# SECTION 4: EXPORT VENN INPUT DATA
# =============================================================================
# For each diagram: one row per miRNA, boolean columns indicating membership.

venn_to_df <- function(lst) {
  all_mirna <- unique(unlist(lst))
  out <- data.frame(miRNA = all_mirna)
  for (nm in names(lst)) {
    out[[nm]] <- out$miRNA %in% lst[[nm]]
  }
  out[do.call(order, lapply(names(lst), function(nm) !out[[nm]])), ]
}

writexl::write_xlsx(
  list(
    M_EEs  = venn_to_df(list_M_EEs),
    C_EEs  = venn_to_df(list_C_EEs),
    M_POPs = venn_to_df(list_M_PCBs),
    C_POPs = venn_to_df(list_C_PCBs),
    M_WQS  = venn_to_df(list_M_WQS),
    C_WQS  = venn_to_df(list_C_WQS)
  ),
  "Figure4_Venn_input.xlsx"
)

cat("Venn input data saved to Figure4_Venn_input.xlsx\n")


# =============================================================================
# SECTION 5: COMBINE AND SAVE
# =============================================================================

fig <- ggarrange(a_p, c_p, e_p, b_p, d_p, f_p,
                 labels        = c("A", "B", "C", "D", "E", "F"),
                 ncol          = 3,
                 nrow          = 2,
                 font.label    = list(size = 30))

label_maternal <- ggdraw() +
  draw_label("Maternal", angle = 90, fontface = "bold", size = 30, x = 0.6) +
  draw_line(x = c(0.9, 0.9), y = c(0.05, 0.95), linewidth = 1.5)

label_cord <- ggdraw() +
  draw_label("Cord", angle = 90, fontface = "bold", size = 30, x = 0.6) +
  draw_line(x = c(0.9, 0.9), y = c(0.05, 0.95), linewidth = 1.5)

left_labels <- plot_grid(label_maternal, label_cord, ncol = 1)
final_fig   <- plot_grid(left_labels, fig, ncol = 2, rel_widths = c(0.05, 1))

tiff("Figure4_Venn.tiff", res = 300, width = 30, height = 20, units = "in")
print(final_fig)
dev.off()

cat("Figure saved to Figure4_Venn.tiff\n")
