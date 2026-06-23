# =============================================================================
# Supplementary Figures — Regression heatmaps and combined WQS + Venn panels
# -----------------------------------------------------------------------------
# Description:
#   Produces four supplementary figures:
#
#   FigS1 :
#     Heatmap of significant linear regression estimates for maternal miRNA
#     vs maternal chemicals.
#   FigS2 :
#     Heatmap of significant linear regression estimates for cordonal miRNA vs 
#     cordonal chemicals
#
#   FigS3:
#     Heatmap of BWQS beta1 (maternal mixtures).
#
#   FigS4:
#     Heatmap of BWQS beta1 (cordonal mixtures).
#
#   
#
#   Also exports Table S3 and Table S4 as Excel files.
#
# Input:
#   - regression_MM.xlsx                          (output of 2_0)
#   - regression_CC.xlsx                          (output of 2_0)
#   - output_bWQS/WQS_all_results.xlsx            (output of 2_1)
#   
#
#
# Author: Ilaria Cosentini
# Date:   June 2026 (last update)
# =============================================================================


# --- Libraries ---------------------------------------------------------------
library(readxl)
library(writexl)
library(ggplot2)
library(scales)
library(ggpubr)
library(data.table)


# --- Load data ---------------------------------------------------------------
df_MM   <- as.data.frame(readxl::read_excel("./regression_MM.xlsx"))
df_CC   <- as.data.frame(readxl::read_excel("./regression_CC.xlsx"))
WQS_all <- as.data.frame(readxl::read_excel("./output_bWQS/WQS_all_results.xlsx"))



# =============================================================================
# SECTION 1: SHARED THEME AND HELPERS
# =============================================================================

heatmap_theme <- theme_minimal(base_size = 15) +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 17),
    axis.text.y        = element_text(size = 17),
    panel.grid.major   = element_line(linetype = "dashed", color = "gray80"),
    axis.line          = element_line(color = "black"),
    axis.ticks         = element_line(color = "black"),
    legend.position    = "right",
    legend.key.height  = unit(0.1, "npc")
  )

fill_scale <- scale_fill_gradient2(
  low      = "#4575b4",
  mid      = "white",
  high     = "#d73027",
  midpoint = 0,
  limits   = c(-0.6, 0.6),
  oob      = scales::squish
)

# Chemical name cleaner
clean_name <- function(met_lip, prefix) {
  nm <- ifelse(grepl(paste0("^", prefix), met_lip),
               gsub(paste0("^", prefix), "", met_lip),
               "PCBs-SUM")
  ifelse(met_lip %in% c("SUM_PCBind", "SUM_C_PCBind"), "PCBs-SUM", nm)
}

chemical_levels <- c("PCBs-SUM", "PCB187", "PCB180", "PCB170",
                     "PCB153", "PCB138", "PCB118", "PCB74",
                     "HCB", "DDE", "Zn", "Se", "Cu")


# =============================================================================
# SECTION 2: LINEAR REGRESSION HEATMAPS (FigS1)
# =============================================================================

# --- Maternal ----------------------------------------------------------------
df_prova_MM <- df_MM[df_MM$p.adjust.FDR < 0.05, ]
df_prova_MM <- df_prova_MM[!df_prova_MM$pollutant %in%
                              c("RAT_Cu_Zn", "M_PCB156", "M_PCB183",
                                "M_PCB187", "M_PCB170", "SUM_PCB"), ]

df_prova_MM$NAME  <- factor(clean_name(df_prova_MM$pollutant, "M_"),
                             levels = chemical_levels)
df_prova_MM$value <- factor(gsub("_CD$", "", df_prova_MM$miRNA),
                             levels = c("hsa-miR-1-3p", "hsa-miR-10b-5p",
                                        "hsa-miR-30d-5p", "hsa-miR-30e-5p",
                                        "hsa-miR-122-5p", "hsa-miR-133a-3p",
                                        "hsa-miR-148a-3p", "hsa-miR-195-5p",
                                        "hsa-miR-205-5p", "hsa-miR-206",
                                        "hsa-miR-214-3p", "hsa-miR-374a-5p"))

fig_m_sng <- ggplot(df_prova_MM, aes(x = value, y = NAME, fill = estimate)) +
  geom_tile(color = "gray80") +
  geom_text(aes(label = ifelse(is.na(estimate), "", sprintf("%.2f", estimate))),
            size = 7) +
  fill_scale +
  labs(x = "MATERNAL-miRNAs", y = "MATERNAL-CHEMICALS", fill = "") +
  heatmap_theme


# --- Cordonal ----------------------------------------------------------------
df_prova_CC <- df_CC[df_CC$p.adjust.FDR < 0.05, ]
df_prova_CC <- df_prova_CC[df_prova_CC$pollutant != "RAT_C_Cu_Zn", ]

df_prova_CC$NAME  <- factor(clean_name(df_prova_CC$pollutant, "C_"),
                             levels = chemical_levels)
df_prova_CC$value <- factor(gsub("_CD$", "", df_prova_CC$miRNA),
                             levels = c("hsa-miR-7-5p", "hsa-miR-10a-5p",
                                        "hsa-miR-16-5p", "hsa-miR-22-3p",
                                        "hsa-miR-27a-3p", "hsa-miR-29a-3p",
                                        "hsa-miR-30d-5p", "hsa-miR-122-5p",
                                        "hsa-miR-133a-3p", "hsa-miR-133b",
                                        "hsa-miR-134-5p", "hsa-miR-195-5p",
                                        "hsa-miR-200c-3p", "hsa-miR-376c-3p",
                                        "hsa-miR-574-3p"))

fig_c_sng <- ggplot(df_prova_CC, aes(x = value, y = NAME, fill = estimate)) +
  geom_tile(color = "gray80") +
  geom_text(aes(label = ifelse(is.na(estimate), "", sprintf("%.2f", estimate))),
            size = 7) +
  fill_scale +
  labs(x = "CORD-miRNAs", y = "CORD-CHEMICALS", fill = "") +
  heatmap_theme


# --- Save FigS1/S2 and tables ---------------------------------------------------
tiff("FigureS1_lm_heatmap.tiff", res = 300, width = 25, height = 20, units = "in")
print(fig_m_sng)
dev.off()


tiff("FigureS2_lm_heatmap.tiff", res = 300, width = 25, height = 20, units = "in")
print(fig_c_sng)
dev.off()

# Table S3 — maternal regression
tableS3 <- df_prova_MM
tableS3$beta   <- paste0(round(tableS3$estimate, 2), " (",
                          round(tableS3$conf.low,  2), "; ",
                          round(tableS3$conf.high, 2), ")")
tableS3$Pvalue <- round(tableS3$p.adjust.FDR, 3)
writexl::write_xlsx(tableS3, "TableS3.xlsx")

# Table S4 — cordonal regression
tableS4 <- df_prova_CC
tableS4$beta   <- paste0(round(tableS4$estimate, 2), " (",
                          round(tableS4$conf.low,  2), "; ",
                          round(tableS4$conf.high, 2), ")")
tableS4$Pvalue <- round(tableS4$p.adjust.FDR, 3)
writexl::write_xlsx(tableS4, "TableS4_lm_cordonal.xlsx")


# =============================================================================
# SECTION 3: BWQS HEATMAPS
# =============================================================================

# Helper: build BWQS heatmap and mark miRNAs also significant in single LM
make_wqs_heatmap <- function(wqs_df, side_prefix, lm_df, ee_cols, x_label, y_label, mirna_levels) {

  df <- wqs_df[grepl(paste0("^", side_prefix, "_WQS_MIX"), wqs_df$MET_LIP), ]
  df <- na.omit(df)
  df$ANALYSIS <- ifelse(df$ANALYSIS == "OCs", "POPs", df$ANALYSIS)
  df$ANALYSIS <- factor(df$ANALYSIS, levels = c("TOT", "POPs", "EEs"))

  # Mark miRNAs also significant in individual LM (* = confirmed by single regression)
  ee_sig  <- intersect(df$miRNA_fin[df$ANALYSIS == "EEs"],
                       lm_df$value[lm_df$pollutant %in% ee_cols])
  pop_sig <- intersect(df$miRNA_fin[df$ANALYSIS == "POPs"],
                       lm_df$value[!lm_df$pollutant %in% ee_cols])

  df$in_sng <- ""
  df$in_sng[df$miRNA_fin %in% ee_sig  & df$ANALYSIS == "EEs"]  <- "*"
  df$in_sng[df$miRNA_fin %in% pop_sig & df$ANALYSIS == "POPs"] <- "*"

  df$miRNA_fin <- factor(df$miRNA_fin, levels = mirna_levels)

  ggplot(df, aes(x = miRNA_fin, y = ANALYSIS, fill = mean)) +
    geom_tile(color = "gray80") +
    geom_text(aes(label = ifelse(is.na(mean), "",
                                 paste0(sprintf("%.2f", mean), in_sng))),
              size = 7) +
    fill_scale +
    labs(x = x_label, y = y_label, fill = "") +
    heatmap_theme
}

fig_m_wqs <- make_wqs_heatmap(
  wqs_df       = WQS_all,
  side_prefix  = "M",
  lm_df        = df_prova_MM,
  ee_cols      = c("M_Se", "M_Cu", "M_Zn"),
  x_label      = "MATERNAL-miRNAs",
  y_label      = "MATERNAL-MIXTURE",
  mirna_levels = c("hsa-miR-1-3p", "hsa-miR-7-5p", "hsa-miR-10b-5p",
                   "hsa-miR-15a-5p", "hsa-miR-17-5p", "hsa-miR-18a-5p",
                   "hsa-miR-20a-5p", "hsa-miR-25-3p", "hsa-miR-27a-3p",
                   "hsa-miR-29a-3p", "hsa-miR-30d-5p", "hsa-miR-31-5p",
                   "hsa-miR-100-5p", "hsa-miR-96-5p", "hsa-miR-122-5p",
                   "hsa-miR-126-3p", "hsa-miR-128-3p", "hsa-miR-130b-3p",
                   "hsa-miR-133a-3p", "hsa-miR-133b", "hsa-miR-134-5p",
                   "hsa-miR-143-3p", "hsa-miR-145-5p", "hsa-miR-193a-5p",
                   "hsa-miR-195-5p", "hsa-miR-200a-3p", "hsa-miR-200b-3p",
                   "hsa-miR-203a-3p", "hsa-miR-205-5p", "hsa-miR-210-3p",
                   "hsa-miR-214-3p", "hsa-miR-296-5p", "hsa-miR-374a-5p",
                   "hsa-miR-375-3p", "hsa-miR-376c-3p")
)

fig_c_wqs <- make_wqs_heatmap(
  wqs_df       = WQS_all,
  side_prefix  = "C",
  lm_df        = df_prova_CC,
  ee_cols      = c("C_Se", "C_Cu", "C_Zn"),
  x_label      = "CORD-miRNAs",
  y_label      = "CORD-MIXTURE",
  mirna_levels = c("hsa-miR-7-5p", "hsa-miR-9-5p", "hsa-miR-16-5p",
                   "hsa-miR-17-3p", "hsa-miR-19b-3p", "hsa-miR-22-3p",
                   "hsa-miR-25-3p", "hsa-miR-27a-3p", "hsa-miR-29a-3p",
                   "hsa-miR-30d-5p", "hsa-miR-34a-5p", "hsa-miR-122-5p",
                   "hsa-miR-124-3p", "hsa-miR-125b-5p", "hsa-miR-133a-3p",
                   "hsa-miR-133b", "hsa-miR-134-5p", "hsa-miR-141-3p",
                   "hsa-miR-143-3p", "hsa-miR-145-5p", "hsa-miR-150-5p",
                   "hsa-miR-191-5p", "hsa-miR-195-5p", "hsa-miR-200a-3p",
                   "hsa-miR-205-5p", "hsa-miR-206", "hsa-miR-374a-5p",
                   "hsa-miR-375-3p", "hsa-miR-376c-3p", "hsa-miR-499a-5p")
)

# FigS3 — maternal WQS heatmap 
tiff("FigureS3_wqs_maternal.tiff", res = 300, width = 30, height = 15, units = "in")
print(fig_m_wqs)
dev.off()

# FigS4 — cordonal WQS heatmap 
tiff("FigureS3_wqs_cordonal.tiff", res = 300, width = 30, height = 15, units = "in")
print(fig_c_wqs)
dev.off()

message("All supplementary figures saved.")
