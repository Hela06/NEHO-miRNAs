# =============================================================================
# Figure 1 — Placental Transfer Boxplot
# -----------------------------------------------------------------------------
# Description:
#   Visualises the placental transfer ratio (cord / maternal) for selected
#   chemicals as a horizontal boxplot with jittered points.
#   Values > 0 indicate higher cordonal than maternal concentration.
#   The x-axis is symmetric around 0 and formatted as percentage.
#
# Input:
#   - INPUT_FIGURE1.xlsx   (one row per subject per chemical;
#                           columns: CHEMICAL, value, centered)
#
# Output:
#   - Figure1_placental_transfer.tiff
#
# Author: Ilaria Cosentini
# Date:   June 2026 (last update)
# =============================================================================


# --- Libraries ---------------------------------------------------------------
library(readxl)
library(ggplot2)
library(dplyr)


# --- Load data ---------------------------------------------------------------
input_plactr <- as.data.frame(readxl::read_excel("INPUT_FIGURE1.xlsx"))


# --- Prepare variables -------------------------------------------------------
chemical_order <- c("Zn", "Se", "Cu", "HCB", "DDE", "PCB138", "PCB153", "PCB180")

input_plactr$CHEMICAL   <- factor(input_plactr$CHEMICAL, levels = chemical_order)
input_plactr$centered_2 <- ifelse(
  input_plactr$centered < 0,
  input_plactr$value * 100 * -1,
  input_plactr$value * 100
)


# --- Plot --------------------------------------------------------------------
p <- ggplot(input_plactr, aes(x = centered_2, y = CHEMICAL)) +

  geom_boxplot(
    fill          = "#4a6fa5",
    color         = "#4a6fa5",
    outlier.shape = NA,
    alpha         = 0.7,
    linewidth     = 0.5,
    width         = 0.6
  ) +

  stat_boxplot(
    geom      = "errorbar",
    width     = 0,
    color     = "white",
    linewidth = 0
  ) +

  geom_jitter(
    color  = "#5ba3d9",
    alpha  = 0.6,
    size   = 1.8,
    height = 0.25,
    width  = 0
  ) +

  geom_vline(xintercept = 0, color = "red", linewidth = 0.4, linetype = "dashed") +

  scale_x_continuous(
    breaks = seq(-100, 120, by = 20),
    labels = function(x) paste0(abs(x), "%"),
    limits = c(-105, 220)
  ) +

  labs(x = "PLACENTAL TRANSFER (CORD/MATERNAL)", y = NULL) +

  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    axis.ticks       = element_blank(),
    plot.margin      = margin(20, 30, 20, 20)
  )


# --- Save --------------------------------------------------------------------
ggsave(
  filename = "Figure1_placental_transfer.tiff",
  plot     = p,
  width    = 10,
  height   = 6,
  dpi      = 330,
  device   = "tiff"
)

message("Figure saved to Figure1_placental_transfer.tiff")
