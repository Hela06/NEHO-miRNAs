# =============================================================================
# Figure 6 — miRNA-Mixture Discrepancy Network
# -----------------------------------------------------------------------------
# Description:
#   Builds a bipartite network connecting mixture nodes (EEs, POPs) to miRNA
#   nodes, separately for maternal and cordonal sides. Edge color indicates
#   association direction (POS = red, NEG = blue). miRNA nodes are positioned
#   centrally; mixture nodes on the left (maternal) and right (cordonal).
#   A header bar labels each side.
#
#   Note: the TOT mixture is excluded from Figure 6 (included in Figure S5).
#   Table data corresponds to Table S6 in the manuscript.
#
# Input:
#   - hard-coded association table (data) derived from BWQS results
#
# Output:
#   - Figure6_network.tiff
#
# Author: Ilaria Cosentini
# Date:   June 2026 (last update)
# =============================================================================


# --- Libraries ---------------------------------------------------------------
library(dplyr)
library(tidyverse)
library(igraph)
library(ggraph)
library(ggforce)
library(data.table)


# =============================================================================
# SECTION 1: ASSOCIATION TABLE
# =============================================================================
# Each row = one miRNA-mixture association (EEs or POPs only, no TOT).
# MAT_MIX / CORD_MIX: mixture the miRNA is significant in (maternal/cordonal).
# MAT_DIR / CORD_DIR: direction of the BWQS beta1 effect (POS / NEG).

data <- tribble(
  ~miRNA,             ~MAT_MIX,  ~MAT_DIR,  ~CORD_MIX,  ~CORD_DIR,
  "hsa-miR-7-5p",     "EEs",     "POS",     "EEs",      "POS",
  "hsa-miR-7-5p",     NA,        NA,        "POPs",     "POS",
  "hsa-miR-25-3p",    "POPs",    "NEG",     "POPs",     "POS",
  "hsa-miR-27a-5p",   "EEs",     "NEG",     "POPs",     "NEG",
  "hsa-miR-29a-3p",   "POPs",    "NEG",     "EEs",      "NEG",
  "hsa-miR-29a-3p",   NA,        NA,        "POPs",     "NEG",
  "hsa-miR-30d-5p",   "EEs",     "NEG",     "EEs",      "POS",
  "hsa-miR-30d-5p",   NA,        NA,        "POPs",     "POS",
  "hsa-miR-122-5p",   "EEs",     "POS",     "EEs",      "POS",
  "hsa-miR-133a-3p",  "POPs",    "POS",     "POPs",     "NEG",
  "hsa-miR-133b",     "EEs",     "POS",     "POPs",     "NEG",
  "hsa-miR-134-5p",   "EEs",     "NEG",     "EEs",      "POS",
  "hsa-miR-143-3p",   "EEs",     "POS",     "POPs",     "NEG",
  "hsa-miR-145-5p",   "EEs",     "POS",     "POPs",     "NEG",
  "hsa-miR-195-5p",   "POPs",    "POS",     "EEs",      "POS",
  "hsa-miR-200a-3p",  "EEs",     "POS",     "POPs",     "NEG",
  "hsa-miR-205-5p",   "EEs",     "POS",     "POPs",     "NEG",
  "hsa-miR-374a-5p",  "EEs",     "NEG",     "POPs",     "POS",
  "hsa-miR-375-3p",   "EEs",     "POS",     "EEs",      "POS",
  "hsa-miR-376c-3p",  "EEs",     "NEG",     "EEs",      "POS"
)


# =============================================================================
# SECTION 2: BUILD EDGES
# =============================================================================

make_edges <- function(df, side) {
  df %>%
    filter(!is.na(.data[[paste0(side, "_MIX")]])) %>%
    transmute(
      from      = paste0(side, "_", .data[[paste0(side, "_MIX")]]),
      to        = paste0(miRNA, "_", side),
      DIRECTION = .data[[paste0(side, "_DIR")]]
    )
}

edges <- bind_rows(
  make_edges(data, "MAT"),
  make_edges(data, "CORD")
)

edges$ASSOCIATION <- ifelse(edges$DIRECTION == "POS", "red", "blue")


# =============================================================================
# SECTION 3: BUILD NODES
# =============================================================================

nodes <- tibble(name = unique(c(edges$from, edges$to))) %>%
  mutate(
    SIDE = case_when(
      str_detect(name, "^MAT_")  ~ "MATERNAL",
      str_detect(name, "^CORD_") ~ "CORD",
      str_detect(name, "_MAT$")  ~ "miRNA_M",
      str_detect(name, "_CORD$") ~ "miRNA_C",
      TRUE                       ~ ""
    ),
    TYPE = case_when(
      str_detect(name, "^MAT_")  ~ "MATERNAL",
      str_detect(name, "^CORD_") ~ "CORD",
      TRUE                       ~ "miRNA"
    ),
    label = gsub("MAT_|CORD_|_MAT$|_CORD$", "", name)
  )

# Attach direction to miRNA nodes (for fill color)
nodes <- nodes %>%
  left_join(edges[, c("to", "DIRECTION")], by = c("name" = "to")) %>%
  distinct()

# Hide label on maternal miRNA nodes (shown only on cordonal side)
nodes$label <- ifelse(nodes$SIDE == "miRNA_M", NA, nodes$label)


# =============================================================================
# SECTION 4: NODE ORDER AND LAYOUT
# =============================================================================

ordered_levels <- c(
  "CORD_EEs", "CORD_POPs",
  "hsa-miR-374a-5p_MAT",  "hsa-miR-374a-5p_CORD",
  "hsa-miR-205-5p_MAT",   "hsa-miR-205-5p_CORD",
  "hsa-miR-200a-3p_MAT",  "hsa-miR-200a-3p_CORD",
  "hsa-miR-145-5p_MAT",   "hsa-miR-145-5p_CORD",
  "hsa-miR-143-3p_MAT",   "hsa-miR-143-3p_CORD",
  "hsa-miR-133b_MAT",     "hsa-miR-133b_CORD",
  "hsa-miR-133a-3p_MAT",  "hsa-miR-133a-3p_CORD",
  "hsa-miR-25-3p_MAT",    "hsa-miR-25-3p_CORD",
  "hsa-miR-30d-5p_MAT",   "hsa-miR-30d-5p_CORD",
  "hsa-miR-376c-3p_MAT",  "hsa-miR-376c-3p_CORD",
  "hsa-miR-134-5p_MAT",   "hsa-miR-134-5p_CORD",
  "hsa-miR-195-5p_MAT",   "hsa-miR-195-5p_CORD",
  "hsa-miR-27a-5p_MAT",   "hsa-miR-27a-5p_CORD",
  "hsa-miR-29a-3p_MAT",   "hsa-miR-29a-3p_CORD",
  "hsa-miR-7-5p_MAT",     "hsa-miR-7-5p_CORD",
  "hsa-miR-375-3p_MAT",   "hsa-miR-375-3p_CORD",
  "hsa-miR-122-5p_MAT",   "hsa-miR-122-5p_CORD",
  "MAT_EEs", "MAT_POPs"
)

nodes <- nodes %>%
  mutate(name = factor(name, levels = ordered_levels)) %>%
  arrange(name) %>%
  mutate(name = as.character(name))


# --- Custom layout -----------------------------------------------------------
miRNA_nodes   <- nodes$name[nodes$TYPE == "miRNA"]
miRNA_M_nodes <- miRNA_nodes[str_detect(miRNA_nodes, "_MAT$")]
miRNA_C_nodes <- miRNA_nodes[str_detect(miRNA_nodes, "_CORD$")]

n_mirna  <- length(miRNA_M_nodes)
y_coords <- seq(-n_mirna / 2, n_mirna / 2, length.out = n_mirna)

layout <- data.frame(name = nodes$name, x = NA_real_, y = NA_real_)

# Central column: miRNA nodes
layout$x[layout$name %in% miRNA_nodes]    <- 0
layout$y[layout$name %in% miRNA_M_nodes]  <- y_coords
layout$y[layout$name %in% miRNA_C_nodes]  <- y_coords

# Left: maternal mixture nodes and miRNA_M offset
mat_nodes  <- nodes$name[str_detect(nodes$name, "MAT")]
mix_M      <- mat_nodes[str_detect(mat_nodes, "^MAT_")]
mirna_M_off <- mat_nodes[!str_detect(mat_nodes, "^MAT_")]
layout$x[layout$name %in% mix_M]       <- -5
layout$x[layout$name %in% mirna_M_off] <- -0.15
layout$y[layout$name %in% mix_M]       <- c(3, -3)[seq_along(mix_M)]

# Right: cordonal mixture nodes and miRNA_C offset
cord_nodes  <- nodes$name[str_detect(nodes$name, "CORD")]
mix_C       <- cord_nodes[str_detect(cord_nodes, "^CORD_")]
mirna_C_off <- cord_nodes[!str_detect(cord_nodes, "^CORD_")]
layout$x[layout$name %in% mix_C]        <- 5
layout$x[layout$name %in% mirna_C_off]  <- 0.15
layout$y[layout$name %in% mix_C]        <- c(3, -3)[seq_along(mix_C)]

# Align layout to vertex order in graph
g              <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
vertex_order   <- V(g)$name
layout_aligned <- layout[match(vertex_order, layout$name), ]


# =============================================================================
# SECTION 5: PLOT
# =============================================================================

nodes_with_coords <- nodes %>% left_join(layout, by = "name")

x_min_mat  <- min(layout$x[nodes$SIDE %in% c("MATERNAL", "miRNA_M")], na.rm = TRUE)
x_max_mat  <- max(layout$x[nodes$SIDE %in% c("MATERNAL", "miRNA_M")], na.rm = TRUE)
x_min_cord <- min(layout$x[nodes$SIDE %in% c("CORD", "miRNA_C")],     na.rm = TRUE)
x_max_cord <- max(layout$x[nodes$SIDE %in% c("CORD", "miRNA_C")],     na.rm = TRUE)
y_top      <- max(layout$y, na.rm = TRUE) + 1.5

graph_dsc <- ggraph(g, layout = "manual",
                    x = layout_aligned$x, y = layout_aligned$y) +
  geom_edge_link(aes(color = DIRECTION), linewidth = 0.8) +
  geom_node_point(data = nodes_with_coords,
                  aes(shape = TYPE,
                      fill  = ifelse(TYPE == "miRNA", DIRECTION, TYPE)),
                  size = 15, color = "black", show.legend = FALSE) +
  geom_node_text(aes(label = label), repel = FALSE, size = 7,
                 nudge_y = ifelse(nodes$TYPE %in% c("CORD", "MATERNAL"), 0.5, 1),
                 nudge_x = ifelse(nodes$TYPE == "CORD",  0.5,
                           ifelse(nodes$TYPE == "miRNA", -0.05, -0.5))) +
  scale_shape_manual(values = c(miRNA = 22, MATERNAL = 21, CORD = 21)) +
  scale_fill_manual(values  = c(miRNA = "white", CORD = "lightblue",
                                MATERNAL = "lightpink",
                                POS = "#d73027", NEG = "#4575b4")) +
  scale_edge_colour_manual(values = c(POS = "#d73027", NEG = "#4575b4")) +
  annotate("segment",
           x = x_min_mat,  xend = x_max_mat,
           y = y_top,      yend = y_top,
           linewidth = 2,  color = "black") +
  annotate("text",
           x = (x_min_mat + x_max_mat) / 2, y = y_top + 0.5,
           label = "MATERNAL", fontface = "bold", color = "black", size = 7) +
  annotate("segment",
           x = x_min_cord, xend = x_max_cord,
           y = y_top,      yend = y_top,
           linewidth = 2,  color = "black") +
  annotate("text",
           x = (x_min_cord + x_max_cord) / 2, y = y_top + 0.5,
           label = "CORD", fontface = "bold", color = "black", size = 7) +
  theme_void() +
  theme(legend.title = element_text(size = 20, face = "bold"),
        legend.text  = element_text(size = 18))


# =============================================================================
# SECTION 6: SAVE
# =============================================================================

tiff("Figure6_network.tiff", res = 300, width = 20, height = 15, units = "in")
print(graph_dsc)
dev.off()

cat("Figure saved to Figure6_network.tiff\n")
