# =============================================================================
# Figure 7 — Pathway Enrichment Analysis and Heatmap
# -----------------------------------------------------------------------------
# Description:
#   Two-step analysis:
#
#   PART 1:
#     - Computes pairwise pathway similarity from enrichment results
#     - Clusters pathways into communities (Louvain) using Jaccard similarity
#     - Maps Reactome hierarchy to identify hub pathways per community
#     -> produces: ptw_hr_3 (pathway groups with hub labels)
#
#   PART 2:
#     - Retrieves validated miRNA targets via multiMiR
#     - Filters targets confirmed in >= 3 databases
#     - Runs compareCluster pathway enrichment (Reactome)
#     - Builds binary heatmap: pathway group x BWQS cluster (UP/DOWN)
#     -> produces: Figure7_enrichment.tiff
#
# Input:
#   - output_bWQS/WQS_all_results.xlsx          (output of 2_1_BWQS.R)
#   - Reactome_pathways_March2026.txt            (Reactome hierarchy file)
#   - ReactomePathways_list_March2026.txt        (Reactome pathway names)
#
# Output:
#   - multimir_results_raw.rds
#   - RESULT_ENRICHMENT_CLUSTERPROF.xlsx
#   - RESULT_ENRICHMENT_CLUSTERPROF_wMiRNA.xlsx
#   - pathway_communities_and_degree.csv
#   - graph_pathways_grouping_wqs.tiff
#   - Figure7_enrichment.tiff
#
# Author: Ilaria Cosentini
# Date:   June 2026 (last update)
# =============================================================================


# --- Libraries ---------------------------------------------------------------
library(readxl)
library(writexl)
library(dplyr)
library(tidyr)
library(data.table)
library(multiMiR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(igraph)
library(ggraph)
library(ComplexHeatmap)
library(circlize)
library(grid)


# --- Load BWQS results -------------------------------------------------------
WQS_all <- as.data.frame(readxl::read_excel("./output_bWQS/WQS_all_results.xlsx"))

# Subset to significant effects (credible interval excludes zero)
WQS_sig <- WQS_all[!is.na(WQS_all$mean_flagged), c("miRNA", "ANALYSIS", "MET_LIP", "miRNA_fin", "mean_flagged")]
WQS_sig$GROUP      <- ifelse(grepl("_CD$", WQS_sig$miRNA), "CORD", "MATERNAL")
WQS_sig$ANALYSIS_2 <- ifelse(WQS_sig$mean_flagged < 0,
                              paste0(WQS_sig$ANALYSIS, "_DOWN"),
                              paste0(WQS_sig$ANALYSIS, "_UP"))
WQS_sig$mean_flagged <- NULL


# =============================================================================
# PART 1: MULTIMIR TARGET RETRIEVAL AND PATHWAY ENRICHMENT
# =============================================================================

# --- 1.1 Retrieve validated targets ------------------------------------------
mir_sign  <- unique(WQS_sig$miRNA_fin)
multi_mir <- get_multimir(org = "hsa", mirna = mir_sign, table = "validated", summary = TRUE)
saveRDS(multi_mir, "multimir_results_raw.rds")
# multi_mir <- readRDS("multimir_results_raw.rds")  # reload if needed

df_mir    <- multi_mir@data
df_mir_sb <- unique(df_mir[, c("mature_mirna_acc", "mature_mirna_id", "target_ensembl")])

universe_target <- unique(df_mir_sb$target_ensembl)


# --- 1.2 Keep targets confirmed in >= 3 databases ----------------------------
database_counts <- df_mir_sb %>%
  group_by(mature_mirna_acc, target_ensembl) %>%
  summarise(n_databases = n_distinct(df_mir$database[
    df_mir$mature_mirna_acc == cur_group()$mature_mirna_acc &
    df_mir$target_ensembl   == cur_group()$target_ensembl
  ]), .groups = "drop")

# Simpler equivalent using the full df_mir
database_counts <- df_mir %>%
  group_by(mature_mirna_acc, target_ensembl) %>%
  summarise(n_databases = n_distinct(database), .groups = "drop")

sel_miR_sb <- subset(database_counts, n_databases >= 3 & target_ensembl != "")

enr_WQS <- merge(sel_miR_sb,
                 unique(df_mir[, c("mature_mirna_acc", "mature_mirna_id",
                                   "target_ensembl", "target_entrez", "target_symbol")]),
                 by = c("mature_mirna_acc", "target_ensembl"), all.x = TRUE)
enr_WQS <- enr_WQS[!duplicated(enr_WQS), ]


# --- 1.3 Map targets to Entrez IDs -------------------------------------------
input_comp_cl <- merge(WQS_sig, enr_WQS,
                       by.x = "miRNA_fin", by.y = "mature_mirna_id", all.x = TRUE)

conversion          <- bitr(input_comp_cl$target_ensembl,
                            fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
conversion_universe <- bitr(universe_target,
                            fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

input_comp_cl <- merge(input_comp_cl, conversion,
                       by.x = "target_ensembl", by.y = "ENSEMBL", all.x = TRUE)

# Use bitr-derived ENTREZID where target_entrez is missing
input_comp_cl$target_entrez <- ifelse(
  is.na(input_comp_cl$target_entrez) | input_comp_cl$target_entrez == "",
  input_comp_cl$ENTREZID,
  input_comp_cl$target_entrez
)
input_comp_cl$ENTREZID <- NULL


# --- 1.4 Pathway enrichment (all BWQS analyses) ------------------------------
formula_res_rev <- clusterProfiler::compareCluster(
  target_entrez ~ GROUP + ANALYSIS_2,
  data     = input_comp_cl,
  fun      = "enrichPathway",
  universe = conversion_universe$ENTREZID
)

pathways_rev <- as.data.frame(formula_res_rev)
writexl::write_xlsx(pathways_rev, "RESULT_ENRICHMENT_CLUSTERPROF.xlsx")


# --- 1.5 Map Entrez IDs back to miRNA names ----------------------------------
replace_entrez_with_mirna <- function(entrez_string, mirna_df) {
  entrez_list <- unlist(strsplit(entrez_string, "/"))
  mirna_list  <- mirna_df$miRNA_fin[match(entrez_list, mirna_df$target_entrez)]
  paste(unique(mirna_list[!is.na(mirna_list)]), collapse = "/")
}

pathways_wqs <- pathways_rev %>%
  mutate(miRNA = sapply(geneID, replace_entrez_with_mirna, mirna_df = input_comp_cl))

pathways_wqs <- separate(pathways_wqs, col = "ANALYSIS_2",
                         into = c("WQS_GROUP", "ASSOC_DIRECT"), sep = "_")
pathways_wqs$ASSOC_NUMB <- ifelse(pathways_wqs$ASSOC_DIRECT == "UP", 1, -1)

writexl::write_xlsx(pathways_wqs, "RESULT_ENRICHMENT_CLUSTERPROF_wMiRNA.xlsx")


# =============================================================================
# PART 2: PATHWAY GROUPING VIA SIMILARITY NETWORK (Louvain)
# =============================================================================

# --- 2.1 Pairwise pathway similarity -----------------------------------------
edox2             <- enrichplot::pairwise_termsim(formula_res_rev)
similarity_matrix <- as.matrix(data.frame(edox2@termsim))
colnames(similarity_matrix) <- row.names(similarity_matrix)

# Build adjacency matrix with similarity cutoff
similarity_cutoff <- 0.2
adjacency_matrix  <- (similarity_matrix >= similarity_cutoff) * 1
diag(adjacency_matrix) <- 0

pathway_graph <- graph_from_adjacency_matrix(adjacency_matrix,
                                             mode = "undirected", diag = FALSE)

# Assign edge weights from similarity matrix
edges      <- as_edgelist(pathway_graph)
edge_idx   <- cbind(edges[, 1], edges[, 2])
E(pathway_graph)$weight <- similarity_matrix[edge_idx]


# --- 2.2 Louvain community detection -----------------------------------------
communities   <- cluster_louvain(pathway_graph, weights = E(pathway_graph)$weight)
pathway_groups <- data.frame(
  Pathway = V(pathway_graph)$name,
  Group   = communities$membership
)

tiff("graph_pathways_grouping_wqs.tiff", res = 300, width = 35, height = 35, units = "in")
ggraph(pathway_graph, layout = "fr") +
  geom_edge_link(aes(width = weight), alpha = 0.6, color = "gray") +
  geom_node_point(aes(color = factor(communities$membership)), size = 5) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  theme_void() +
  ggtitle("Jaccard-Based Pathway Similarity Network (Louvain Clustering)") +
  labs(color = "Community", width = "Jaccard Similarity")
dev.off()

node_degree <- igraph::degree(pathway_graph, mode = "all")
pathway_info <- data.frame(
  Pathway   = V(pathway_graph)$name,
  Community = communities$membership,
  Degree    = node_degree
) %>% arrange(Community, desc(Degree))

write.csv(pathway_info, "pathway_communities_and_degree.csv", row.names = FALSE)

hub_pathways <- pathway_info %>%
  group_by(Community) %>%
  slice_max(order_by = Degree, n = 1, with_ties = FALSE) %>%
  select(Community, Hub_Pathway = Pathway)

pathway_info <- pathway_info %>%
  left_join(hub_pathways, by = "Community")


# --- 2.3 Reactome hierarchy grouping -----------------------------------------
ptw_hr  <- read.delim("./Reactome_pathways_March2026.txt",    header = FALSE)
all_ptw <- read.delim("./ReactomePathways_list_March2026.txt", header = FALSE)

ptw_2 <- merge(pathway_info,
               unique(pathways_wqs[, c("Description", "ID")]),
               by.x = "Pathway", by.y = "Description", all.x = TRUE)
ptw_2 <- ptw_2[!duplicated(ptw_2), ]

ptw_hr_1 <- ptw_hr[ptw_hr$V1 %in% ptw_2$ID | ptw_hr$V2 %in% ptw_2$ID, ]

grp         <- graph_from_data_frame(ptw_hr_1, directed = FALSE)
comm        <- cluster_louvain(grp, weights = NULL)
ptw_degree  <- igraph::degree(grp, mode = "all")

ptw_hr_2 <- data.frame(
  Pathway   = V(grp)$name,
  Group_2   = comm$membership,
  Degree    = ptw_degree
)

ptw_hr_2$HIERARCHY <- ifelse(
  ptw_hr_2$Pathway %in% ptw_hr_1$V1 & !(ptw_hr_2$Pathway %in% ptw_hr_1$V2),
  "HUB", "NO"
)

ptw_hr_2 <- ptw_hr_2 %>%
  arrange(Group_2, desc(Degree)) %>%
  group_by(Group_2) %>%
  mutate(Hub_Value = if ("HUB" %in% HIERARCHY) {
    Pathway[HIERARCHY == "HUB"][1]
  } else {
    Pathway[which.max(Degree)]
  }) %>%
  ungroup()

ptw_hr_3 <- merge(ptw_hr_2, ptw_2, by.x = "Pathway", by.y = "ID")
ptw_hr_3 <- merge(ptw_hr_3, all_ptw, by.x = "Hub_Value", by.y = "V1", all.x = TRUE)
ptw_hr_3 <- ptw_hr_3 %>% arrange(Group_2)


# =============================================================================
# PART 3: BINARY HEATMAP — pathway group x BWQS cluster
# =============================================================================

# --- 3.1 Merge enrichment with pathway groups --------------------------------
pathways_wqs_sb <- pathways_wqs[pathways_wqs$Count >= 4, ]

pathway_wqs_macro <- merge(pathways_wqs_sb, ptw_hr_3[, c("Pathway", "Hub_Value", "V2")],
                           by.x = "ID", by.y = "Pathway", all.x = TRUE)

# Keep one pathway per Cluster x direction x group (lowest q-value)
pathway_wqs_macro_sb <- pathway_wqs_macro %>%
  arrange(qvalue) %>%
  filter(!duplicated(cbind(Cluster, ASSOC_NUMB, V2)))


# --- 3.2 Wide binary matrix --------------------------------------------------
binary_enr <- pathway_wqs_macro_sb %>%
  dplyr::select(Cluster, ASSOC_NUMB, V2) %>%
  pivot_wider(names_from = Cluster, values_from = ASSOC_NUMB, values_fill = 0) %>%
  as.data.frame()

row.names(binary_enr) <- binary_enr$V2
binary_enr$V2         <- NULL


# --- 3.3 Split into maternal and cordonal ------------------------------------
binary_M <- binary_enr[, grepl("MATERNAL", colnames(binary_enr)), drop = FALSE]
binary_C <- binary_enr[, grepl("CORD",     colnames(binary_enr)), drop = FALSE]

# Enforce column order
cols_M <- c("MATERNAL.EEs_DOWN", "MATERNAL.EEs_UP",
            "MATERNAL.OCs_DOWN", "MATERNAL.OCs_UP",
            "MATERNAL.TOT_DOWN", "MATERNAL.TOT_UP")
cols_C <- c("CORD.EEs_DOWN", "CORD.EEs_UP",
            "CORD.OCs_DOWN", "CORD.OCs_UP",
            "CORD.TOT_DOWN", "CORD.TOT_UP")

binary_M <- binary_M[, intersect(cols_M, colnames(binary_M)), drop = FALSE]
binary_C <- binary_C[, intersect(cols_C, colnames(binary_C)), drop = FALSE]


# --- 3.4 Row and column annotations ------------------------------------------
col_split   <- rep(c("EEs", "POPs", "TOT"), each = 2)
col_split_C <- col_split

# Row grouping from Reactome hierarchy community
split_Mr_2 <- ptw_hr_3$Group_2[match(row.names(binary_M), ptw_hr_3$V2)]
split_Cr_2 <- ptw_hr_3$Group_2[match(row.names(binary_C), ptw_hr_3$V2)]

col_fun <- colorRamp2(c(1, 0, -1), c("red", "white", "blue"))

ha_col <- HeatmapAnnotation(
  ANALYSIS = anno_block(
    gp     = gpar(fill = "white"),
    labels = c("EEs", "POPs", "TOT"),
    labels_gp = gpar(col = "black", fontsize = 35)
  )
)


# --- 3.5 Draw heatmaps -------------------------------------------------------
ph_M <- Heatmap(
  as.matrix(binary_M),
  col                  = col_fun,
  cluster_rows         = FALSE,
  cluster_columns      = FALSE,
  cluster_column_slices = FALSE,
  column_title         = " ",
  column_title_gp      = gpar(fontsize = 30),
  show_column_names    = FALSE,
  column_split         = col_split,
  row_split            = split_Mr_2,
  show_row_names       = FALSE,
  row_names_gp         = gpar(fontsize = 30),
  border               = TRUE,
  row_title            = NULL,
  bottom_annotation    = ha_col,
  rect_gp              = gpar(col = "gray", lwd = 1),
  row_names_max_width  = max_text_width(rownames(binary_M), gp = gpar(fontsize = 20)),
  show_heatmap_legend  = FALSE
)

ph_C <- Heatmap(
  as.matrix(binary_C),
  col                  = col_fun,
  cluster_rows         = FALSE,
  cluster_columns      = FALSE,
  cluster_column_slices = FALSE,
  column_title         = " ",
  column_title_gp      = gpar(fontsize = 30),
  show_column_names    = FALSE,
  column_split         = col_split_C,
  row_split            = split_Cr_2,
  show_row_names       = TRUE,
  row_names_gp         = gpar(fontsize = 30),
  border               = TRUE,
  row_title            = NULL,
  bottom_annotation    = ha_col,
  rect_gp              = gpar(col = "gray", lwd = 1),
  row_names_max_width  = max_text_width(rownames(binary_C), gp = gpar(fontsize = 35)),
  show_heatmap_legend  = FALSE
)


# --- 3.6 Combine and save ----------------------------------------------------
tiff("Figure7_enrichment.tiff", res = 300, width = 40, height = 35, units = "in")

draw(ph_M + ph_C,
     padding              = unit(c(2, 2, 2, 2), "cm"),
     heatmap_legend_side  = "left")

seekViewport("global")

# MATERNAL label box
grid.rect(x = unit(0.133, "npc"), y = unit(0.97, "npc"),
          width = unit(0.215, "npc"), height = unit(0.03, "npc"),
          gp = gpar(fill = "white", col = "black", lwd = 1.5), just = "top")
grid.text("MATERNAL", x = unit(0.133, "npc"), y = unit(0.955, "npc"),
          gp = gpar(fontsize = 35, fontface = "bold", col = "black"))

# CORD label box
grid.rect(x = unit(0.35, "npc"), y = unit(0.97, "npc"),
          width = unit(0.215, "npc"), height = unit(0.03, "npc"),
          gp = gpar(fill = "white", col = "black", lwd = 1.5), just = "top")
grid.text("CORD", x = unit(0.35, "npc"), y = unit(0.955, "npc"),
          gp = gpar(fontsize = 35, fontface = "bold", col = "black"))

dev.off()

cat("Figure saved to Figure7_enrichment.tiff\n")
