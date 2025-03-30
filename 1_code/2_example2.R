# ===================================================
# GENE ONTOLOGY (GO) TERM SIMILARITY NETWORK ANALYSIS
# This script analyzes and visualizes relationships between
# enriched GO terms from pathway analysis results
# ===================================================

# Load required packages
library(r4projects)      # For project management functions
library(simplifyEnrichment)  # For enrichment analysis simplification
library(tidyverse)       # Collection of data manipulation packages
library(tibble)          # Enhanced data frames
library(igraph)          # Network analysis package
library(readr)           # For reading csv files
library(tidygraph)       # Tidy interface to graph manipulation
library(ggsci)           # Scientific journal color palettes
library(ggraph)          # Grammar of graphics for networks
library(shadowtext)      # For text with shadows/outlines
library(extrafont)       # For additional fonts


# Set working directory to the root of the project
# Using r4projects functionality to ensure consistent path handling
setwd(r4projects::get_project_wd())
rm(list = ls())

####################################################
# DATA IMPORT AND EXPLORATION
####################################################

# Read pathway analysis results from GitHub repository
# Contains Gene Ontology (GO) enrichment analysis results from different ontologies:
# BP (Biological Process), MF (Molecular Function), CC (Cellular Component)
pathway_result <-
  readr::read_csv(
    "https://raw.githubusercontent.com/jaspershen-lab/MD7115-NTU/refs/heads/main/2_demo_data/pathway_result.csv"
  )

# Display the first few rows of the pathway results to examine the structure
head(pathway_result)

####################################################
# COMPUTE GO TERM SIMILARITY MATRICES
####################################################

# Process Biological Process (BP) GO terms
# Calculate semantic similarity between BP GO terms using Wang's method
# Wang's method considers the topology of the GO graph structure and
# semantic contribution of each GO term to the similarity calculation
bp_sim_matrix <-
  simplifyEnrichment::GO_similarity(go_id = pathway_result$ID[pathway_result$ONTOLOGY == "BP"],
                                    ont = "BP",
                                    measure = "Wang") %>%
  as.data.frame() %>%
  # Convert row names to a column for easier manipulation
  tibble::rownames_to_column(var = "name1") %>%
  # Convert from wide to long format for pairwise comparison
  tidyr::pivot_longer(cols = -name1,
                      names_to = "name2",
                      values_to = "sim") %>%
  # Remove self-comparisons (diagonal of similarity matrix)
  dplyr::filter(name1 != name2) %>%
  # Filter for meaningful similarities (above 0.5 threshold)
  # Higher values indicate greater semantic similarity
  dplyr::filter(sim > 0.5)

# Create a unique identifier for each GO term pair
# This helps remove duplicate pairs (A-B and B-A are the same relationship)
# by sorting the term names and creating a consistent identifier
name <- apply(bp_sim_matrix, 1, function(x) {
  paste(sort(x[1:2]), collapse = "_")
})

# Remove duplicate GO term pairs while keeping the data structure
# We only need one entry for each pair of GO terms
bp_sim_matrix <-
  bp_sim_matrix %>%
  dplyr::mutate(name = name) %>%
  dplyr::arrange(name) %>%
  dplyr::distinct(name, .keep_all = TRUE) %>%
  dplyr::select(-name)

# Process Molecular Function (MF) GO terms
# Same approach as with BP terms, but using MF ontology
mf_sim_matrix <-
  simplifyEnrichment::GO_similarity(go_id = pathway_result$ID[pathway_result$ONTOLOGY == "MF"],
                                    ont = "MF",
                                    measure = "Wang") %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "name1") %>%
  tidyr::pivot_longer(cols = -name1,
                      names_to = "name2",
                      values_to = "sim") %>%
  dplyr::filter(name1 != name2) %>%
  dplyr::filter(sim > 0.5)

# Create unique identifiers for MF term pairs
name <- apply(mf_sim_matrix, 1, function(x) {
  paste(sort(x[1:2]), collapse = "_")
})

# Remove duplicate MF term pairs
mf_sim_matrix <-
  mf_sim_matrix %>%
  dplyr::mutate(name = name) %>%
  dplyr::arrange(name) %>%
  dplyr::distinct(name, .keep_all = TRUE) %>%
  dplyr::select(-name)

# Process Cellular Component (CC) GO terms
# Same approach as with BP and MF terms, but using CC ontology
cc_sim_matrix <-
  simplifyEnrichment::GO_similarity(go_id = pathway_result$ID[pathway_result$ONTOLOGY == "CC"],
                                    ont = "CC",
                                    measure = "Wang") %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "name1") %>%
  tidyr::pivot_longer(cols = -name1,
                      names_to = "name2",
                      values_to = "sim") %>%
  dplyr::filter(name1 != name2) %>%
  dplyr::filter(sim > 0.5)

# Create unique identifiers for CC term pairs
name <- apply(cc_sim_matrix, 1, function(x) {
  paste(sort(x[1:2]), collapse = "_")
})

# Remove duplicate CC term pairs
cc_sim_matrix <-
  cc_sim_matrix %>%
  dplyr::mutate(name = name) %>%
  dplyr::arrange(name) %>%
  dplyr::distinct(name, .keep_all = TRUE) %>%
  dplyr::select(-name)

# Combine all three ontology similarity matrices into one
# This allows us to analyze relationships across all GO domains
sim_matrix <-
  rbind(bp_sim_matrix, mf_sim_matrix, cc_sim_matrix) %>%
  tibble::as_tibble()

# Preview the combined similarity matrix
head(sim_matrix)

####################################################
# PREPARE GRAPH DATA FOR NETWORK VISUALIZATION
####################################################

# Create edge data from similarity matrix
# Edges represent semantic relationships between GO terms
edge_data <-
  sim_matrix %>%
  # Rename columns to standard network terminology (from/to)
  dplyr::rename(from = name1, to = name2) %>%
  # Keep only strong relationships (similarity > 0.5)
  dplyr::filter(sim > 0.5)

# Display sample of edge data
head(edge_data)

# Create node data from pathway results
# Nodes represent individual GO terms
node_data <-
  pathway_result %>%
  dplyr::select(ID, everything()) %>%
  dplyr::rename(node = ID)  # Rename ID column to 'node' for network compatibility

# Display sample of node data
head(node_data)

# Create a graph data structure using tidygraph
# Combines node and edge information into a single object
graph_data <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = FALSE) %>%  # Undirected graph as similarity is symmetric
  # Calculate node degree (number of connections) as a centrality measure
  # Higher degree indicates GO terms that share similarity with many other terms
  dplyr::mutate(degree = tidygraph::centrality_degree())

####################################################
# COMMUNITY DETECTION AND MODULE ASSIGNMENT
####################################################

# Detect subnetworks/communities within the graph
# Edge betweenness clustering algorithm identifies modular structures
# This helps group related GO terms into functional modules
subnetwork <-
  igraph::cluster_edge_betweenness(graph = graph_data, weights = abs(edge_attr(graph_data, "sim")))

# Create module labels for each detected community
# Each GO term will be assigned to a specific module
cluster <-
  paste("GO_Module", as.character(igraph::membership(subnetwork)), sep = "_")

# Add module information to the graph data
graph_data <-
  graph_data %>%
  tidygraph::mutate(module = cluster)

# Extract node attributes with module assignments
# This creates a table of GO terms with their module assignments
result_with_module <-
  igraph::vertex_attr(graph_data) %>%
  do.call(cbind, .) %>%
  tibble::as_tibble() %>%
  # Convert p.adjust to numeric for proper sorting
  dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
  # Sort by module and significance (most significant first)
  dplyr::arrange(module, p.adjust)

# Further organize results by ontology type, module, and significance
# This helps identify patterns within each ontology category
result_with_module <-
  result_with_module %>%
  dplyr::arrange(ONTOLOGY, module, p.adjust)

# Calculate the number of GO terms in each module
# This helps assess module size and importance
module_content_number <-
  result_with_module %>%
  dplyr::count(module) %>%
  dplyr::rename(module_content_number = n)

# Add module size information to results
result_with_module <-
  result_with_module %>%
  dplyr::left_join(module_content_number, by = "module")

# Add module size information to graph data
# This will be used for visualization
graph_data <-
  graph_data %>%
  # Activate nodes for manipulation
  activate(what = "nodes") %>%
  dplyr::left_join(module_content_number, by = "module")

# Display the final graph data structure
graph_data

####################################################
# VISUALIZATION OF GO TERM CLUSTERS
####################################################

# Select representative labels for each module
# Choose the most significant GO term with highest count in each module
# This identifies the most important term to label each module in the visualization
cluster_label_module <-
  igraph::as_data_frame(graph_data, what = "vertices") %>%
  dplyr::group_by(module) %>%
  dplyr::filter(p.adjust == min(p.adjust) &
                  Count == max(Count)) %>%
  dplyr::slice_head(n = 1) %>%  # Take the first if there are ties
  pull(Description)  # Extract just the description text

# Get all GO term descriptions for reference
cluster_label_all <-
  igraph::as_data_frame(graph_data, what = "vertices")$Description

# Create network visualization
# Initialize the plot with a force-directed layout (Fruchterman-Reingold)
plot =
  graph_data %>%
  ggraph(layout = 'fr', circular = FALSE) +
  # Add edges with width based on similarity strength
  geom_edge_link(
    aes(width = sim),
    # Edge width represents similarity strength
    strength = 1,
    color = "black",
    alpha = 1,
    show.legend = TRUE
  ) +
  # Add nodes with size based on significance and color by module
  geom_node_point(
    aes(fill = module, # Color nodes by module assignment
        size = -log(p.adjust, 10)),
    # Size nodes by significance (-log10 of adjusted p-value)
    shape = 21,
    # Circle with border
    alpha = 1,
    show.legend = TRUE
  )

# Add text labels for the representative terms in each module
plot =
  plot +
  geom_node_text(aes(
    x = x,
    y = y,
    # Only label the representative terms for each module
    label = ifelse(Description %in% cluster_label_module, Description, NA)
  ),
  size = 3,
  repel = TRUE)  # Repel text to avoid overlap

# Configure plot scales and guides
plot =
  plot +
  guides(fill = guide_legend(ncol = 1)) +  # Single column for module colors
  scale_edge_width_continuous(range = c(0.1, 2)) +  # Edge width range
  scale_size_continuous(range = c(1, 7))  # Node size range

# Set theme and appearance
plot =
  plot +
  ggraph::theme_graph() +  # Use the graph theme
  theme(
    # Make background transparent for better export
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  )

# Display the plot
plot

# Load font libraries for PDF output
extrafont::loadfonts()

# Save the plot as a PDF file
ggplot2::ggsave(plot,
                filename = "similarity_network_plot.pdf",
                width = 9,
                height = 7)