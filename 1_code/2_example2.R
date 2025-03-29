####################################################
# SETUP AND INSTALLATION OF REQUIRED PACKAGES
####################################################

# Install and load remotes package if not already installed
# remotes is used for installing packages from GitHub
if (!require(remotes)) {
  install.packages("remotes")
}

# Install and load r4projects package from GitHub if not available
# This package helps with project workflow management
if (!require(r4projects)) {
  remotes::install_github("jaspershen/r4projects")
}

# Install required packages if not available
# simplifyEnrichment is a Bioconductor package for simplifying enrichment results
if (!require(simplifyEnrichment)) {
  BiocManager::install("simplifyEnrichment")
}

# tibble provides enhanced data frame functionality
if (!require(tibble)) {
  install.packages("tibble")
}


# Set working directory to the root of the project
# Using r4projects functionality to ensure consistent path handling
setwd(r4projects::get_project_wd())

####################################################
# DATA IMPORT AND EXPLORATION
####################################################

# Read pathway analysis results from GitHub repository
# Contains Gene Ontology (GO) enrichment analysis results
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
# Calculate similarity between BP GO terms using Wang's method
bp_sim_matrix <-
  simplifyEnrichment::GO_similarity(go_id = pathway_result$ID[pathway_result$ONTOLOGY == "BP"],
                                    ont = "BP",
                                    measure = "Wang") %>%
  as.data.frame() %>%
  # Convert row names to a column for easier manipulation
  tibble::rownames_to_column(var = "name1") %>%
  # Convert from wide to long format
  tidyr::pivot_longer(cols = -name1,
                      names_to = "name2",
                      values_to = "sim") %>%
  # Remove self-comparisons (diagonal of similarity matrix)
  dplyr::filter(name1 != name2) %>%
  # Filter for meaningful similarities (above 0.5 threshold)
  dplyr::filter(sim > 0.5)

# Create a unique identifier for each GO term pair
# This helps remove duplicate pairs (A-B and B-A are the same relationship)
name <- apply(bp_sim_matrix, 1, function(x) {
  paste(sort(x[1:2]), collapse = "_")
})

# Remove duplicate GO term pairs while keeping the data structure
bp_sim_matrix <-
  bp_sim_matrix %>%
  dplyr::mutate(name = name) %>%
  dplyr::arrange(name) %>%
  dplyr::distinct(name, .keep_all = TRUE) %>%
  dplyr::select(-name)

# Process Molecular Function (MF) GO terms
# Similar approach as with BP terms
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

name <- apply(mf_sim_matrix, 1, function(x) {
  paste(sort(x[1:2]), collapse = "_")
})

mf_sim_matrix <-
  mf_sim_matrix %>%
  dplyr::mutate(name = name) %>%
  dplyr::arrange(name) %>%
  dplyr::distinct(name, .keep_all = TRUE) %>%
  dplyr::select(-name)

# Process Cellular Component (CC) GO terms
# Similar approach as with BP and MF terms
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

name <- apply(cc_sim_matrix, 1, function(x) {
  paste(sort(x[1:2]), collapse = "_")
})

cc_sim_matrix <-
  cc_sim_matrix %>%
  dplyr::mutate(name = name) %>%
  dplyr::arrange(name) %>%
  dplyr::distinct(name, .keep_all = TRUE) %>%
  dplyr::select(-name)

# Combine all three ontology similarity matrices into one
sim_matrix <-
  rbind(bp_sim_matrix, mf_sim_matrix, cc_sim_matrix) %>%
  tibble::as_tibble()

# Preview the combined similarity matrix
head(sim_matrix)

####################################################
# PREPARE GRAPH DATA FOR NETWORK VISUALIZATION
####################################################

# Create edge data from similarity matrix
# Edges represent relationships between GO terms
edge_data <-
  sim_matrix %>%
  # Rename columns to standard network terminology
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
  dplyr::rename(node = ID)

# Display sample of node data
head(node_data)

# Create a graph data structure using tidygraph
# Combines node and edge information into a single object
graph_data <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = FALSE) %>%
  # Calculate node degree (number of connections) as a centrality measure
  dplyr::mutate(degree = tidygraph::centrality_degree())

####################################################
# COMMUNITY DETECTION AND MODULE ASSIGNMENT
####################################################

# Detect subnetworks/communities within the graph
# Edge betweenness clustering algorithm identifies modular structures
subnetwork <-
  igraph::cluster_edge_betweenness(graph = graph_data, weights = abs(edge_attr(graph_data, "sim")))

# Create module labels for each detected community
cluster <-
  paste("GO_Module", as.character(igraph::membership(subnetwork)), sep = "_")

# Add module information to the graph data
graph_data <-
  graph_data %>%
  tidygraph::mutate(module = cluster)

# Extract node attributes with module assignments
result_with_module <-
  igraph::vertex_attr(graph_data) %>%
  do.call(cbind, .) %>%
  tibble::as_tibble() %>%
  # Convert p.adjust to numeric for proper sorting
  dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
  # Sort by module and significance
  dplyr::arrange(module, p.adjust)

# Further organize results by ontology type, module, and significance
result_with_module <-
  result_with_module %>%
  dplyr::arrange(ONTOLOGY, module, p.adjust)

# Calculate the number of GO terms in each module
module_content_number <-
  result_with_module %>%
  dplyr::count(module) %>%
  dplyr::rename(module_content_number = n)

# Add module size information to results
result_with_module <-
  result_with_module %>%
  dplyr::left_join(module_content_number, by = "module")

# Add module size information to graph data
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
cluster_label_module <-
  igraph::as_data_frame(graph_data, what = "vertices") %>%
  dplyr::group_by(module) %>%
  dplyr::filter(p.adjust == min(p.adjust) &
                  Count == max(Count)) %>%
  dplyr::slice_head(n = 1) %>%
  pull(Description)

# Get all GO term descriptions for reference
cluster_label_all <-
  igraph::as_data_frame(graph_data, what = "vertices")$Description

# Create network visualization
plot =
  graph_data %>%
  ggraph(layout = 'fr', circular = FALSE) +
  geom_edge_link(
    aes(width = sim),
    strength = 1,
    color = "black",
    alpha = 1,
    show.legend = TRUE
  ) +
  geom_node_point(
    aes(fill = module, size = -log(p.adjust, 10)),
    shape = 21,
    alpha = 1,
    show.legend = TRUE
  )

plot =
  plot +
  geom_node_text(aes(
    x = x,
    y = y,
    label = ifelse(Description %in% cluster_label_module, Description, NA)
  ),
  size = 3,
  repel = TRUE)

plot =
  plot +
  guides(fill = guide_legend(ncol = 1)) +
  scale_edge_width_continuous(range = c(0.1, 2)) +
  scale_size_continuous(range = c(1, 7))

plot =
  plot +
  ggraph::theme_graph() +
  theme(
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
