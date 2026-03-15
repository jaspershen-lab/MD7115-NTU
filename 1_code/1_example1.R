# ===================================================
# NETWORK VISUALIZATION FOR MULTI-OMICS AND WEARABLE DATA
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

# Set working directory to project root
# This helps maintain consistent file paths across different systems
setwd(r4projects::get_project_wd())
rm(list = ls())

# ===================================================
# DATA IMPORT
# ===================================================

# Import node data from GitHub repository
# This contains information about each node in our network
node_data <-
  readr::read_csv(
    "https://raw.githubusercontent.com/jaspershen-lab/MD7115-NTU/main/2_demo_data/example_node_data.csv"
  )

# Import edge data from GitHub repository
# This contains information about connections between nodes
edge_data <-
  readr::read_csv(
    "https://raw.githubusercontent.com/jaspershen-lab/MD7115-NTU/main/2_demo_data/example_edge_data.csv"
  )

# Display the first few rows of each dataset to inspect
head(node_data)
head(edge_data)

# ===================================================
# NETWORK CREATION
# ===================================================

# Create a tidygraph object from our nodes and edges
# This creates an undirected graph (directed = FALSE)
# Also calculate the degree centrality for each node
total_graph <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = FALSE) %>%
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

# ===================================================
# NETWORK LAYOUT PREPARATION
# ===================================================

# Create a bipartite graph layout (nodes separated into two types)
g <- total_graph

# Set node types for bipartite layout
V(g)$type <- bipartite_mapping(g)$type

# Create initial bipartite layout and extract coordinates
coords <-
  create_layout(g, layout = "bipartite") %>%
  dplyr::select(node, class, mol_name, class2, x, y)

# Adjust y-coordinates based on data type
# Initially setting wearable and internal-omics to different levels
coords$y[coords$class2 == "wearable"] = 0.3
coords$y[coords$class2 == "internal-omics"] = 1

# Check distribution of node classes
table(coords$class)

# Further refine y-coordinates to separate different omics types
coords$y[coords$class2 == "wearable"] = -0.5    # Wearable data at bottom
coords$y[coords$class == "metabolomics"] = 1    # Metabolomics at top
coords$y[coords$class == "lipidomics"] = 0.5    # Lipidomics in middle
# Several omics types at same level
coords$y[coords$class == "proteomics"] = 0
coords$y[coords$class == "cytokine"] = 0
coords$y[coords$class == "metabolic_panel"] = 0
coords$y[coords$class == "total_protein"] = 0

# Convert layout to circular/radial coordinates
# This transforms the bipartite layout into a circular layout
coords <-
  coords %>%
  dplyr::select(x, y) %>%
  dplyr::mutate(
    theta = x / (max(x) + 1) * 2 * pi,  # Convert x to angle (theta)
    r = y + 1,                          # Convert y to radius
    x = r * cos(theta),                 # Calculate new x position
    y = r * sin(theta)                  # Calculate new y position
  )

# Create new graph layout with manual positioning
my_graph <-
  create_layout(
    graph = g,
    layout = "manual",
    x = coords$x,
    y = coords$y
  )

# ===================================================
# COLOR PALETTE DEFINITION
# ===================================================

# Define color palette for internal omics data classes
# Using AAAS journal color palette
class_color =
  c(
    "lipidomics" = ggsci::pal_aaas()(10)[1],
    "metabolomics" = ggsci::pal_aaas()(10)[3],
    "cytokine" = ggsci::pal_aaas()(10)[4],
    "total_protein" = ggsci::pal_aaas()(10)[5],
    "cortisol" = ggsci::pal_aaas()(10)[6],
    "metabolic_panel" = ggsci::pal_aaas()(10)[7],
    "proteomics" = ggsci::pal_aaas()(10)[8]
  )

# Define color palette for wearable data types
# Using D3.js color palette
wearable_color =
  c(
    "sleep" = ggsci::pal_d3()(n = 10)[2],
    "cgm"  = ggsci::pal_d3()(n = 10)[6],      # Continuous glucose monitoring
    "hr" = ggsci::pal_d3()(n = 10)[7],        # Heart rate
    "step"  = ggsci::pal_d3()(n = 10)[9],     # Step count
    "food" = ggsci::pal_d3()(n = 10)[10]      # Food intake
  )

# ===================================================
# NETWORK VISUALIZATION
# ===================================================

# Create the network visualization using ggraph
plot <-
  ggraph(my_graph, layout = 'bipartite') +
  # Add edges as diagonal lines with color based on correlation
  # and width based on adjusted p-value significance
  geom_edge_diagonal(aes(
    color = lagged_cor,                      # Color by correlation value
    width = -log(lagged_cor_p_adjust, 10)    # Width by significance (-log10 p-value)
  ),
  alpha = 0.5,
  show.legend = TRUE) +
  
  # Add nodes as points with fill color by class and size by degree
  geom_node_point(
    aes(fill = class, size = Degree),
    shape = 21,                             # Circle with border
    alpha = 0.5,
    show.legend = TRUE
  ) +
  
  # Add text labels for wearable nodes with shadow
  shadowtext::geom_shadowtext(
    aes(
      x = x,
      y = y,
      label = ifelse(class2 == "wearable", mol_name, NA),  # Only label wearable nodes
      color = class
    ),
    bg.color = "white",                     # White background for visibility
    size = 5,
    show.legend = FALSE
  ) +
  
  # Add text labels for omics nodes with shadow
  # with angle adjustment based on position in the circle
  shadowtext::geom_shadowtext(
    aes(
      x = x,
      y = y,
      label = ifelse(class2 == "wearable", NA, mol_name),  # Only label omics nodes
      color = class,
      nudge_y = ifelse(class2 == "wearable", "inward", 'outward'),
      angle = -((-node_angle(x, y) + 90) %% 180) + 90      # Adjust text angle based on position
    ),
    bg.color = "white",
    size = 2,
    show.legend = FALSE,
    check_overlap = TRUE                    # Avoid text overlap
  ) +
  
  # Set scales for visual elements
  scale_size_continuous(range = c(3, 10)) +                      # Node size range
  scale_fill_manual(values = c(class_color, wearable_color)) +   # Node fill colors
  scale_color_manual(values = c(class_color, wearable_color)) +  # Text colors
  
  # Configure legend appearance
  guides(
    linetype = "none",
    color = guide_colorbar(title = "Lagged correlation", override.aes = list(linetype = "none")),
    size = guide_legend(
      title = "Degree",
      override.aes = list(
        linetype = NA,
        fill = "transparent",
        shape = 21,
        color = "black"
      )
    ),
    fill = guide_legend(
      title = "Class",
      override.aes = list(
        shape = 21,
        size = 3,
        alpha = 1
      )
    )
  ) +
  
  # Set edge width range
  scale_edge_width_continuous(range = c(0.3, 2)) +
  
  # Set edge color gradient (blue for negative, red for positive correlations)
  ggraph::scale_edge_color_gradientn(colours = c(alpha("#3B4992FF", 0.7), "white", alpha("#EE0000FF", 0.7))) +
  
  # Use ggraph's graph theme
  ggraph::theme_graph() +
  
  # Customize theme for transparent background
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  )

# Display the plot
plot

# Load additional fonts for better text rendering
extrafont::loadfonts()

# Save the plot as a PDF file
ggsave(plot,
       filename = "example_network.pdf",
       width = 9.1,
       height = 7)
