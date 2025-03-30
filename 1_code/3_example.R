# ===================================================
# MICROBIOME TAXONOMIC TREE VISUALIZATION WITH PERSONALIZED SCORES
# This script creates a comprehensive circular phylogenetic tree 
# visualization of microbiome data with body site-specific enrichment scores
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
library(tidyverse)       # Data manipulation and visualization
library(ggtreeExtra)     # Extensions for ggtree
library(ggtree)          # Visualization of phylogenetic trees
library(microbiomeViz)   # Visualization for microbiome data
library(phyloseq)        # Analysis of microbiome census data

# Set working directory to the root of the project
setwd(r4projects::get_project_wd())

# Clear the workspace
rm(list = ls())

# ===================================================
# DATA IMPORT
# ===================================================

# Import personalized microbiome scores from GitHub repository
# This contains enrichment values (fc1, fc2) and significance metrics for genera across body sites
personalized_score =
  readr::read_csv(
    "https://raw.githubusercontent.com/jaspershen-lab/MD7115-NTU/refs/heads/main/2_demo_data/personalized_score.csv"
  )

# Display the first few rows to examine the structure
head(personalized_score)

# Load pre-processed microbiome data (phyloseq object)
# This contains the taxonomic hierarchy and abundance information
load(
  url(
    "https://github.com/jaspershen-lab/MD7115-NTU/raw/main/2_demo_data/microbiome_data.RData"
  )
)

# View the microbiome data structure
microbiome_data

# ===================================================
# TREE CREATION AND BASIC VISUALIZATION
# ===================================================

# Create a taxonomic tree from the phyloseq object
# use_abundance=FALSE means we're only using taxonomy, not abundance values
tree =
  microbiomeViz::parsePhyloseq(
    physeq = microbiome_data,
    use_abundance = FALSE,
    node.size.scale = 0,   # Don't scale node size based on abundance
    node.size.offset = 0
  )

# Create a basic circular tree visualization
raw_plot <-
  tree.backbone(
    tree = tree,
    size = 0.3,         # Line size
    shape = 16,         # Node shape (solid circle)
    layout = "circular", # Circular tree layout
    fill = "black",
    color = "black"
  )

# Display the basic tree
raw_plot

# ===================================================
# TREE NODE ANNOTATION AND ENHANCEMENT
# ===================================================

# Convert abbreviated taxonomic level codes to full names
# This makes the visualization more interpretable
raw_plot$data$nodeClass2 =
  as.character(raw_plot$data$nodeClass)

# Assign taxonomic levels to each node
raw_plot$data$nodeClass2[is.na(raw_plot$data$nodeClass2)] = "Root"
raw_plot$data$nodeClass2[raw_plot$data$nodeClass2 == "f"] = "Family"
raw_plot$data$nodeClass2[raw_plot$data$nodeClass2 == "c"] = "Class"
raw_plot$data$nodeClass2[raw_plot$data$nodeClass2 == "o"] = "Order"
raw_plot$data$nodeClass2[raw_plot$data$nodeClass2 == "p"] = "Phylum"
raw_plot$data$nodeClass2[raw_plot$data$nodeClass2 == "g"] = "Genus"
raw_plot$data$nodeClass2[raw_plot$data$nodeClass2 == "k"] = "Kingdom"

# Set node sizes based on taxonomic level
# Higher taxonomic levels (e.g., Kingdom) get larger nodes than lower levels (e.g., Genus)
raw_plot$data$nodeSize2 = raw_plot$data$nodeSize
raw_plot$data$nodeSize2[raw_plot$data$nodeClass2 == "Root"] = 3.5
raw_plot$data$nodeSize2[raw_plot$data$nodeClass2 == "Kingdom"] = 3
raw_plot$data$nodeSize2[raw_plot$data$nodeClass2 == "Phylum"] = 2.5
raw_plot$data$nodeSize2[raw_plot$data$nodeClass2 == "Class"] = 2
raw_plot$data$nodeSize2[raw_plot$data$nodeClass2 == "Order"] = 1.5
raw_plot$data$nodeSize2[raw_plot$data$nodeClass2 == "Family"] = 1
raw_plot$data$nodeSize2[raw_plot$data$nodeClass2 == "Genus"] = 0.5

# Add colored nodes to the tree based on taxonomic level
# Each taxonomic level will have a different color and size
raw_plot =
  raw_plot +
  ggtree::geom_point2(aes(color = nodeClass2, size = nodeSize2), show.legend = FALSE) +
  ggnewscale::new_scale(new_aes = "fill")  # Prepare for next layer with new fill scale

# Display the enhanced tree
raw_plot

# Extract genus names from the node labels
# Split the labels at "__" and take the second part (the actual taxon name)
raw_plot$data$label2 =
  raw_plot$data$label %>%
  stringr::str_split("__") %>%
  purrr::map(function(x) {
    x[2]  # Extract the second element (taxon name)
  }) %>%
  unlist()

# Only keep genus names for genus-level nodes, set others to NA
raw_plot$data$label2[as.character(raw_plot$data$nodeClass) != "g"] = NA

# ===================================================
# DATA PREPARATION FOR ENRICHMENT VISUALIZATION
# ===================================================

# Reshape personalized score data for first fold change (fc1)
# Convert from long to wide format with genera as rows and body sites as columns
personalized_score_fc1 =
  personalized_score %>%
  dplyr::select(genus, class, fc1) %>%
  tidyr::pivot_wider(names_from = class, values_from = "fc1")

# Add suffix to column names to indicate they contain fc1 values
colnames(personalized_score_fc1)[-1] = paste(colnames(personalized_score_fc1)[-1], "fc1", sep = "_")

# Reshape data for second fold change (fc2)
personalized_score_fc2 =
  personalized_score %>%
  dplyr::select(genus, class, fc2) %>%
  tidyr::pivot_wider(names_from = class, values_from = "fc2")

# Add suffix to column names to indicate they contain fc2 values
colnames(personalized_score_fc2)[-1] = paste(colnames(personalized_score_fc2)[-1], "fc2", sep = "_")

# Reshape data for adjusted p-values of fc1
personalized_score_fc1_p_adjust =
  personalized_score %>%
  dplyr::select(genus, class, fc1_p_adjust) %>%
  tidyr::pivot_wider(names_from = class, values_from = "fc1_p_adjust")

# Add suffix to column names to indicate they contain adjusted p-values
colnames(personalized_score_fc1_p_adjust)[-1] =
  paste(colnames(personalized_score_fc1_p_adjust)[-1],
        "fc1_p_adjust",
        sep = "_")

# Extract taxonomic information from the phyloseq object
variable_info <-
  tax_table(microbiome_data)

# Combine all data into a comprehensive information table
new_info =
  data.frame(Genus = raw_plot$data$label2) %>%
  # Add taxonomic hierarchy information
  dplyr::left_join(as.data.frame(variable_info)[, c("Genus", "Kingdom", "Phylum", "Class", "Order", "Family")], by = "Genus") %>%
  # Add fold change data
  dplyr::left_join(personalized_score_fc1, by = c("Genus" = "genus")) %>%
  # Add adjusted p-value data
  dplyr::left_join(personalized_score_fc1_p_adjust, by = c("Genus" = "genus")) %>%
  # Add second fold change data
  dplyr::left_join(personalized_score_fc2, by = c("Genus" = "genus")) %>%
  # Create body site presence indicators
  dplyr::mutate(Stool = case_when(!is.na(Stool_fc1) ~ "Stool", TRUE ~ "no")) %>%
  dplyr::mutate(Skin = case_when(!is.na(Skin_fc1) ~ "Skin", TRUE ~ "no")) %>%
  dplyr::mutate(Oral = case_when(!is.na(Oral_fc1) ~ "Oral", TRUE ~ "no")) %>%
  dplyr::mutate(Nasal = case_when(!is.na(Nasal_fc1) ~ "Nasal", TRUE ~ "no")) %>%
  # Add significance markers (* for p < 0.05)
  dplyr::mutate(Stool_fc1_star = case_when(
    !is.na(Stool_fc1_p_adjust) & Stool_fc1_p_adjust < 0.05 ~ "*",
    TRUE ~ ""
  )) %>%
  dplyr::mutate(Skin_fc1_star = case_when(
    !is.na(Skin_fc1_p_adjust) & Skin_fc1_p_adjust < 0.05 ~ "*",
    TRUE ~ ""
  )) %>%
  dplyr::mutate(Oral_fc1_star = case_when(
    !is.na(Oral_fc1_p_adjust) & Oral_fc1_p_adjust < 0.05 ~ "*",
    TRUE ~ ""
  )) %>%
  dplyr::mutate(Nasal_fc1_star = case_when(
    !is.na(Nasal_fc1_p_adjust) & Nasal_fc1_p_adjust < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Add the combined information to the tree data
raw_plot$data =
  cbind(raw_plot$data, new_info)

# ===================================================
# PHYLUM HIGHLIGHTING
# ===================================================

# Extract phylum information for highlighting different phyla on the tree
hight_data =
  raw_plot$data %>%
  # Filter for phylum-level nodes
  dplyr::filter(stringr::str_detect(label, "p__")) %>%
  dplyr::select(node, label) %>%
  # Remove the "p__" prefix for cleaner display
  dplyr::mutate(label = stringr::str_replace_all(label, "p__", "")) %>%
  dplyr::rename(id = node, type = label)

# Add phylum highlighting to the tree
# Each phylum will be highlighted with a different color
plot1 =
  raw_plot +
  ggtree::geom_hilight(
    data = hight_data,
    mapping = aes(node = id, fill = type),
    alpha = .4  # Semi-transparent highlighting
  )

# Display the tree with phylum highlighting
plot1

# ===================================================
# ADDING BODY SITE ENRICHMENT HEATMAPS
# ===================================================

# Set negative fold changes to zero for visualization purposes
# This focuses on enrichment rather than depletion
plot1$data$Stool_fc1[which(plot1$data$Stool_fc1 < 0)] = 0
plot1$data$Skin_fc1[which(plot1$data$Skin_fc1 < 0)] = 0
plot1$data$Oral_fc1[which(plot1$data$Oral_fc1 < 0)] = 0
plot1$data$Nasal_fc1[which(plot1$data$Nasal_fc1 < 0)] = 0

# Prepare data for stool enrichment heatmap
stool_fc1_info =
  plot1$data[, c("Stool_fc1", "isTip")]
rownames(stool_fc1_info) = plot1$data$label
stool_fc1_info = stool_fc1_info %>%
  dplyr::filter(isTip) %>%  # Only include leaf nodes (tips)
  dplyr::select(-isTip)

# Prepare data for skin enrichment heatmap
skin_fc1_info =
  plot1$data[, c("Skin_fc1", "isTip")]
rownames(skin_fc1_info) = plot1$data$label
skin_fc1_info = skin_fc1_info %>%
  dplyr::filter(isTip) %>%
  dplyr::select(-isTip)

# Prepare data for oral enrichment heatmap
oral_fc1_info =
  plot1$data[, c("Oral_fc1", "isTip")]
rownames(oral_fc1_info) = plot1$data$label
oral_fc1_info = oral_fc1_info %>%
  dplyr::filter(isTip) %>%
  dplyr::select(-isTip)

# Prepare data for nasal enrichment heatmap
nasal_fc1_info =
  plot1$data[, c("Nasal_fc1", "isTip")]
rownames(nasal_fc1_info) = plot1$data$label
nasal_fc1_info = nasal_fc1_info %>%
  dplyr::filter(isTip) %>%
  dplyr::select(-isTip)

# Determine the range of fold change values for consistent color scale
range(
  c(
    plot1$data$Stool_fc1,
    plot1$data$Skin_fc1,
    plot1$data$Nasal_fc1,
    plot1$data$Oral_fc1
  ),
  na.rm = TRUE
)

# Prepare for new fill scale for heatmaps
plot1 = plot1 +
  ggnewscale::new_scale_fill()

# Define color palette for different body sites
body_site_color = c(
  "Stool" = ggsci::pal_jama()(n = 7)[2],
  "Skin" = ggsci::pal_jama()(n = 7)[3],
  "Oral" = ggsci::pal_jama()(n = 7)[4],
  "Nasal" = ggsci::pal_jama()(n = 7)[5]
)

# ===================================================
# BUILDING LAYERED VISUALIZATION WITH HEATMAPS
# ===================================================

# Add stool enrichment heatmap layer
plot2 =
  gheatmap(
    p = plot1,
    data = stool_fc1_info,
    offset = -0.1,        # Position of the heatmap
    width = .08,          # Width of the heatmap
    colnames_angle = 95,  # Angle of column names
    colnames_offset_y = .25,
    colnames = FALSE,     # Don't show column names
    color = alpha(body_site_color["Stool"], 1),  # Border color for the heatmap
    legend_title = "Index1"
  ) +
  # Add color gradient for heatmap values
  scale_fill_gradientn(
    colours = c(RColorBrewer::brewer.pal(n = 11, name = "BrBG"))[-c(5:7)],
    na.value = "white",
    limits = c(0, 0.8)    # Set consistent color scale
  ) +
  # Prepare for adding significance markers
  ggnewscale::new_scale(new_aes = "fill") +
  # Add asterisks for significant genera
  geom_tippoint(
    mapping = aes(shape = Stool_fc1_star),
    x = 6.35,             # Position of the markers
    size = 2,
    show.legend = FALSE
  ) +
  scale_shape_manual(values = c("*" = "*")) +
  ggnewscale::new_scale(new_aes = "shape")

# Display the tree with stool heatmap
plot2

# Add skin enrichment heatmap layer
plot3 =
  gheatmap(
    p = plot2,
    data = skin_fc1_info,
    offset = 0.5,          # Different position for skin heatmap
    width = .08,
    colnames_angle = 95,
    colnames_offset_y = .25,
    colnames = FALSE,
    color = alpha(body_site_color["Skin"], 1),
    legend_title = "Index1"
  ) +
  scale_fill_gradientn(
    colours = c(RColorBrewer::brewer.pal(n = 11, name = "BrBG"))[-c(5:7)],
    na.value = "white",
    limits = c(0, 0.8)
  ) +
  ggnewscale::new_scale(new_aes = "fill") +
  # Add significance markers for skin
  geom_tippoint(
    mapping = aes(shape = Skin_fc1_star),
    x = 6.95,
    size = 2,
    show.legend = FALSE
  ) +
  scale_shape_manual(values = c("*" = "*")) +
  ggnewscale::new_scale(new_aes = "shape")

# Display the tree with stool and skin heatmaps
plot3

# Add oral enrichment heatmap layer
plot4 =
  gheatmap(
    p = plot3,
    data = oral_fc1_info,
    offset = 1.1,         # Different position for oral heatmap
    width = .08,
    colnames_angle = 95,
    colnames_offset_y = .25,
    colnames = FALSE,
    color = alpha(body_site_color["Oral"], 1),
    legend_title = "Index1"
  ) +
  scale_fill_gradientn(
    colours = c(RColorBrewer::brewer.pal(n = 11, name = "BrBG"))[-c(5:7)],
    na.value = "white",
    limits = c(0, 0.8)
  ) +
  ggnewscale::new_scale(new_aes = "fill") +
  # Add significance markers for oral
  geom_tippoint(
    mapping = aes(shape = Oral_fc1_star),
    x = 7.55,
    size = 2,
    show.legend = FALSE
  ) +
  scale_shape_manual(values = c("*" = "*")) +
  ggnewscale::new_scale(new_aes = "shape")

# Display the tree with stool, skin, and oral heatmaps
plot4

# Add nasal enrichment heatmap layer
plot5 =
  gheatmap(
    p = plot4,
    data = nasal_fc1_info,
    offset = 1.7,         # Different position for nasal heatmap
    width = .08,
    colnames_angle = 95,
    colnames_offset_y = .25,
    colnames = FALSE,
    color = alpha(body_site_color["Nasal"], 1),
    legend_title = "Index1"
  ) +
  scale_fill_gradientn(
    colours = c(RColorBrewer::brewer.pal(n = 11, name = "BrBG"))[-c(5:7)],
    na.value = "white",
    limits = c(0, 0.8)
  ) +
  ggnewscale::new_scale(new_aes = "fill") +
  # Add significance markers for nasal
  geom_tippoint(
    mapping = aes(shape = Nasal_fc1_star),
    x = 8.15,
    size = 2,
    show.legend = FALSE
  ) +
  scale_shape_manual(values = c("*" = "*")) +
  ggnewscale::new_scale(new_aes = "shape")

# Display the tree with all four body site heatmaps
plot5

# ===================================================
# FINISHING THE VISUALIZATION
# ===================================================

# Add genus labels to the tree tips
# Labels are colored by phylum for additional taxonomic context
plot6 =
  plot5 +
  ggnewscale::new_scale(new_aes = "color") +
  geom_tiplab(
    aes(label = label2, color = Phylum),  # Color labels by phylum
    offset = 2.5,                         # Position of the labels
    size = 2,
    show.legend = FALSE
  )

# Display the final tree visualization
plot6

# Save the final visualization as a PDF file
ggsave(plot6,
       filename = "microbiome_tree.pdf",
       width = 14,
       height = 14)
