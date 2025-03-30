# BiocManager is a package manager for Bioconductor
if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

####load required packages
###remotes
if (!require(remotes)) {
  install.packages("remotes")
}

if (!require("tidyverse", quietly = TRUE)){
  install.packages("tidyverse")
}

##readr to read data (csv, tsv, etc.)
if (!require("readr", quietly = TRUE)){
  install.packages("readr")
}
##tibble to create data.frame
if (!require("tibble", quietly = TRUE)){
  install.packages("tibble")
}
# shatdowtext is a package for adding shadow to text
if (!require("shadowtext", quietly = TRUE)){
  install.packages("shadowtext")
}
# ggsci is a package for scientific journal color palettes
if (!require("ggsci", quietly = TRUE)){
  install.packages("ggsci")
}

# igraph is a package for creating, manipulating, and analyzing network data
if (!require("igraph", quietly = TRUE)){
  install.packages("igraph")
}
# ggraph is a package for creating network visualizations using ggplot2
if(!require("ggraph", quietly = TRUE)){
  install.packages("ggraph")
}
#tidygraph is a package for tidy manipulation of graph data
if(!require("tidygraph", quietly = TRUE)){
  install.packages("tidygraph")
}
# ggnetword is another package for creating network visualizations using ggplot2
if(!require("ggnetwork", quietly = TRUE)){
  install.packages("ggnetwork")
}
# extrafont is a package for using extra fonts in ggplot2
if (!require(extrafont)) {
  install.packages("extrafont")
}

# microbiomeViz is a package for visualizing microbiome data
if (!require("microbiomeViz", quietly = TRUE)) {
  remotes::install_github("lch14forever/microbiomeViz")
}
# phyloseq is a package for microbiome data analysis
if (!require("phyloseq", quietly = TRUE)) {
  BiocManager::install("phyloseq")
}
# ggtree is a package for visualizing phylogenetic trees
if(!require("ggtree", quietly = TRUE)){
  BiocManager::install("ggtree")  
}
# treeio is a package for reading and writing phylogenetic tree data
if(!require("treeio", quietly = TRUE)){
  BiocManager::install("treeio")  
}
# ggtreeExtra is a package for adding extra features to ggtree
if(!require("ggtreeExtra", quietly = TRUE)){
  BiocManager::install("ggtreeExtra")  
}

####r4projects
if (!require(r4projects)) {
  remotes::install_github("jaspershen/r4projects")
}
