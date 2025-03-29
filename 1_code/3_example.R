library(tidyverse)
library(phyloseq)

setwd(r4projects::get_project_wd())
rm(list = ls())

personalized_score = readr::read_csv("2_demo_data/personalized_score.csv")

head(personalized_score)

load("2_demo_data/microbiome_data.RData")

microbiome_data

####HMP
if (!requireNamespace("microbiomeViz", quietly = TRUE)) {
  remotes::install_github("lch14forever/microbiomeViz")
}
library(ggtreeExtra)
library(ggtree)
if (!requireNamespace("phyloseq", quietly = TRUE)) {
  BiocManager::install("phyloseq")
}
library(phyloseq)
library(tidyverse)

###create tree
tree =
  microbiomeViz::parsePhyloseq(
    physeq = microbiome_data,
    use_abundance = FALSE,
    node.size.scale = 0,
    node.size.offset = 0
  )

raw_plot <-
  tree.backbone(
    tree = tree,
    size = 0.3,
    shape = 16,
    layout = "circular",
    fill = "black",
    color = "black"
  )

raw_plot

raw_plot$data$nodeClass2 =
  as.character(raw_plot$data$nodeClass)

raw_plot$data$nodeClass2[is.na(raw_plot$data$nodeClass2)] = "Root"
raw_plot$data$nodeClass2[raw_plot$data$nodeClass2 == "f"] = "Family"
raw_plot$data$nodeClass2[raw_plot$data$nodeClass2 == "c"] = "Class"
raw_plot$data$nodeClass2[raw_plot$data$nodeClass2 == "o"] = "Order"
raw_plot$data$nodeClass2[raw_plot$data$nodeClass2 == "p"] = "Phylum"
raw_plot$data$nodeClass2[raw_plot$data$nodeClass2 == "g"] = "Genus"
raw_plot$data$nodeClass2[raw_plot$data$nodeClass2 == "k"] = "Kingdom"

raw_plot$data$nodeSize2 = raw_plot$data$nodeSize
raw_plot$data$nodeSize2[raw_plot$data$nodeClass2 == "Root"] = 3.5
raw_plot$data$nodeSize2[raw_plot$data$nodeClass2 == "Kingdom"] = 3
raw_plot$data$nodeSize2[raw_plot$data$nodeClass2 == "Phylum"] = 2.5
raw_plot$data$nodeSize2[raw_plot$data$nodeClass2 == "Class"] = 2
raw_plot$data$nodeSize2[raw_plot$data$nodeClass2 == "Order"] = 1.5
raw_plot$data$nodeSize2[raw_plot$data$nodeClass2 == "Family"] = 1
raw_plot$data$nodeSize2[raw_plot$data$nodeClass2 == "Genus"] = 0.5

# #####add node point
raw_plot =
  raw_plot +
  ggtree::geom_point2(aes(color = nodeClass2, size = nodeSize2), show.legend = FALSE) +
  ggnewscale::new_scale(new_aes = "fill")

raw_plot

###label 2 is the name of taxa
raw_plot$data$label2 =
  raw_plot$data$label %>%
  stringr::str_split("__") %>%
  purrr::map(function(x) {
    x[2]
  }) %>%
  unlist()

raw_plot$data$label2[as.character(raw_plot$data$nodeClass) != "g"] = NA

#####add new information
personalized_score_fc1 =
  personalized_score %>%
  dplyr::select(genus, class, fc1) %>%
  tidyr::pivot_wider(names_from = class, values_from = "fc1")

colnames(personalized_score_fc1)[-1] = paste(colnames(personalized_score_fc1)[-1], "fc1", sep = "_")

personalized_score_fc2 =
  personalized_score %>%
  dplyr::select(genus, class, fc2) %>%
  tidyr::pivot_wider(names_from = class, values_from = "fc2")

colnames(personalized_score_fc2)[-1] = paste(colnames(personalized_score_fc2)[-1], "fc2", sep = "_")

personalized_score_fc1_p_adjust =
  personalized_score %>%
  dplyr::select(genus, class, fc1_p_adjust) %>%
  tidyr::pivot_wider(names_from = class, values_from = "fc1_p_adjust")

colnames(personalized_score_fc1_p_adjust)[-1] =
  paste(colnames(personalized_score_fc1_p_adjust)[-1],
        "fc1_p_adjust",
        sep = "_")

variable_info <-
  tax_table(microbiome_data)

new_info =
  data.frame(Genus = raw_plot$data$label2) %>%
  dplyr::left_join(as.data.frame(variable_info)[, c("Genus", "Kingdom", "Phylum", "Class", "Order", "Family")], by = "Genus") %>%
  dplyr::left_join(personalized_score_fc1, by = c("Genus" = "genus")) %>%
  dplyr::left_join(personalized_score_fc1_p_adjust, by = c("Genus" = "genus")) %>%
  dplyr::left_join(personalized_score_fc2, by = c("Genus" = "genus")) %>%
  dplyr::mutate(Stool = case_when(!is.na(Stool_fc1) ~ "Stool", TRUE ~ "no")) %>%
  dplyr::mutate(Skin = case_when(!is.na(Skin_fc1) ~ "Skin", TRUE ~ "no")) %>%
  dplyr::mutate(Oral = case_when(!is.na(Oral_fc1) ~ "Oral", TRUE ~ "no")) %>%
  dplyr::mutate(Nasal = case_when(!is.na(Nasal_fc1) ~ "Nasal", TRUE ~ "no")) %>%
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

raw_plot$data =
  cbind(raw_plot$data, new_info)

######add highlight
hight_data =
  raw_plot$data %>%
  dplyr::filter(stringr::str_detect(label, "p__")) %>%
  dplyr::select(node, label) %>%
  dplyr::mutate(label = stringr::str_replace_all(label, "p__", "")) %>%
  dplyr::rename(id = node, type = label)

plot1 =
  raw_plot +
  ggtree::geom_hilight(
    data = hight_data,
    mapping = aes(node = id, fill = type),
    alpha = .4
  )

plot1

######add heatmap
##add fc1 information
plot1$data$Stool_fc1[which(plot1$data$Stool_fc1 < 0)] = 0
plot1$data$Skin_fc1[which(plot1$data$Skin_fc1 < 0)] = 0
plot1$data$Oral_fc1[which(plot1$data$Oral_fc1 < 0)] = 0
plot1$data$Nasal_fc1[which(plot1$data$Nasal_fc1 < 0)] = 0

##add multiple tip information
stool_fc1_info =
  plot1$data[, c("Stool_fc1", "isTip")]
rownames(stool_fc1_info) = plot1$data$label
stool_fc1_info = stool_fc1_info %>%
  dplyr::filter(isTip) %>%
  dplyr::select(-isTip)

skin_fc1_info =
  plot1$data[, c("Skin_fc1", "isTip")]
rownames(skin_fc1_info) = plot1$data$label
skin_fc1_info = skin_fc1_info %>%
  dplyr::filter(isTip) %>%
  dplyr::select(-isTip)

oral_fc1_info =
  plot1$data[, c("Oral_fc1", "isTip")]
rownames(oral_fc1_info) = plot1$data$label
oral_fc1_info = oral_fc1_info %>%
  dplyr::filter(isTip) %>%
  dplyr::select(-isTip)

nasal_fc1_info =
  plot1$data[, c("Nasal_fc1", "isTip")]
rownames(nasal_fc1_info) = plot1$data$label
nasal_fc1_info = nasal_fc1_info %>%
  dplyr::filter(isTip) %>%
  dplyr::select(-isTip)

range(
  c(
    plot1$data$Stool_fc1,
    plot1$data$Skin_fc1,
    plot1$data$Nasal_fc1,
    plot1$data$Oral_fc1
  ),
  na.rm = TRUE
)

plot1 = plot1 +
  ggnewscale::new_scale_fill()

body_site_color = c(
  "Stool" = ggsci::pal_jama()(n = 7)[2],
  "Skin" = ggsci::pal_jama()(n = 7)[3],
  "Oral" = ggsci::pal_jama()(n = 7)[4],
  "Nasal" = ggsci::pal_jama()(n = 7)[5]
)

plot2 =
  gheatmap(
    p = plot1,
    data = stool_fc1_info,
    offset = -0.1,
    width = .08,
    colnames_angle = 95,
    colnames_offset_y = .25,
    colnames = FALSE,
    color = alpha(body_site_color["Stool"], 1),
    legend_title = "Index1"
  ) +
  scale_fill_gradientn(
    colours = c(RColorBrewer::brewer.pal(n = 11, name = "BrBG"))[-c(5:7)],
    na.value = "white",
    limits = c(0, 0.8)
  ) +
  ggnewscale::new_scale(new_aes = "fill") +
  geom_tippoint(
    mapping = aes(shape = Stool_fc1_star),
    x = 6.35,
    size = 2,
    show.legend = FALSE
  ) +
  scale_shape_manual(values = c("*" = "*")) +
  ggnewscale::new_scale(new_aes = "shape")

plot2

plot3 =
  gheatmap(
    p = plot2,
    data = skin_fc1_info,
    offset = 0.5,
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
  geom_tippoint(
    mapping = aes(shape = Skin_fc1_star),
    x = 6.95,
    size = 2,
    show.legend = FALSE
  ) +
  scale_shape_manual(values = c("*" = "*")) +
  ggnewscale::new_scale(new_aes = "shape")

plot3

plot4 =
  gheatmap(
    p = plot3,
    data = oral_fc1_info,
    offset = 1.1,
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
  geom_tippoint(
    mapping = aes(shape = Oral_fc1_star),
    x = 7.55,
    size = 2,
    show.legend = FALSE
  ) +
  scale_shape_manual(values = c("*" = "*")) +
  ggnewscale::new_scale(new_aes = "shape")

plot4

plot5 =
  gheatmap(
    p = plot4,
    data = nasal_fc1_info,
    offset = 1.7,
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
  geom_tippoint(
    mapping = aes(shape = Nasal_fc1_star),
    x = 8.15,
    size = 2,
    show.legend = FALSE
  ) +
  scale_shape_manual(values = c("*" = "*")) +
  ggnewscale::new_scale(new_aes = "shape")

plot5

##add tip lab
##only add some tip points
# idx1 = which(!is.na(raw_plot$data$Stool_fc1_p_adjust) & raw_plot$data$Stool_fc1_p_adjust < 0.001)
# idx1 = sample(idx1, 10)
# raw_plot$data$label2[-idx1] = NA
plot6 =
  plot5 +
  ggnewscale::new_scale(new_aes = "color") +
  geom_tiplab(
    aes(label = label2, color = Phylum),
    offset = 2.5,
    size = 2,
    show.legend = FALSE
  )

plot6

ggsave(plot6,
       filename = "microbiome_tree.pdf",
       width = 14,
       height = 14)
