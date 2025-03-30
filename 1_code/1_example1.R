####load required packages
###remotes
if (!require(remotes)) {
  install.packages("remotes")
}

####r4projects
if (!require(r4projects)) {
  remotes::install_github("jaspershen/r4projects")
}

#######set to root directory
setwd(r4projects::get_project_wd())

###load data
if (!require(readr)) {
  install.packages("readr")
}

node_data <-
  readr::read_csv(
    "https://raw.githubusercontent.com/jaspershen-lab/MD7115-NTU/main/2_demo_data/example_node_data.csv"
  )

edge_data <-
  readr::read_csv(
    "https://raw.githubusercontent.com/jaspershen-lab/MD7115-NTU/main/2_demo_data/example_edge_data.csv"
  )

head(node_data)

head(edge_data)

total_graph <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = FALSE) %>%
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))


#####up-down
g <- total_graph
library(igraph)
library(ggraph)

V(g)$type <- bipartite_mapping(g)$type

coords <-
  create_layout(g, layout = "bipartite") %>%
  dplyr::select(node, class, mol_name, class2, x, y)

coords$y[coords$class2 == "wearable"] = 0.3
coords$y[coords$class2 == "internal-omics"] = 1

table(coords$class)

coords$y[coords$class2 == "wearable"] = -0.5
coords$y[coords$class == "metabolomics"] = 1
coords$y[coords$class == "lipidomics"] = 0.5
coords$y[coords$class == "proteomics"] = 0
coords$y[coords$class == "cytokine"] = 0
coords$y[coords$class == "metabolic_panel"] = 0
coords$y[coords$class == "total_protein"] = 0

coords <-
  coords %>%
  dplyr::select(x, y) %>%
  dplyr::mutate(
    theta = x / (max(x) + 1) * 2 * pi,
    r = y + 1,
    x = r * cos(theta),
    y = r * sin(theta)
  )

my_graph <-
  create_layout(
    graph = g,
    layout = "manual",
    x = coords$x,
    y = coords$y
    # node.position = coords
  )

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

wearable_color =
  c(
    "sleep" = ggsci::pal_d3()(n = 10)[2],
    "cgm"  = ggsci::pal_d3()(n = 10)[6],
    "hr" = ggsci::pal_d3()(n = 10)[7],
    "step"  = ggsci::pal_d3()(n = 10)[9],
    "food" = ggsci::pal_d3()(n = 10)[10]
  )

plot <-
  ggraph(my_graph, layout = 'bipartite') +
  geom_edge_diagonal(aes(
    color = lagged_cor,
    width = -log(lagged_cor_p_adjust, 10)
  ),
  alpha = 0.5,
  show.legend = TRUE) +
  geom_node_point(
    aes(fill = class, size = Degree),
    shape = 21,
    alpha = 0.5,
    show.legend = TRUE
  ) +
  shadowtext::geom_shadowtext(
    aes(
      x = x,
      y = y,
      label = ifelse(class2 == "wearable", mol_name, NA),
      color = class
    ),
    bg.color = "white",
    size = 5,
    show.legend = FALSE
  ) +
  shadowtext::geom_shadowtext(
    aes(
      x = x,
      y = y,
      label = ifelse(class2 == "wearable", NA, mol_name),
      color = class,
      nudge_y = ifelse(class2 == "wearable", "inward", 'outward'),
      angle = -((-node_angle(x, y) + 90) %% 180) + 90
    ),
    bg.color = "white",
    size = 2,
    show.legend = FALSE,
    check_overlap = TRUE
  ) +
  scale_size_continuous(range = c(3, 10)) +
  scale_fill_manual(values = c(class_color, wearable_color)) +
  scale_color_manual(values = c(class_color, wearable_color)) +
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
  scale_edge_width_continuous(range = c(0.3, 2)) +
  ggraph::scale_edge_color_gradientn(colours = c(alpha("#3B4992FF", 0.7), "white", alpha("#EE0000FF", 0.7))) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  )

plot

extrafont::loadfonts()

ggsave(plot,
       filename = "example_network.pdf",
       width = 9.1,
       height = 7)
