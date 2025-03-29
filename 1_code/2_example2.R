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

###read data
pathway_result <-
  readr::read_csv(
    "https://github.com/jaspershen-lab/MD7115-NTU/blob/main/2_demo_data/pathway_result.csv"
  )

##get the similartiy matrix
######get the similarity between GO terms

message("Calculating similartiy matrix...")

if (any(dir(file.path(path, "intermediate_data")) == "sim_matrix")) {
  load(file.path(path, "intermediate_data/sim_matrix"))
} else{
  if (database == "go") {
    sim_matrix <-
      get_go_result_sim(
        result = result,
        sim.cutoff = sim.cutoff,
        measure_method = measure_method
      )
  }
  
  if (database == "kegg") {
    sim_matrix <-
      tryCatch(
        sim_matrix <-
          simplifyEnrichment::term_similarity_from_KEGG(term_id = c(result$ID), method = "jaccard") %>%
          as.data.frame() %>%
          tibble::rownames_to_column(var = "name1") %>%
          tidyr::pivot_longer(
            cols = -name1,
            names_to = "name2",
            values_to = "sim"
          ) %>%
          dplyr::filter(name1 != name2),
        error = function(x) {
          data.frame(name1 = character(),
                     name2 = character(),
                     sim = numeric())
        }
      )
  }
  
  
  if (database == "reactome") {
    sim_matrix <-
      tryCatch(
        sim_matrix <-
          simplifyEnrichment::term_similarity_from_Reactome(term_id = c(result$ID), method = "jaccard") %>%
          as.data.frame() %>%
          tibble::rownames_to_column(var = "name1") %>%
          tidyr::pivot_longer(
            cols = -name1,
            names_to = "name2",
            values_to = "sim"
          ) %>%
          dplyr::filter(name1 != name2),
        error = function(x) {
          data.frame(name1 = character(),
                     name2 = character(),
                     sim = numeric())
        }
      )
  }
  
  save(sim_matrix, file = file.path(path, "intermediate_data/sim_matrix"))
}

####module detection
message("Identifying modules...")

if (any(dir(file.path(path, "intermediate_data")) == "graph_data")) {
  load(file.path(path, "intermediate_data/graph_data"))
} else{
  edge_data <-
    rbind(sim_matrix) %>%
    dplyr::rename(from = name1, to = name2) %>%
    dplyr::filter(sim > sim.cutoff)
  
  node_data <-
    rbind(result) %>%
    as.data.frame() %>%
    dplyr::select(ID, everything()) %>%
    dplyr::rename(node = ID)
  
  graph_data <-
    tidygraph::tbl_graph(nodes = node_data,
                         edges = edge_data,
                         directed = FALSE) %>%
    dplyr::mutate(degree = tidygraph::centrality_degree())
  
  subnetwork <-
    igraph::cluster_edge_betweenness(graph = graph_data, weights = abs(edge_attr(graph_data, "sim")))
  
  # save(subnetwork, file = file.path(path, "subnetwork"))
  cluster <-
    paste(database, "Module", as.character(igraph::membership(subnetwork)), sep = "_")
  
  graph_data <-
    graph_data %>%
    tidygraph::mutate(module = cluster)
  
  ###clustered different GO terms
  result_with_module <-
    igraph::vertex_attr(graph_data) %>%
    do.call(cbind, .) %>%
    as.data.frame() %>%
    dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
    dplyr::arrange(module, p.adjust)
  
  if (database == "go") {
    result_with_module <-
      result_with_module %>%
      dplyr::arrange(ONTOLOGY, module, p.adjust)
  }
  
  ###add module content number
  module_content_number <-
    result_with_module %>%
    dplyr::count(module) %>%
    dplyr::rename(module_content_number = n)
  
  result_with_module <-
    result_with_module %>%
    dplyr::left_join(module_content_number, by = "module")
  
  save(result_with_module,
       file = file.path(path, "intermediate_data/result_with_module"))
  graph_data <-
    graph_data %>%
    activate(what = "nodes") %>%
    dplyr::left_join(module_content_number, by = "module")
  
  save(graph_data, file = file.path(path, "intermediate_data/graph_data"))
}

###plot to show the clusters of GO terms
cluster_label_module <-
  igraph::as_data_frame(graph_data, what = "vertices") %>%
  dplyr::group_by(module) %>%
  dplyr::filter(p.adjust == min(p.adjust) &
                  Count == max(Count)) %>%
  dplyr::slice_head(n = 1) %>%
  pull(Description)

cluster_label_all <-
  igraph::as_data_frame(graph_data, what = "vertices")$Description

plot <-
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
  ) +
  geom_node_text(aes(
    x = x,
    y = y,
    label = ifelse(Description %in% cluster_label_module, Description, NA)
  ),
  size = 3,
  repel = TRUE) +
  guides(fill = guide_legend(ncol = 1)) +
  scale_edge_width_continuous(range = c(0.1, 2)) +
  scale_size_continuous(range = c(1, 7)) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  )

plot

extrafont::loadfonts()

ggplot2::ggsave(
  plot,
  filename =
    file.path(path, "similarity_network_plot.pdf"),
  width = 9,
  height = 7
)

###output some files
result_with_module <-
  igraph::vertex_attr(graph_data) %>%
  do.call(cbind, .) %>%
  as.data.frame() %>%
  dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
  dplyr::arrange(module, p.adjust)

if (database == "go") {
  result_with_module <-
    result_with_module %>%
    dplyr::arrange(ONTOLOGY, module, p.adjust)
}

module_result <-
  result_with_module %>%
  plyr::dlply(.variables = .(module)) %>%
  purrr::map(function(x) {
    # cat(unique(x$module), " ")
    if (nrow(x) == 1) {
      return(x)
    }
    
    x =
      x %>%
      dplyr::arrange(p.adjust)
    
    x$node <-
      paste(x$node, collapse = ";")
    
    x$Description <-
      paste(x$Description, collapse = ";")
    
    x$BgRatio <-
      paste(x$BgRatio, collapse = ";")
    
    x$pvalue <- min(as.numeric(x$pvalue))
    x$p.adjust <- min(as.numeric(x$p.adjust))
    x$qvalue <- min(as.numeric(x$qvalue))
    x$geneID =
      x$geneID %>%
      stringr::str_split(pattern = "/") %>%
      unlist() %>%
      unique() %>%
      paste(collapse = '/')
    
    x$Count <-
      length(stringr::str_split(x$geneID[1], pattern = "/")[[1]])
    
    x =
      x %>%
      dplyr::select(module, everything()) %>%
      dplyr::distinct(module, .keep_all = TRUE)
    
    x
    
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  dplyr::mutate(
    module_annotation = case_when(
      module == "Other" ~ Description,
      module != "Other" ~ stringr::str_split(Description, ";")[[1]][1]
    )
  ) %>%
  dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::select(module_annotation, everything())

module_result$module_annotation <-
  stringr::str_split(module_result$Description, ";") %>%
  purrr::map(function(x) {
    x[1]
  }) %>%
  unlist()

save(module_result,
     file = file.path(path, "intermediate_data/module_result"))

wb = openxlsx::createWorkbook()
openxlsx::modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roma")
addWorksheet(wb, sheetName = "enriched_pathway_result", gridLines = TRUE)
addWorksheet(wb, sheetName = "enriched_module_result", gridLines = TRUE)
freezePane(wb,
           sheet = 1,
           firstRow = TRUE,
           firstCol = TRUE)
freezePane(wb,
           sheet = 2,
           firstRow = TRUE,
           firstCol = TRUE)
writeDataTable(
  wb,
  sheet = 1,
  x = result_with_module,
  colNames = TRUE,
  rowNames = FALSE
)

writeDataTable(
  wb,
  sheet = 2,
  x = module_result,
  colNames = TRUE,
  rowNames = FALSE
)

saveWorkbook(wb,
             file = file.path(path, "enriched_result.xlsx"),
             overwrite = TRUE)


####output some results
dir.create(file.path(path, "Similarity_plot"),
           recursive = TRUE,
           showWarnings = FALSE)

message("Output module plot...")

for (temp_cluster in unique(result_with_module$module)) {
  cat(temp_cluster, " ")
  
  if (sum(result_with_module$module == temp_cluster) == 1) {
    next()
  }
  
  plot1 <-
    graph_data %>%
    tidygraph::filter(module == temp_cluster) %>%
    ggraph(layout = 'fr', circular = FALSE) +
    geom_edge_link(
      aes(width = sim),
      strength = 1,
      color = "black",
      alpha = 1,
      show.legend = TRUE
    ) +
    geom_node_point(
      aes(fill = -log(p.adjust, 10), size = Count),
      shape = 21,
      alpha = 1,
      show.legend = TRUE
    ) +
    geom_node_text(aes(x = x, y = y, label = Description), check_overlap = TRUE) +
    guides(fill = guide_legend(ncol = 1)) +
    scale_edge_width_continuous(range = c(0.1, 2)) +
    scale_size_continuous(range = c(3, 10)) +
    ggraph::theme_graph() +
    theme(
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      legend.position = "left",
      legend.background = element_rect(fill = "transparent", color = NA)
    )
  
  plot1
  
  plot2 <-
    result_with_module %>%
    dplyr::filter(module == temp_cluster) %>%
    dplyr::mutate(p.adjust = -log(as.numeric(p.adjust, 10))) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::mutate(Description = factor(Description, levels = Description)) %>%
    ggplot(aes(p.adjust, Description)) +
    geom_bar(stat = "identity", fill = "black") +
    geom_text(
      aes(x = 0, Description, label = Description),
      hjust = 0,
      size = 5,
      color = "red"
    ) +
    theme_bw() +
    labs(y = "", x = "-log10(FDR adjusted P value)") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  temp_data =
    result_with_module %>%
    dplyr::filter(module == temp_cluster) %>%
    dplyr::mutate(p.adjust = -log(as.numeric(p.adjust, 10))) %>%
    dplyr::select(Description, p.adjust) %>%
    dplyr::mutate(Description = stringr::str_replace_all(Description, ",", "")) %>%
    plyr::dlply(.variables = .(Description)) %>%
    purrr::map(function(x) {
      data.frame(word = stringr::str_split(x$Description, " ")[[1]],
                 p.adjust = x$p.adjust)
    }) %>%
    do.call(rbind, .) %>%
    as.data.frame() %>%
    plyr::dlply(.variables = .(word)) %>%
    purrr::map(function(x) {
      x$p.adjust <- sum(x$p.adjust)
      x %>%
        dplyr::distinct(word, .keep_all = TRUE)
    }) %>%
    do.call(rbind, .) %>%
    as.data.frame() %>%
    dplyr::filter(!word %in% remove_words)
  
  plot3 =
    temp_data %>%
    ggplot(aes(label = word, size = p.adjust)) +
    geom_text_wordcloud() +
    scale_radius(range = c(5, 15), limits = c(0, NA)) +
    theme_minimal()
  
  plot =
    plot1 + plot2 + plot3 + patchwork::plot_layout(nrow = 1)
  
  ggsave(
    plot,
    filename = file.path(
      path,
      "Similarity_plot",
      paste(temp_cluster, "sim_plot.pdf", sep = "_")
    ),
    width = 21,
    height = 7
  )
}


# ##matrix tow show the cluster GO terms
# result_with_module %>%
#   dplyr::group_by(ONTOLOGY) %>%
#   dplyr::summarise(n = n()) %>%
#   dplyr::mutate(n = n * 10 / max(n) + 2)
#output the correlation matrix

# if(database == "go"){
#   message("Output correlation matrix plot...")
#   for(ont in c('MF', "BP", "CC")) {
#     cat(ont, " ")
#     show_matrix_cluster(
#       result = result_with_module %>% dplyr::mutate(Direction = "UP"),
#       ont = ont,
#       measure = "Wang",
#       remove_words = remove_words,
#       margin = 15,
#       width = 14,
#       height = 8,
#       path = path,
#       top = 15
#     )
#   }
# }

###output the cluster annotation for each cluster
if (database == "go") {
  dir.create(file.path(path, "GO_module_graph"), showWarnings = FALSE)
  
  unique(module_result$module) %>%
    purrr::map(
      .f = function(x) {
        cat(x, " ")
        number <- module_result %>%
          dplyr::filter(module == x) %>%
          pull(module_content_number) %>%
          as.numeric()
        if (number == 1) {
          return(NULL)
        }
        
        temp_id <-
          module_result %>%
          dplyr::filter(module == x) %>%
          dplyr::pull(node) %>%
          stringr::str_split(";") %>%
          `[[`(1) %>%
          pRoloc::goIdToTerm(keepNA = FALSE) %>%
          data.frame(id = ., class = "YES") %>%
          tibble::rownames_to_column(var = "name")
        
        temp_plot =
          GOSim::getGOGraph(term = temp_id$name, prune = Inf) %>%
          igraph::igraph.from.graphNEL() %>%
          tidygraph::as_tbl_graph() %>%
          tidygraph::left_join(temp_id, by = "name") %>%
          dplyr::mutate(class = case_when(is.na(class) ~ "NO", TRUE ~ class))
        
        plot =
          temp_plot %>%
          ggraph(layout = 'kk', circular = FALSE) +
          geom_edge_link(
            color = ggsci::pal_aaas()(n = 10)[1],
            alpha = 1,
            arrow = grid::arrow(
              angle = 10,
              length = unit(0.2, "inches"),
              type = "closed"
            ),
            show.legend = FALSE
          ) +
          geom_node_point(
            aes(fill = class),
            shape = 21,
            alpha = 1,
            size = 6,
            show.legend = FALSE
          ) +
          geom_node_text(aes(
            x = x,
            y = y,
            label = ifelse(class == "YES", id, NA)
          ),
          size = 3,
          repel = TRUE) +
          scale_fill_manual(values = c('YES' = "red", 'NO' = "white")) +
          ggraph::theme_graph() +
          theme(
            plot.background = element_rect(fill = "transparent", color = NA),
            panel.background = element_rect(fill = "transparent", color = NA),
            legend.position = "left",
            legend.background = element_rect(fill = "transparent", color = NA)
          )
        # plot
        
        ggsave(
          plot,
          filename = file.path(path, "GO_module_graph", paste(x, "_GO graph.pdf", sep = "")),
          width = 7,
          height = 7
        )
      }
    )
  
}
}

combine_database <-
  function(module_result_go,
           module_result_kegg,
           module_result_reactome,
           variable_info,
           sim.cutoff = 0.5,
           measure_method = c("jaccard"),
           path = ".") {
    ######calculate the similarity (jaccard index) between all the pathways
    
    if (!is.null(module_result_go)) {
      module_result_go <-
        module_result_go %>%
        dplyr::filter(ONTOLOGY != "CC") %>%
        dplyr::arrange(p.adjust) %>%
        dplyr::mutate(database = "GO") %>%
        dplyr::select(
          module_annotation,
          module,
          Description,
          BgRatio,
          pvalue,
          qvalue,
          p.adjust,
          Count,
          database,
          geneID,
          pathway_id = node
        ) %>%
        dplyr::mutate(Count = as.numeric(Count)) %>%
        dplyr::filter(!is.na(module_annotation))
    }
    
    if (!is.null(module_result_kegg)) {
      module_result_kegg <-
        module_result_kegg %>%
        dplyr::arrange(p.adjust) %>%
        dplyr::mutate(database = "KEGG") %>%
        dplyr::select(
          module_annotation,
          module,
          Description,
          BgRatio,
          pvalue,
          qvalue,
          p.adjust,
          Count,
          database,
          geneID,
          pathway_id = node
        ) %>%
        dplyr::mutate(Count = as.numeric(Count)) %>%
        dplyr::filter(!is.na(module_annotation))
    }
    
    if (!is.null(module_result_reactome)) {
      module_result_reactome <-
        module_result_reactome %>%
        dplyr::arrange(p.adjust) %>%
        dplyr::mutate(database = "Reactome") %>%
        dplyr::select(
          module_annotation,
          module,
          Description,
          BgRatio,
          pvalue,
          qvalue,
          p.adjust,
          Count,
          database,
          geneID,
          pathway_id = node
        ) %>%
        dplyr::mutate(Count = as.numeric(Count)) %>%
        dplyr::filter(!is.na(module_annotation))
    }
    
    jaccard_index <-
      get_jaccard_index_for_three_databases(
        module_result_go = module_result_go,
        module_result_kegg = module_result_kegg,
        module_result_reactome = module_result_reactome,
        variable_info = variable_info
      )
    
    edge_data =
      jaccard_index %>%
      dplyr::filter(value > 0.5) %>%
      dplyr::rename(from = name1, to = name2, sim = value)
    
    node_data <-
      rbind(module_result_go,
            module_result_kegg,
            module_result_reactome) %>%
      dplyr::select(module, everything()) %>%
      dplyr::rename(node = module)
    
    graph_data <-
      tidygraph::tbl_graph(nodes = node_data,
                           edges = edge_data,
                           directed = FALSE) %>%
      dplyr::mutate(degree = tidygraph::centrality_degree())
    
    subnetwork <-
      igraph::cluster_edge_betweenness(graph = graph_data, weights = abs(edge_attr(graph_data, "sim")))
    cluster <-
      paste("Function_module", as.character(igraph::membership(subnetwork)), sep = "_")
    
    graph_data <-
      graph_data %>%
      tidygraph::mutate(module = cluster)
    
    ###clustered different GO terms
    result_with_module <-
      igraph::vertex_attr(graph_data) %>%
      do.call(cbind, .) %>%
      as.data.frame() %>%
      dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
      dplyr::arrange(module, p.adjust)
    
    ##add module content number
    module_content_number <-
      result_with_module %>%
      dplyr::count(module) %>%
      dplyr::rename(module_content_number = n)
    
    result_with_module <-
      result_with_module %>%
      dplyr::left_join(module_content_number, by = "module")
    
    dir.create(
      file.path(path, "intermediate_data"),
      showWarnings = FALSE,
      recursive = TRUE
    )
    save(result_with_module,
         file = file.path(path, "intermediate_data/result_with_module"))
    
    graph_data <-
      graph_data %>%
      activate(what = "nodes") %>%
      dplyr::left_join(module_content_number, by = "module")
    
    save(graph_data, file = file.path(path, "intermediate_data/graph_data"))
    
    module_result <-
      result_with_module %>%
      plyr::dlply(.variables = .(module)) %>%
      purrr::map(function(x) {
        # cat(unique(x$module), " ")
        if (nrow(x) == 1) {
          x$module_content <-
            paste(x$node, collapse = ";")
          x <-
            x %>%
            dplyr::select(module, everything()) %>%
            dplyr::distinct(module, .keep_all = TRUE) %>%
            dplyr::select(-node)
          return(x)
        }
        
        x =
          x %>%
          dplyr::arrange(p.adjust)
        
        x$module_content <-
          paste(x$node, collapse = ";")
        
        x$Description <-
          paste(x$Description, collapse = ";")
        
        x$BgRatio <-
          paste(x$BgRatio, collapse = ";")
        
        x$pvalue <- min(as.numeric(x$pvalue))
        x$p.adjust <- min(as.numeric(x$p.adjust))
        x$qvalue <- min(as.numeric(x$qvalue))
        x$geneID =
          x$geneID %>%
          stringr::str_split(pattern = "/") %>%
          unlist() %>%
          unique() %>%
          paste(collapse = '/')
        
        x$Count <-
          length(stringr::str_split(x$geneID[1], pattern = "/")[[1]])
        
        x <-
          x %>%
          dplyr::select(module, everything()) %>%
          dplyr::distinct(module, .keep_all = TRUE) %>%
          dplyr::select(-node)
        x
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame() %>%
      dplyr::mutate(
        module_annotation = case_when(
          module == "Other" ~ Description,
          module != "Other" ~ stringr::str_split(Description, ";")[[1]][1]
        )
      ) %>%
      dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::select(module_annotation, everything())
    
    module_result$module_annotation <-
      stringr::str_split(module_result$Description, ";") %>%
      purrr::map(function(x) {
        x[1]
      }) %>%
      unlist()
    
    save(module_result,
         file = file.path(path, "intermediate_data/module_result"))
    
    wb = openxlsx::createWorkbook()
    openxlsx::modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roma")
    addWorksheet(wb, sheetName = "enriched_pathway_result", gridLines = TRUE)
    addWorksheet(wb, sheetName = "enriched_module_result", gridLines = TRUE)
    freezePane(wb,
               sheet = 1,
               firstRow = TRUE,
               firstCol = TRUE)
    freezePane(wb,
               sheet = 2,
               firstRow = TRUE,
               firstCol = TRUE)
    writeDataTable(
      wb,
      sheet = 1,
      x = result_with_module,
      colNames = TRUE,
      rowNames = FALSE
    )
    
    writeDataTable(
      wb,
      sheet = 2,
      x = module_result,
      colNames = TRUE,
      rowNames = FALSE
    )
    
    saveWorkbook(wb,
                 file = file.path(path, "enriched_result.xlsx"),
                 overwrite = TRUE)
    
    ###plot to show the clusters of GO terms
    cluster_label_module <-
      igraph::as_data_frame(graph_data, what = "vertices") %>%
      dplyr::group_by(module) %>%
      dplyr::filter(p.adjust == min(p.adjust) &
                      Count == max(Count)) %>%
      dplyr::slice_head(n = 1) %>%
      pull(Description)
    
    cluster_label_all <-
      node_data$Description
    
    plot <-
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
        aes(fill = database, size = -log(p.adjust, 10)),
        shape = 21,
        alpha = 1,
        show.legend = TRUE
      ) +
      geom_node_text(aes(
        x = x,
        y = y,
        label = ifelse(Description %in% cluster_label_module, Description, NA)
      ),
      size = 3,
      repel = TRUE) +
      guides(fill = guide_legend(ncol = 1)) +
      scale_edge_width_continuous(range = c(0.1, 2)) +
      scale_size_continuous(range = c(1, 7)) +
      scale_fill_manual(values = database_color) +
      ggraph::theme_graph() +
      theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "right",
        legend.background = element_rect(fill = "transparent", color = NA),
      )
    
    plot
    
    extrafont::loadfonts()
    
    ggsave(
      plot,
      filename = file.path(path, "/combined_similarity_plot.pdf"),
      width = 9,
      height = 7
    )
  }

plot_module_bar <- function(module_result) {
  plot =
    module_result %>%
    dplyr::mutate(log.p = -log(p.adjust, 10)) %>%
    dplyr::mutate(Count = as.numeric(Count)) %>%
    dplyr::arrange(log.p) %>%
    dplyr::mutate(module_annotation = factor(module_annotation, levels = module_annotation)) %>%
    ggplot(aes(log.p, module_annotation)) +
    scale_y_discrete(
      labels = function(x)
        str_wrap(x, width = 50)
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
    geom_segment(aes(
      x = 0,
      y = module_annotation,
      xend = log.p,
      yend = module_annotation
    )) +
    geom_point(
      aes(size = Count),
      fill = "black",
      shape = 21,
      alpha = 1
    ) +
    scale_size_continuous(range = c(3, 7)) +
    theme_bw() +
    labs(y = "", x = "-log10(FDR adjusted P-values)") +
    geom_vline(xintercept = 0) +
    theme(panel.grid.minor = element_blank())
  
  plot
}

###overlap between differnt pathways
library(UpSetR)
library(ComplexHeatmap)

####ggraph to show the enriched pathways
plot_pathway_gene_network =
  function(module_result, marker_info) {
    library(tidygraph)
    library(ggraph)
    library(igraph)
    temp_data <-
      module_result %>%
      dplyr::arrange(p.adjust) %>%
      plyr::dlply(.variables = .(module_annotation)) %>%
      purrr::map(
        .f = function(x) {
          gene_id <- stringr::str_split(x$geneID, pattern = "/")[[1]]
          
          x <-
            data.frame(x %>% dplyr::select(-geneID), gene_id, stringsAsFactors = FALSE)
          
          gene_id1 =
            marker_info$variable_id[match(x$gene_id, marker_info$ensembl)]
          gene_id2 =
            marker_info$variable_id[match(x$gene_id, marker_info$uniprot)]
          gene_id3 =
            marker_info$variable_id[match(x$gene_id, marker_info$entrezid)]
          
          gene_id =
            data.frame(gene_id1, gene_id2, gene_id3) %>%
            apply(1, function(x) {
              x[!is.na(x)][1]
            })
          
          x$gene_id <-
            gene_id
          
          x <-
            x %>%
            dplyr::filter(!is.na(gene_id))
          
          x$gene_score <-
            marker_info$score[match(x$gene_id, marker_info$variable_id)]
          x$gene_fdr <-
            marker_info$fdr[match(x$gene_id, marker_info$variable_id)]
          x$recover_score <-
            marker_info$recover_score[match(x$gene_id, marker_info$variable_id)]
          x$Count <- nrow(x)
          
          x
        }
      ) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    edge_data <- temp_data
    node_data <- temp_data
    
    node_data1 <-
      node_data %>%
      dplyr::select(module_annotation, p.adjust, Count) %>%
      dplyr::rename(name = module_annotation,
                    fdr = p.adjust,
                    count = Count) %>%
      dplyr::mutate(score = NA) %>%
      dplyr::distinct(name, .keep_all = TRUE) %>%
      dplyr::mutate(class = "pathway") %>%
      dplyr::arrange(fdr) %>%
      dplyr::mutate(path_name = name)
    
    node_data2 <-
      node_data %>%
      dplyr::select(gene_id, gene_fdr, Count, gene_score) %>%
      dplyr::rename(
        name = gene_id,
        fdr = gene_fdr,
        score = gene_score,
        count = Count
      ) %>%
      dplyr::mutate(count = NA) %>%
      dplyr::distinct(name, .keep_all = TRUE) %>%
      dplyr::mutate(class = "gene") %>%
      dplyr::mutate(path_name = NA)
    
    node_data2$count <-
      purrr::map(node_data2$name, function(x) {
        sum(x == edge_data$gene_id)
      }) %>%
      unlist()
    
    node_data <- rbind(node_data1, node_data2)
    
    edge_data <-
      temp_data
    
    rownames(node_data) <-
      rownames(edge_data) <-
      NULL
    
    edge_data <-
      edge_data %>%
      dplyr::select(from = module_annotation, to = gene_id) %>%
      dplyr::mutate(path = from)
    
    node_data$symbol <-
      marker_info$symbol[match(node_data$name, marker_info$variable_id)]
    
    node_data <-
      node_data %>%
      mutate(symbol = case_when(is.na(symbol) ~ name, TRUE ~ symbol))
    
    total_graph <-
      tidygraph::tbl_graph(nodes = node_data,
                           edges = edge_data,
                           directed = FALSE) %>%
      dplyr::mutate(Degree = centrality_degree(mode = 'all'))
    
    g <- total_graph
    
    V(g)$type <- bipartite_mapping(g)$type
    
    coords <-
      ggraph::create_layout(g, layout = "bipartite") %>%
      dplyr::select(name, class, x, y)
    
    coords$index = 1:nrow(coords)
    
    coords$y[coords$y == 0] <- 1.6
    
    coords =
      coords %>%
      plyr::dlply(.variables = .(class)) %>%
      purrr::map(function(x) {
        x =
          x %>%
          dplyr::arrange(x)
        
        x$x = seq(0, 100, length.out = nrow(x))
        x
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame() %>%
      dplyr::arrange(index)
    
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
      ggraph::create_layout(
        graph = g,
        layout = "manual",
        x = coords$x,
        y = coords$y
        # node.position = coords
      )
    
    # RColorBrewer::display.brewer.all()
    path_col <-
      colorRampPalette(colors =
                         RColorBrewer::brewer.pal(n = 11, name = "Spectral"))(length(unique(edge_data$path)))
    
    names(path_col) <- unique(edge_data$path)
    
    plot <-
      ggraph(my_graph, layout = 'bipartite') +
      geom_edge_diagonal(
        strength = 1,
        aes(color = path),
        edge_width = 0.5,
        alpha = 1,
        show.legend = FALSE
      ) +
      geom_node_point(
        aes(fill = score, size = count),
        shape = 21,
        alpha = 1,
        show.legend = TRUE
      ) +
      geom_node_text(
        aes(
          x = x * 1.03,
          y = y * 1.03,
          hjust = ifelse(class == "pathway", "inward", 'outward'),
          size = ifelse(class == "pathway", 3, 2),
          angle = -((-node_angle(x, y) + 90) %% 180) + 90,
          label = symbol
        ),
        show.legend = FALSE
      ) +
      ggraph::scale_edge_color_manual(values = path_col) +
      scale_size_continuous(range = c(1, 8)) +
      scale_fill_gradient2(
        low = "blue",
        mid = "white",
        high = "red",
        midpoint = 0
      ) +
      ggraph::theme_graph() +
      theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "right",
        legend.background = element_rect(fill = "transparent", color = NA)
      )
    
    plot
    
  }

get_go_result_sim <-
  function(result,
           sim.cutoff = 0,
           measure_method = c("Wang", "Resnik", "Rel", "Jiang", "Lin", "TCSS")) {
    measure_method <-
      match.arg(measure_method)
    
    if (nrow(result) == 0) {
      return(data.frame(
        name1 = character(),
        name2 = character(),
        sim = numeric()
      ))
    }
    
    bp_sim_matrix <-
      simplifyEnrichment::GO_similarity(go_id = result$ID[result$ONTOLOGY == "BP"],
                                        ont = "BP",
                                        measure = measure_method) %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "name1") %>%
      tidyr::pivot_longer(cols = -name1,
                          names_to = "name2",
                          values_to = "sim") %>%
      dplyr::filter(name1 != name2) %>%
      dplyr::filter(sim > sim.cutoff)
    
    name <- apply(bp_sim_matrix, 1, function(x) {
      paste(sort(x[1:2]), collapse = "_")
    })
    
    bp_sim_matrix <-
      bp_sim_matrix %>%
      dplyr::mutate(name = name) %>%
      dplyr::arrange(name) %>%
      dplyr::distinct(name, .keep_all = TRUE) %>%
      dplyr::select(-name)
    
    mf_sim_matrix <-
      simplifyEnrichment::GO_similarity(go_id = result$ID[result$ONTOLOGY == "MF"],
                                        ont = "MF",
                                        measure = measure_method) %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "name1") %>%
      tidyr::pivot_longer(cols = -name1,
                          names_to = "name2",
                          values_to = "sim") %>%
      dplyr::filter(name1 != name2) %>%
      dplyr::filter(sim > sim.cutoff)
    
    name <- apply(mf_sim_matrix, 1, function(x) {
      paste(sort(x[1:2]), collapse = "_")
    })
    
    mf_sim_matrix <-
      mf_sim_matrix %>%
      dplyr::mutate(name = name) %>%
      dplyr::arrange(name) %>%
      dplyr::distinct(name, .keep_all = TRUE) %>%
      dplyr::select(-name)
    
    cc_sim_matrix <-
      simplifyEnrichment::GO_similarity(go_id = result$ID[result$ONTOLOGY == "CC"],
                                        ont = "CC",
                                        measure = measure_method) %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "name1") %>%
      tidyr::pivot_longer(cols = -name1,
                          names_to = "name2",
                          values_to = "sim") %>%
      dplyr::filter(name1 != name2) %>%
      dplyr::filter(sim > sim.cutoff)
    
    name <- apply(cc_sim_matrix, 1, function(x) {
      paste(sort(x[1:2]), collapse = "_")
    })
    
    cc_sim_matrix <-
      cc_sim_matrix %>%
      dplyr::mutate(name = name) %>%
      dplyr::arrange(name) %>%
      dplyr::distinct(name, .keep_all = TRUE) %>%
      dplyr::select(-name)
    
    sim_matrix <-
      rbind(bp_sim_matrix, mf_sim_matrix, cc_sim_matrix) %>%
      as.data.frame()
    sim_matrix
  }


show_matrix_cluster <-
  function(result,
           ont = NULL,
           measure = "Wang",
           remove_words,
           margin = 15,
           width = 14,
           height = 8,
           path = "GO_result",
           top = 20) {
    result <-
      result %>%
      dplyr::mutate(module_content_number = as.numeric(module_content_number))
    
    result <-
      result %>%
      dplyr::mutate(module_content_number = as.numeric(module_content_number)) %>%
      dplyr::filter(module_content_number > 1)
    
    if (nrow(result) == 0) {
      return(NULL)
    }
    
    if (!is.null(ont)) {
      mat <-
        simplifyEnrichment::GO_similarity(
          go_id = result %>%
            dplyr::filter(ONTOLOGY == ont) %>%
            dplyr::pull(node),
          measure = measure,
          ont = ont
        )
      
      if (is.null(nrow(mat))) {
        return(NULL)
      }
    } else{
      mat =
        simplifyEnrichment::term_similarity_from_Reactome(term_id = result %>%
                                                            dplyr::pull(node),
                                                          method = "jaccard")
    }
    
    if (!is.null(ont)) {
      cc_order <-
        result %>%
        dplyr::arrange(module) %>%
        dplyr::filter(ONTOLOGY == ont)
    } else{
      cc_order <-
        result %>%
        dplyr::arrange(module)
    }
    
    cc_order =
      cc_order[stringr::str_order(cc_order$module, numeric = TRUE), ]
    
    ####only remain the top 20 module 2 with smallest P values
    temp_p <-
      cc_order %>%
      plyr::dlply(.variables = .(Direction)) %>%
      purrr::map(
        .f = function(x) {
          x %>%
            plyr::dlply(.variables = .(module)) %>%
            purrr::map(function(x) {
              min(x$p.adjust)
            }) %>%
            unlist() %>%
            sort() %>%
            head(top)
        }
      )
    
    temp_p =
      temp_p[c("UP", 'DOWN')] %>%
      unlist()
    
    names(temp_p) =
      names(temp_p) %>%
      stringr::str_replace_all("UP\\.", "") %>%
      stringr::str_replace_all("DOWN\\.", "")
    
    
    cc_order =
      cc_order %>%
      dplyr::filter(module %in% names(temp_p))
    
    mat = mat[cc_order$node, cc_order$node]
    
    keywords <-
      cc_order %>%
      plyr::dlply(.variables = .(module)) %>%
      purrr::map(
        .f = function(x) {
          x <-
            x %>%
            dplyr::select(Description, p.adjust) %>%
            dplyr::mutate(p.adjust = -log(as.numeric(p.adjust), 10)) %>%
            plyr::dlply(.variables = .(Description)) %>%
            purrr::map(
              .f = function(y) {
                data.frame(
                  word = stringr::str_split(y$Description, " ")[[1]],
                  p.adjust = y$p.adjust
                )
              }
            ) %>%
            do.call(rbind, .) %>%
            as.data.frame() %>%
            dplyr::arrange(word)
          
          rownames(x) <- NULL
          
          x <-
            x %>%
            dplyr::group_by(word) %>%
            dplyr::summarise(freq = sum(as.numeric(p.adjust))) %>%
            dplyr::ungroup() %>%
            dplyr::filter(!word %in% remove_words)
          
        }
      )
    
    keywords =
      keywords[unique(cc_order$module)]
    
    
    name <-
      lapply(keywords, function(x) {
        !is.null(x)
      }) %>%
      unlist() %>% which() %>%
      names()
    
    keywords =
      keywords[names(keywords) %in% name]
    
    library(ComplexHeatmap)
    library(circlize)
    col_fun = colorRamp2(c(0, 1), c("white", "red"))
    
    module_col <-
      unique(result$module)[unique(result$module) != 'Other']
    
    module_col <-
      colorRampPalette(colors = RColorBrewer::brewer.pal(n = 10, name = "Spectral")[c(7, 2, 6, 10, 9, 1, 5, 8, 4, 3)])(n = length(module_col))
    
    names(module_col) <-
      unique(result$module)[unique(result$module) != 'Other']
    
    direction_col = module_col[unique(cc_order$module)]
    direction_col[names(direction_col) %in% cc_order$module[cc_order$Direction == "UP"]] =
      ggsci::pal_aaas()(n = 10)[2]
    
    direction_col[names(direction_col) %in% cc_order$module[cc_order$Direction == "DOWN"]] =
      ggsci::pal_aaas()(n = 10)[1]
    
    ha1 =
      rowAnnotation(
        foo = anno_block(
          gp = gpar(fill = direction_col),
          labels = NULL,
          labels_gp = gpar(col = "white", fontsize = 10),
          width = unit(5, "mm")
        ),
        foo2 = anno_block(
          gp = gpar(fill = unique(module_col[cc_order$module])),
          labels = unique(stringr::str_replace(cc_order$module, "Module ", "")),
          labels_gp = gpar(col = "white", fontsize = 10),
          width = unit(5, "mm")
        )
      )
    
    ha2 =
      HeatmapAnnotation(foo = anno_block(
        gp = gpar(fill = unique(module_col[cc_order$module])),
        labels = unique(stringr::str_replace(cc_order$module, "Module ", "")),
        labels_gp = gpar(col = "white", fontsize = 10),
        height = unit(5, "mm")
      ))
    
    ht = Heatmap(
      mat,
      cluster_columns = FALSE,
      cluster_rows = FALSE,
      show_row_names = FALSE,
      show_column_names = FALSE,
      name = "GO similarity",
      col = col_fun,
      right_annotation = ha1,
      top_annotation = ha2,
      border = FALSE,
      row_split = as.numeric(factor(
        cc_order$module, levels = unique(cc_order$module)
      )),
      column_split = as.numeric(factor(
        cc_order$module, levels = unique(cc_order$module)
      )),
      column_title = NULL,
      row_title = NULL,
      row_gap = unit(0.5, "mm"),
      column_gap = unit(0.5, "mm")
    )
    
    cl =
      factor(cc_order$module, levels = unique(cc_order$module))
    
    align_to = split(seq_len(nrow(mat)), cl)
    align_to = align_to[names(align_to) %in% names(keywords)]
    align_to
    
    fontsize_range = c(8, 16)
    
    gbl = lapply(names(align_to), function(nm) {
      kw = keywords[[nm]]$word
      freq = keywords[[nm]]$freq
      fontsize = scale_fontsize(freq, rg = c(1, max(10, freq)), fs = fontsize_range)
      word_cloud_grob(
        text = kw,
        fontsize = fontsize,
        col = colorRampPalette(colors = ggsci::pal_aaas()(n =
                                                            10))(length(kw))
      )
    })
    
    names(gbl) = names(align_to)
    
    # margin = unit(8, "pt")
    margin = unit(margin, "pt")
    gbl_h = lapply(gbl, function(x)
      convertHeight(grobHeight(x), "cm") + margin)
    
    gbl_h = do.call(unit.c, gbl_h)
    
    gbl_w = lapply(gbl, function(x)
      convertWidth(grobWidth(x), "cm"))
    
    gbl_w = do.call(unit.c, gbl_w)
    
    gbl_w = max(gbl_w) + margin
    
    panel_fun = function(index, nm) {
      # background
      grid.rect(gp = gpar(fill = "#DDDDDD", col = NA))
      # border
      grid.lines(
        c(0, 1, 1, 0),
        c(0, 0, 1, 1),
        gp = gpar(col = "#AAAAAA"),
        default.units = "npc"
      )
      gb = gbl[[nm]]
      # a viewport within the margins
      pushViewport(
        viewport(
          x = margin / 2,
          y = margin / 2,
          width = grobWidth(gb),
          height = grobHeight(gb),
          just = c("left", "bottom")
        )
      )
      grid.draw(gb)
      popViewport()
    }
    
    ht = ht + rowAnnotation(
      keywords = anno_link(
        align_to = align_to,
        which = "row",
        panel_fun = panel_fun,
        size = gbl_h,
        gap = unit(2, "mm"),
        width = gbl_w + unit(5, "mm"),
        # 5mm for the link
        link_gp = gpar(fill = "#DDDDDD", col = "#AAAAAA"),
        internal_line = FALSE
      )
    ) # you can set it to TRUE to see what happens
    
    ht
    
    if (!is.null(ont)) {
      height = case_when(ont == "BP" ~ 12, ont == "MF" ~ 5, ont == "CC" ~ 6)
    } else{
      height = height
    }
    
    ht
    
    pdf(file = file.path(path, paste(ont, "_sim_matrix.pdf", sep = "")),
        width = width,
        height = height)
    draw(ht, ht_gap = unit(2, "pt"))
    dev.off()
  }


get_jaccard_index_for_three_databases <-
  function(module_result_go,
           module_result_kegg,
           module_result_reactome,
           variable_info) {
    met_data <-
      rbind(module_result_go,
            module_result_kegg,
            module_result_reactome)
    
    if (nrow(met_data) == 0 | nrow(met_data) == 1) {
      return(data.frame(
        name1 = character(),
        name2 = character(),
        value = numeric()
      ))
    }
    
    temp_data <-
      met_data$geneID %>%
      stringr::str_split("/") %>%
      purrr::map(function(x) {
        if (stringr::str_detect(x[1], "ENSG")) {
          return(x)
        }
        
        if (stringr::str_detect(x[1], "[A-Za-z]")) {
          return(variable_info$ensembl[match(x, variable_info$uniprot)])
        }
        
        return(variable_info$ensembl[match(x, variable_info$entrezid)])
        
      })
    
    names(temp_data) =
      met_data$module
    
    ##calculate jaccard index
    jaccard_index =
      purrr::map(
        1:(length(temp_data) - 1),
        .f = function(idx) {
          purrr::map(
            temp_data[(idx + 1):length(temp_data)],
            .f = function(y) {
              length(intersect(temp_data[[idx]], y)) / length(union(temp_data[[idx]], y))
            }
          ) %>%
            unlist() %>%
            data.frame(value = .) %>%
            dplyr::mutate(name2 = names(temp_data)[(idx + 1):length(temp_data)]) %>%
            # tibble::rownames_to_column(var = "name2") %>%
            data.frame(name1 = names(temp_data)[idx], .) %>%
            dplyr::select(name1, name2, value)
        }
      ) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    return(jaccard_index)
  }

pathway_gene_level_heatmap <-
  function(module_result_all,
           module_result_go,
           module_result_kegg,
           module_result_reactome,
           subject_data,
           marker_info,
           sample_info) {
    library(plyr)
    library(ComplexHeatmap)
    library(circlize)
    col_fun = circlize::colorRamp2(seq(-3, 3, length.out = 9), rev(RColorBrewer::brewer.pal(n = 9, name = "BrBG")))
    ##heatmap
    temp_sample_info <-
      sample_info %>%
      dplyr::mutate(g_stage2 = g_stage)
    
    temp_sample_info$g_stage2[temp_sample_info$postdelivery] <- 50
    
    temp_sample_info <-
      temp_sample_info %>%
      dplyr::arrange(g_stage2)
    
    subject_data <-
      subject_data[, temp_sample_info$sample_id]
    
    pathway <-
      module_result_all %>%
      plyr::dlply(.variables = .(module_annotation)) %>%
      purrr::map(
        .f = function(x) {
          # cat(x$module_annotation[1], " ")
          name <- x$module_annotation[1]
          gene <- x$geneID %>%
            stringr::str_split("/") %>%
            `[[`(1)
          
          gene1 <-
            marker_info$variable_id[match(gene, marker_info$ensembl)]
          gene2 <-
            marker_info$variable_id[match(gene, marker_info$uniprot)]
          gene3 <-
            marker_info$variable_id[match(gene, marker_info$entrezid)]
          
          gene =
            data.frame(gene1, gene2, gene3) %>%
            dplyr::filter(!is.na(gene1) |
                            !is.na(gene2) | !is.na(gene3)) %>%
            apply(1, function(x) {
              x[!is.na(x)][1]
            })
          
          data <-
            subject_data[gene, , drop = FALSE] %>%
            as.data.frame()
          
          data <-
            data.frame(
              name,
              gene = gene,
              data,
              stringsAsFactors = FALSE,
              check.names = FALSE
            )
          data
        }
      ) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    rownames(pathway) <- NULL
    
    range(temp_sample_info$g_stage2)
    ga <- temp_sample_info$g_stage2
    ga[ga == 50] <- NA
    
    ha1 = HeatmapAnnotation(
      ga = ga,
      na_col = "Black",
      col = list(ga =
                   circlize::colorRamp2(
                     breaks = seq(
                       min(ga, na.rm = TRUE),
                       to = max(ga, na.rm = TRUE),
                       length.out = 11
                     ),
                     colors = rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral"))
                   )),
      annotation_name_side = c("left"),
      show_legend = FALSE
    )
    
    path_col <-
      colorRampPalette(colors =
                         RColorBrewer::brewer.pal(n = 11, name = "Spectral"))(length(unique(pathway$name)))
    
    names(path_col) <- unique(pathway$name)
    
    path_level <- unique(pathway$name)
    
    ha2 =
      rowAnnotation(foo2 = anno_block(
        gp = gpar(fill = path_col),
        labels = rev(path_level),
        labels_gp = gpar(col = "black", fontsize = 10),
        width = unit(2, "mm"),
        labels_rot = 0
      ))
    
    pathway_data <- pathway %>%
      dplyr::select(-c(name, gene))
    
    pathway_data <-
      pathway_data %>%
      dplyr::select(temp_sample_info$sample_id)
    
    range(pathway_data)
    
    if (abs(range(pathway_data)[1]) > range(pathway_data)[2]) {
      pathway_data[pathway_data < -range(pathway_data)[2]] <-
        -range(pathway_data)[2]
    } else{
      pathway_data[pathway_data > -range(pathway_data)[1]] <-
        -range(pathway_data)[1]
    }
    
    plot <-
      Heatmap(
        as.matrix(pathway_data),
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        border = TRUE,
        col = col_fun,
        name = "z-score",
        clustering_method_rows = "ward.D",
        top_annotation = ha1,
        left_annotation = ha2,
        row_split = as.numeric(factor(pathway$name, levels = unique(pathway$name))),
        row_title = NULL,
        column_split = c(rep(1, sum(!is.na(
          ga
        ))), rep(2, sum(is.na(
          ga
        )))),
        column_title = NULL
      )
    # plot
    plot <- as.ggplot(plot)
    
    plot
  }

##ggraph to show the enriched pathways
graph_pathway_heatmap <-
  function(module_result,
           marker_info,
           subject_data,
           show_column_names = TRUE) {
    library(plyr)
    temp_data <-
      module_result %>%
      dplyr::arrange(p.adjust) %>%
      plyr::dlply(.variables = .(module_annotation)) %>%
      purrr::map(
        .f = function(x) {
          # message(x$module_annotation)
          gene_id <-
            stringr::str_split(x$geneID, pattern = "/")[[1]]
          
          x <-
            data.frame(x %>% dplyr::select(-geneID), gene_id, stringsAsFactors = FALSE)
          # x$gene_id <-
          #   marker_info$gene_id[match(x$gene_id, marker_info$ensembl)]
          
          gene_id1 =
            marker_info$variable_id[match(x$gene_id, marker_info$ensembl)]
          gene_id2 =
            marker_info$variable_id[match(x$gene_id, marker_info$uniprot)]
          gene_id3 =
            marker_info$variable_id[match(x$gene_id, marker_info$entrezid)]
          
          gene_id =
            data.frame(gene_id1, gene_id2, gene_id3) %>%
            # dplyr::filter(!is.na(gene_id1) | !is.na(gene_id2) | !is.na(gene_id3)) %>%
            apply(1, function(x) {
              x[!is.na(x)][1]
            })
          
          x$gene_id <-
            gene_id
          
          x <-
            x %>%
            dplyr::filter(!is.na(gene_id))
          
          x$gene_score <-
            marker_info$score[match(x$gene_id, marker_info$variable_id)]
          x$gene_fdr <-
            marker_info$fdr[match(x$gene_id, marker_info$variable_id)]
          x$recover_score <-
            marker_info$recover_score[match(x$gene_id, marker_info$variable_id)]
          x$Count <- nrow(x)
          x
        }
      ) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    edge_data <- temp_data
    node_data <- temp_data
    
    node_data1 <-
      node_data %>%
      dplyr::select(module_annotation, p.adjust, Count) %>%
      dplyr::rename(name = module_annotation,
                    fdr = p.adjust,
                    count = Count) %>%
      dplyr::mutate(score = NA) %>%
      dplyr::distinct(name, .keep_all = TRUE) %>%
      dplyr::mutate(class = "pathway") %>%
      dplyr::arrange(fdr) %>%
      dplyr::mutate(path_name = name)
    
    node_data2 <-
      node_data %>%
      dplyr::select(gene_id, gene_fdr, Count, gene_score) %>%
      dplyr::rename(
        name = gene_id,
        fdr = gene_fdr,
        score = gene_score,
        count = Count
      ) %>%
      dplyr::mutate(count = NA) %>%
      dplyr::distinct(name, .keep_all = TRUE) %>%
      dplyr::mutate(class = "gene") %>%
      dplyr::mutate(path_name = NA)
    
    node_data2$count <-
      purrr::map(node_data2$name, function(x) {
        sum(x == edge_data$gene_id)
      }) %>%
      unlist()
    
    node_data <- rbind(node_data1, node_data2)
    
    edge_data <-
      temp_data
    
    edge_data <-
      edge_data %>%
      dplyr::select(from = module_annotation, to = gene_id) %>%
      dplyr::mutate(path = from)
    
    node_data$symbol <-
      marker_info$symbol[match(node_data$name, marker_info$variable_id)]
    
    node_data <-
      node_data %>%
      mutate(symbol = case_when(is.na(symbol) ~ name, TRUE ~ symbol))
    
    total_graph <-
      tidygraph::tbl_graph(nodes = node_data,
                           edges = edge_data,
                           directed = FALSE) %>%
      dplyr::mutate(Degree = centrality_degree(mode = 'all'))
    
    library(igraph)
    
    g <- total_graph
    
    V(g)$type <- bipartite_mapping(g)$type
    
    coords <-
      create_layout(g, layout = "bipartite") %>%
      dplyr::select(name, class, x, y)
    
    coords$index = 1:nrow(coords)
    
    coords =
      coords %>%
      plyr::dlply(.variables = .(class)) %>%
      purrr::map(
        .f = function(x) {
          x =
            x %>%
            dplyr::arrange(x)
          x$x =
            seq(from = 1,
                to = 100,
                length.out = nrow(x))
          
          x %>%
            dplyr::arrange(index)
          
        }
      ) %>%
      do.call(rbind, .) %>%
      as.data.frame() %>%
      dplyr::arrange(index)
    
    coords <-
      coords %>%
      dplyr::mutate(x1 = y, y1 = x) %>%
      dplyr::select(-c(x, y)) %>%
      dplyr::rename(x = x1, y = y1) %>%
      dplyr::left_join(marker_info[, c("variable_id", "symbol")], by = c("name" = "variable_id"))
    
    coords$x[coords$x == 0] <- 2
    
    my_graph <-
      create_layout(
        graph = g,
        layout = "manual",
        x = coords$x,
        y = coords$y
        # node.position = coords
      )
    
    library(ggraph)
    
    # RColorBrewer::display.brewer.all()
    path_col <-
      colorRampPalette(colors =
                         RColorBrewer::brewer.pal(n = 11, name = "Spectral"))(length(unique(edge_data$path)))
    
    names(path_col) <- unique(edge_data$path)
    
    plot1 <-
      ggraph(my_graph, layout = 'bipartite') +
      geom_edge_diagonal(
        strength = 1,
        aes(color = path),
        edge_width = 0.5,
        alpha = 1,
        show.legend = FALSE
      ) +
      geom_node_point(
        aes(fill = score, size = count),
        shape = 21,
        alpha = 1,
        show.legend = FALSE
      ) +
      geom_node_text(aes(
        x = x,
        y = y,
        hjust = 1,
        size = 3,
        label = ifelse(class == "pathway", symbol, NA)
      ), show.legend = FALSE) +
      ggraph::scale_edge_color_manual(values = path_col) +
      scale_size_continuous(range = c(1, 8)) +
      scale_fill_gradient2(
        low = "blue",
        mid = "white",
        high = "red",
        midpoint = 0
      ) +
      ggraph::theme_graph() +
      theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "right",
        legend.background = element_rect(fill = "transparent", color = NA)
      )
    
    # plot1
    
    ###important genes heatmap
    col_fun <-
      colorRamp2(seq(-3, 3, length.out = 9), rev(RColorBrewer::brewer.pal(n = 9, name = "BrBG")))
    ###quantitative pathways
    ##heatmap
    library(plyr)
    
    pathway <-
      module_result %>%
      plyr::dlply(.variables = .(module_annotation)) %>%
      purrr::map(
        .f = function(x) {
          name <- x$module_annotation[1]
          gene <- x$geneID %>%
            stringr::str_split("/") %>%
            `[[`(1)
          
          gene1 =
            marker_info$variable_id[match(gene, marker_info$ensembl)]
          gene2 =
            marker_info$variable_id[match(gene, marker_info$entrezid)]
          gene3 =
            marker_info$variable_id[match(gene, marker_info$uniprot)]
          
          gene =
            data.frame(gene1, gene2, gene3) %>%
            apply(1, function(x) {
              x[!is.na(x)][1]
            })
          
          gene <- gene[!is.na(gene)]
          
          data <- subject_data[gene, ] %>%
            as.data.frame()
          
          data <-
            data.frame(
              name,
              gene = gene,
              data,
              stringsAsFactors = FALSE,
              check.names = FALSE
            )
          
          data
        }
      ) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    pathway <-
      pathway %>%
      dplyr::distinct(gene, .keep_all = TRUE)
    
    pathway =
      coords %>%
      dplyr::filter(class != "pathway") %>%
      dplyr::arrange(desc(y)) %>%
      dplyr::pull(name) %>%
      match(., pathway$gene) %>%
      `[`(pathway, ., )
    rownames(pathway) <- NULL
    
    rownames(pathway) <-
      marker_info$symbol[match(pathway$gene, marker_info$variable_id)]
    
    pathway_data <- pathway %>%
      dplyr::select(-c(name, gene))
    
    # range(pathway_data)
    
    if (abs(range(pathway_data)[1]) > range(pathway_data)[2]) {
      pathway_data[pathway_data < -range(pathway_data)[2]] <-
        -range(pathway_data)[2]
    } else{
      pathway_data[pathway_data > -range(pathway_data)[1]] <-
        -range(pathway_data)[1]
    }
    
    score <-
      marker_info$score[match(pathway$gene, marker_info$variable_id)]
    recover_score <-
      marker_info$recover_score[match(pathway$gene, marker_info$variable_id)]
    recover_score[is.infinite(recover_score) &
                    recover_score > 0] <- 2
    recover_score[is.infinite(recover_score) &
                    recover_score < 0] <- -2
    
    # lm_reg <-
    #   lm(formula = y ~ x,
    #      data = data.frame(y = font_size, x = word_count_breaks))
    #
    # text <- list(
    #   data.frame(text = word) %>%
    #     dplyr::group_by(text) %>%
    #     dplyr::summarise(fontsize = dplyr::n()) %>%
    #     as.data.frame()
    # ) %>%
    #   purrr::map(function(z) {
    #     z$fontsize <-
    #       predict(lm_reg, newdata = data.frame(x = z$fontsize))
    #     z
    #   })
    #
    # names(text) = "a"
    
    ha_right1 <-
      rowAnnotation(score = anno_points(x = score),
                    recover_score = anno_points(x = recover_score))
    
    plot2 <-
      Heatmap(
        as.matrix(pathway_data),
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        show_row_names = FALSE,
        show_column_names = TRUE,
        border = TRUE,
        col = col_fun,
        name = "z-score",
        clustering_method_rows = "ward.D",
        column_names_rot = 45,
        column_split = c(rep(1, sum(
          colnames(pathway_data) != "PP"
        )), rep(2, 1)),
        column_title = NULL,
        right_annotation = ha_right1
      )
    
    library(ggplotify)
    
    plot2 <- as.ggplot(plot2)
    
    library(patchwork)
    
    plot =
      plot1 + plot2 + patchwork::plot_layout(widths = c(1, 4))
    plot
    
  }


##ggraph to show the enriched pathways
pathway_heatmap <-
  function(module_result,
           marker_info,
           subject_data,
           show_column_names = TRUE,
           row_names_gp = gpar(cex = 0.5)) {
    library(plyr)
    ###important genes heatmap
    col_fun <-
      colorRamp2(seq(-3, 3, length.out = 9), rev(RColorBrewer::brewer.pal(n = 9, name = "BrBG")))
    ###quantitative pathways
    ##heatmap
    library(plyr)
    
    pathway_data <-
      module_result %>%
      plyr::dlply(.variables = .(module_annotation)) %>%
      purrr::map(
        .f = function(x) {
          name <- x$module_annotation[1]
          gene <- x$geneID %>%
            stringr::str_split("/") %>%
            `[[`(1)
          
          gene1 =
            marker_info$variable_id[match(gene, marker_info$ensembl)]
          gene2 =
            marker_info$variable_id[match(gene, marker_info$entrezid)]
          gene3 =
            marker_info$variable_id[match(gene, marker_info$uniprot)]
          
          gene =
            data.frame(gene1, gene2, gene3) %>%
            apply(1, function(x) {
              x[!is.na(x)][1]
            })
          
          gene <- gene[!is.na(gene)]
          
          data <- subject_data[gene, ] %>%
            as.data.frame()
          
          data <-
            data.frame(
              name,
              gene = gene,
              data,
              stringsAsFactors = FALSE,
              check.names = FALSE
            )
          
          data %>%
            dplyr::select(-c(name, gene)) %>%
            colMeans()
        }
      ) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    if (abs(range(pathway_data)[1]) > range(pathway_data)[2]) {
      pathway_data[pathway_data < -range(pathway_data)[2]] <-
        -range(pathway_data)[2]
    } else{
      pathway_data[pathway_data > -range(pathway_data)[1]] <-
        -range(pathway_data)[1]
    }
    
    ha_right1 <-
      rowAnnotation(
        count = anno_points(x = as.numeric(module_result$Count), ylim = c(0, 450)),
        "-log(p.adjust)" = anno_points(
          x = -log(module_result$p.adjust, 10),
          ylim = c(0, 30)
        )
      )
    
    plot <-
      Heatmap(
        as.matrix(pathway_data),
        cluster_columns = FALSE,
        cluster_rows = TRUE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_names_gp = row_names_gp,
        border = TRUE,
        col = col_fun,
        name = "z-score",
        clustering_method_rows = "ward.D",
        column_names_rot = 45,
        right_annotation = ha_right1,
        column_split = c(rep(1, sum(
          colnames(pathway_data) != "PP"
        )), rep(2, 1)),
        column_title = NULL
      )
    
    library(ggplotify)
    
    # plot <- as.ggplot(plot)
    
    library(patchwork)
    
    plot
  }


transcriptome_gene_ga_heatmap2 <-
  function(result_all, gene_marker, subject_data2) {
    path_level =
      unique(result_all$module_annotation)
    col_fun = colorRamp2(c(-3, 0, 3), c("#4292C6", "white", "red"))
    library(plyr)
    pathway <-
      result_all %>%
      plyr::dlply(.variables = .(module_annotation)) %>%
      purrr::map(
        .f = function(x) {
          name <- x$module_annotation[1]
          gene <- x$geneID %>%
            stringr::str_split("/") %>%
            `[[`(1)
          
          gene1 <-
            gene_marker$gene_id[match(gene, gene_marker$ensembl)]
          gene2 <-
            gene_marker$gene_id[match(gene, gene_marker$uniprot)]
          gene3 <-
            gene_marker$gene_id[match(gene, gene_marker$entrezid)]
          
          gene =
            data.frame(gene1, gene2, gene3) %>%
            apply(1, function(x) {
              x[!is.na(x)][1]
            })
          
          data <- subject_data_mean[gene, ] %>%
            apply(2, mean)
          data
        }
      ) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    library(ComplexHeatmap)
    
    path_col <-
      colorRampPalette(colors =
                         RColorBrewer::brewer.pal(n = 11, name = "Spectral"))(length(path_level))
    
    names(path_col) <- unique(path_level)
    
    ha2 = rowAnnotation(
      pathway = factor(rownames(pathway), levels = rownames(pathway)),
      col = list(pathway = path_col),
      show_legend = FALSE,
      gp = gpar(col = "black")
      # annotation_name_side = c("right")
    )
    
    pathway_data <- pathway %>%
      apply(1, function(x) {
        (x - mean(x)) / sd(x)
      }) %>%
      t()
    
    if (abs(range(pathway_data)[1]) > range(pathway_data)[2]) {
      pathway_data[pathway_data < -range(pathway_data)[2]] <-
        -range(pathway_data)[2]
    } else{
      pathway_data[pathway_data > -range(pathway_data)[1]] <-
        -range(pathway_data)[1]
    }
    
    ###calculate p-values for each pathway
    temp_p_value <-
      result_all %>%
      plyr::dlply(.variables = .(module_annotation))  %>%
      purrr::map(
        .f = function(x) {
          name <- x$module_annotation[1]
          gene <- x$geneID %>%
            stringr::str_split("/") %>%
            `[[`(1)
          
          gene1 <-
            gene_marker$gene_id[match(gene, gene_marker$ensembl)]
          gene2 <-
            gene_marker$gene_id[match(gene, gene_marker$uniprot)]
          gene3 <-
            gene_marker$gene_id[match(gene, gene_marker$entrezid)]
          
          gene =
            data.frame(gene1, gene2, gene3) %>%
            apply(1, function(x) {
              x[!is.na(x)][1]
            })
          
          temp <-
            lapply(subject_data2, function(y) {
              y[, gene] %>%
                apply(1, mean)
            })
          
          temp_p <-
            c(1, lapply(temp[-c(1, length(temp))], function(z) {
              wilcox.test(z, temp[[1]])$p.value
            }) %>%
              unlist())
          
          temp_p <- c(temp_p, wilcox.test(temp[[length(temp)]], temp[[length(temp) - 1]])$p.value)
          
          temp_p %>% p.adjust(method = "fdr")
        }
      ) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    temp_p_value <-
      temp_p_value[match(rev(path_level), rownames(temp_p_value)), ]
    
    temp_p_value <-
      apply(temp_p_value, 2, function(x) {
        case_when(x > 0.05 ~ "",
                  x <= 0.05 & x > 0.01 ~ "*",
                  x <= 0.01 & x > 0.001 ~ "**",
                  x <= 0.001 ~ "***")
      })
    
    plot <-
      Heatmap(
        as.matrix(pathway_data),
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        border = TRUE,
        col = col_fun,
        column_title_gp = gpar(angle = 45),
        name = "z-score",
        clustering_method_rows = "ward.D",
        rect_gp = gpar(col = "white"),
        right_annotation = ha2,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(label = temp_p_value[i, j], x, y, gp = gpar(fontsize = 10))
        }
      )
    
    library(ggplotify)
    
    plot <- as.ggplot(plot)
    
    plot
  }



transcriptome_boxplot_for_pathway <-
  function(module_result,
           path,
           subject_data,
           marker_info) {
    data <- module_result
    
    data %>%
      plyr::dlply(.variables = .(module_annotation))  %>%
      purrr::map(
        .f = function(x) {
          name <- x$module_annotation[1] %>%
            stringr::str_replace_all("/", "_")
          cat(name, " ")
          gene <- x$geneID %>%
            stringr::str_split("/") %>%
            `[[`(1)
          
          gene1 <-
            marker_info$variable_id[match(gene, marker_info$ensembl)]
          gene2 <-
            marker_info$variable_id[match(gene, marker_info$uniprot)]
          gene3 <-
            marker_info$variable_id[match(gene, marker_info$entrezid)]
          
          
          gene =
            data.frame(gene1, gene2, gene3) %>%
            apply(1, function(x) {
              x[!is.na(x)][1]
            })
          
          temp <-
            subject_data[gene, ] %>%
            tibble::rownames_to_column(var = "gene_id") %>%
            tidyr::pivot_longer(cols = -gene_id, names_to = 'ga')
          
          temp %>%
            dplyr::mutate(ga = factor(ga, levels = unique(ga))) %>%
            ggplot(aes(ga, value)) +
            geom_violin() +
            geom_boxplot(
              fill = NA,
              width = 0.3,
              colour = ggsci::pal_aaas()(n = 10)[2],
              outlier.shape = NA
            ) +
            geom_jitter(shape = 21, fill = "black")
          
          
          plot <-
            temp %>%
            dplyr::mutate(ga = factor(ga, levels = unique(ga))) %>%
            ggplot(aes(ga, x)) +
            geom_violin() +
            geom_boxplot(
              fill = NA,
              width = 0.3,
              colour = ggsci::pal_aaas()(n = 10)[2],
              outlier.shape = NA
            ) +
            geom_jitter(shape = 21, fill = "black") +
            geom_smooth(method = "loess", se = TRUE, aes(group = 1)) +
            labs(x = "", y = "Scaled intensity") +
            theme_bw() +
            theme(
              axis.title = element_text(size = 13),
              axis.text = element_text(size = 12),
              axis.text.x = element_text(
                size = 12,
                angle = 45,
                hjust = 1,
                vjust = 1
              ),
              legend.title = element_text(size = 13),
              legend.text = element_text(size = 12),
              legend.position = c(0, 1),
              legend.justification = c(0, 1),
              panel.background = element_rect(fill = "transparent", color = NA),
              plot.background = element_rect(fill = "transparent", color = NA),
              legend.background = element_rect(fill = "transparent", color = NA),
              strip.background = element_rect(fill = "#0099B47F"),
              strip.text = element_text(color = "white", size = 15),
              panel.grid.minor = element_blank()
            )
          
          ggsave(
            plot,
            filename = file.path(path, paste(name, ".pdf", sep = "")),
            width = 8,
            height = 7
          )
          
        }
      )
  }






module_pathway_network_plot <- function(module_result = module_result_all %>% head(3),
                                        result_with_module = result_with_module,
                                        enrichment_go = enrichment_go,
                                        enrichment_kegg = enrichment_kegg,
                                        enrichment_reactome = enrichment_reactome,
                                        marker_info = marker,
                                        including_molecule = TRUE,
                                        circle = FALSE,
                                        node_color = c(
                                          "Function_module" = "#F05C3BFF",
                                          "Module" = "#46732EFF",
                                          "Pathway" = "#197EC0FF",
                                          "Molecule" = "#3B4992FF"
                                        ),
                                        text_size = 3,
                                        arrange_position = TRUE,
                                        position_ratio = 0.95) {
  edge_color <-
    c(
      "Function_module-Module" = unname(node_color["Function_module"]),
      "Module-Pathway" = unname(node_color["Module"]),
      "Pathway-Molecule" = unname(node_color["Pathway"])
    )
  
  ####function_module vs module
  edge_data1 <-
    data.frame(from = module_result$module, to = module_result$module_content) %>%
    apply(1, function(x) {
      data.frame(from = as.character(x[1]),
                 to = as.character(stringr::str_split(x[2], ";")[[1]]))
    }) %>%
    do.call(rbind, .) %>%
    as.data.frame() %>%
    dplyr::mutate(class = "Function_module-Module")
  
  node_data1 <-
    module_result %>%
    dplyr::select(node = module, annotation = module_annotation, p.adjust, Count) %>%
    dplyr::mutate(Count = as.numeric(Count), class = "Function_module")
  
  ###module vs pathway
  edge_data2 <-
    result_with_module %>%
    dplyr::select(from = node, to = pathway_id) %>%
    apply(1, function(x) {
      data.frame(from = as.character(x[1]),
                 to = as.character(stringr::str_split(x[2], ";")[[1]]))
    }) %>%
    do.call(rbind, .) %>%
    as.data.frame() %>%
    dplyr::mutate(class = "Module-Pathway") %>%
    dplyr::filter(from %in% edge_data1$to)
  
  node_data2 <-
    result_with_module %>%
    dplyr::select(node, annotation = module_annotation, p.adjust, Count, database) %>%
    dplyr::mutate(Count = as.numeric(Count), class = "Module") %>%
    dplyr::filter(node %in% edge_data2$from)
  
  
  #####pathway vs molecule
  edge_data3_go <-
    enrichment_go@result %>%
    dplyr::filter(p.adjust < 0.05 & Count > 5) %>%
    dplyr::filter(ONTOLOGY != "CC") %>%
    dplyr::select(from = ID, to = geneID) %>%
    apply(1, function(x) {
      data.frame(from = as.character(x[1]),
                 to = as.character(stringr::str_split(x[2], "/")[[1]]))
    }) %>%
    do.call(rbind, .) %>%
    as.data.frame()
  
  if (nrow(edge_data3_go) > 0) {
    edge_data3_go <-
      edge_data3_go %>%
      dplyr::filter(from %in% edge_data2$to)
  }
  
  node_data3_go <-
    enrichment_go@result %>%
    dplyr::filter(p.adjust < 0.05 & Count > 5) %>%
    dplyr::filter(ONTOLOGY != "CC") %>%
    dplyr::select(node = ID, annotation = Description, p.adjust, Count) %>%
    dplyr::mutate(database = "GO")
  
  if (nrow(node_data3_go) > 0) {
    node_data3_go <-
      node_data3_go %>%
      dplyr::filter(node %in% edge_data3_go$from)
  }
  
  edge_data3_kegg <-
    enrichment_kegg@result %>%
    dplyr::filter(p.adjust < 0.05 & Count > 5) %>%
    dplyr::select(from = ID, to = geneID) %>%
    apply(1, function(x) {
      data.frame(from = as.character(x[1]),
                 to = as.character(stringr::str_split(x[2], "/")[[1]]))
    }) %>%
    do.call(rbind, .) %>%
    as.data.frame()
  
  if (nrow(edge_data3_kegg) > 0) {
    edge_data3_kegg <-
      edge_data3_kegg %>%
      dplyr::filter(from %in% edge_data2$to)
  }
  
  node_data3_kegg <-
    enrichment_kegg@result %>%
    dplyr::filter(p.adjust < 0.05 & Count > 5) %>%
    dplyr::select(node = ID, annotation = Description, p.adjust, Count) %>%
    dplyr::mutate(database = "KEGG")
  
  if (nrow(node_data3_kegg) > 0) {
    node_data3_kegg <-
      node_data3_kegg %>%
      dplyr::filter(node %in% edge_data3_kegg$from)
  }
  
  edge_data3_reactome <-
    enrichment_reactome@result %>%
    dplyr::filter(p.adjust < 0.05 & Count > 5) %>%
    dplyr::select(from = ID, to = geneID) %>%
    apply(1, function(x) {
      data.frame(from = as.character(x[1]),
                 to = as.character(stringr::str_split(x[2], "/")[[1]]))
    }) %>%
    do.call(rbind, .) %>%
    as.data.frame()
  
  if (nrow(edge_data3_reactome) > 0) {
    edge_data3_reactome <-
      edge_data3_reactome %>%
      dplyr::filter(from %in% edge_data2$to)
  }
  
  node_data3_reactome <-
    enrichment_reactome@result %>%
    dplyr::filter(p.adjust < 0.05 & Count > 5) %>%
    dplyr::select(node = ID, annotation = Description, p.adjust, Count) %>%
    dplyr::mutate(database = "Reactome")
  
  if (nrow(node_data3_reactome) > 0) {
    node_data3_reactome <-
      node_data3_reactome %>%
      dplyr::filter(node %in% edge_data3_reactome$from)
  }
  
  node_data3_molecule <-
    marker_info %>%
    dplyr::select(
      node = ensembl,
      uniprot,
      symbol,
      entrezid,
      annotation = genename,
      p.adjust = fdr,
      SAM_score = score
    ) %>%
    dplyr::mutate(Count = 1, class = "Molecule")
  
  edge_data3_kegg$to <-
    marker_info$ensembl[match(edge_data3_kegg$to, marker_info$uniprot)]
  
  edge_data3_reactome$to <-
    marker_info$ensembl[match(edge_data3_reactome$to, marker_info$entrezid)]
  
  edge_data3 <-
    rbind(edge_data3_go, edge_data3_kegg, edge_data3_reactome) %>%
    dplyr::filter(!is.na(to)) %>%
    dplyr::distinct(from, to, .keep_all = TRUE) %>%
    dplyr::mutate(class = "Pathway-Molecule")
  
  node_data3_pathway <-
    rbind(node_data3_go, node_data3_kegg, node_data3_reactome) %>%
    dplyr::mutate(class = "Pathway")
  
  node_data3_molecule <-
    node_data3_molecule %>%
    dplyr::filter(node %in% edge_data3$to)
  
  node_data3 <-
    node_data3_pathway %>%
    dplyr::full_join(node_data3_molecule, by = intersect(
      colnames(node_data3_pathway),
      colnames(node_data3_molecule)
    ))
  edge_data <-
    rbind(edge_data1, edge_data2, edge_data3)
  
  node_data <-
    node_data1 %>%
    dplyr::full_join(node_data2, by = intersect(colnames(.), colnames(node_data2))) %>%
    dplyr::full_join(node_data3, by = intersect(colnames(.), colnames(node_data3)))
  
  node_data <-
    node_data %>%
    dplyr::filter(node %in% c(edge_data$from, edge_data$to))
  
  edge_data <-
    edge_data %>%
    dplyr::filter(from %in% node_data$node & to %in% node_data$node)
  
  # table(node_data$class)
  #
  # sum(!edge_data$from %in% node_data$node)
  # sum(!edge_data$to %in% node_data$node)
  # sum(!node_data$node %in% c(edge_data$from,
  #                            edge_data$to))
  
  if (!including_molecule) {
    temp_node_data <-
      node_data %>%
      dplyr::filter(class != "Molecule")
    temp_edge_data <-
      edge_data %>%
      dplyr::filter(class != "Pathway-Molecule")
    
    total_graph <-
      tidygraph::tbl_graph(nodes = temp_node_data,
                           edges = temp_edge_data,
                           directed = FALSE)
    
    g <- total_graph
    V(g)$type <- bipartite_mapping(g)$type
    coords <-
      ggraph::create_layout(g, layout = "bipartite")
    coords$index = 1:nrow(coords)
    
    coords$x <-
      coords$x + 1
    
    if (circle) {
      coords$y[coords$class == "Function_module"] <- 0
      coords$y[coords$class == "Module"] <- 1
      coords$y[coords$class == "Pathway"] <- 2
      coords$y[coords$class == "Molecule"] <- 3
      
      if (arrange_position) {
        coords <-
          arrange_coords(coords = coords, ratio = position_ratio)
      }
      
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
        ggraph::create_layout(
          graph = g,
          layout = "manual",
          x = coords$x,
          y = coords$y
          # node.position = coords
        )
      
      plot <-
        ggraph(my_graph, layout = 'bipartite') +
        geom_edge_diagonal(
          strength = 1,
          aes(color = class),
          edge_width = 0.5,
          alpha = 1,
          show.legend = FALSE
        ) +
        geom_node_point(
          aes(fill = class, size = Count),
          shape = 21,
          alpha = 1,
          show.legend = TRUE
        ) +
        geom_node_text(
          aes(
            x = x * 1.03,
            y = y * 1.03,
            hjust = ifelse(class == "Pathway", "outward", 'inward'),
            angle = -((-node_angle(x, y) + 90) %% 180) + 90,
            label = annotation
          ),
          size = text_size,
          show.legend = FALSE
        ) +
        scale_fill_manual(values = node_color) +
        ggraph::scale_edge_color_manual(values = edge_color) +
        scale_size_continuous(range = c(1, 8)) +
        ggraph::theme_graph() +
        theme(
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent", color = NA),
          legend.position = "right",
          legend.background = element_rect(fill = "transparent", color = NA)
        )
      
    } else{
      coords$y[coords$class == "Function_module"] <- 3
      coords$y[coords$class == "Module"] <- 2
      coords$y[coords$class == "Pathway"] <- 1
      coords$y[coords$class == "Molecule"] <- 0
      
      if (arrange_position) {
        coords <-
          arrange_coords(coords = coords, ratio = position_ratio)
      }
      
      my_graph <-
        ggraph::create_layout(
          graph = g,
          layout = "manual",
          x = coords$x,
          y = coords$y
          # node.position = coords
        )
      
      plot <-
        ggraph(my_graph, layout = 'bipartite') +
        geom_edge_diagonal(
          strength = 1,
          aes(color = class),
          edge_width = 0.5,
          alpha = 1,
          show.legend = FALSE
        ) +
        geom_node_point(
          aes(fill = class, size = Count),
          shape = 21,
          alpha = 1,
          show.legend = TRUE
        ) +
        geom_node_text(
          aes(x = x, y = y, label = annotation),
          hjust = 1,
          angle = 90,
          size = text_size,
          show.legend = FALSE
        ) +
        scale_fill_manual(values = node_color) +
        ggraph::scale_edge_color_manual(values = edge_color) +
        scale_size_continuous(range = c(1, 8)) +
        ggraph::theme_graph() +
        theme(
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent", color = NA),
          legend.position = "right",
          legend.background = element_rect(fill = "transparent", color = NA)
        )
      
    }
  } else{
    total_graph <-
      tidygraph::tbl_graph(nodes = node_data,
                           edges = edge_data,
                           directed = FALSE)
    
    g <- total_graph
    V(g)$type <- bipartite_mapping(g)$type
    coords <-
      ggraph::create_layout(g, layout = "bipartite")
    coords$index = 1:nrow(coords)
    
    if (circle) {
      coords$y[coords$class == "Function_module"] <- 0
      coords$y[coords$class == "Module"] <- 1
      coords$y[coords$class == "Pathway"] <- 2
      coords$y[coords$class == "Molecule"] <- 3
      
      if (arrange_position) {
        coords <-
          arrange_coords(coords = coords, ratio = position_ratio)
      }
      
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
        ggraph::create_layout(
          graph = g,
          layout = "manual",
          x = coords$x,
          y = coords$y
          # node.position = coords
        )
      
      plot <-
        ggraph(my_graph, layout = 'bipartite') +
        geom_edge_diagonal(
          strength = 1,
          aes(color = class),
          edge_width = 0.5,
          alpha = 1,
          show.legend = FALSE
        ) +
        geom_node_point(
          aes(fill = class, size = Count),
          shape = 21,
          alpha = 1,
          show.legend = TRUE
        ) +
        geom_node_text(
          aes(
            x = x * 1.03,
            y = y * 1.03,
            hjust = ifelse(class == "Pathway", "outward", 'inward'),
            angle = -((-node_angle(x, y) + 90) %% 180) + 90,
            label = ifelse(class == "Molecule", NA, annotation)
          ),
          size = text_size,
          show.legend = FALSE
        ) +
        scale_fill_manual(values = node_color) +
        ggraph::scale_edge_color_manual(values = edge_color) +
        scale_size_continuous(range = c(1, 8)) +
        ggraph::theme_graph() +
        theme(
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent", color = NA),
          legend.position = "right",
          legend.background = element_rect(fill = "transparent", color = NA)
        )
      
    } else{
      coords$y[coords$class == "Function_module"] <- 3
      coords$y[coords$class == "Module"] <- 2
      coords$y[coords$class == "Pathway"] <- 1
      coords$y[coords$class == "Molecule"] <- 0
      
      if (arrange_position) {
        coords <-
          arrange_coords(coords = coords, ratio = position_ratio)
      }
      
      my_graph <-
        ggraph::create_layout(
          graph = g,
          layout = "manual",
          x = coords$x,
          y = coords$y
          # node.position = coords
        )
      
      plot <-
        ggraph(my_graph, layout = 'bipartite') +
        geom_edge_diagonal(
          strength = 1,
          aes(color = class),
          edge_width = 0.5,
          alpha = 1,
          show.legend = FALSE
        ) +
        geom_node_point(
          aes(fill = class, size = Count),
          shape = 21,
          alpha = 1,
          show.legend = TRUE
        ) +
        geom_node_text(
          aes(
            x = x,
            y = y,
            label = ifelse(class == "Molecule", NA, annotation)
          ),
          hjust = 1,
          angle = 90,
          size = text_size,
          show.legend = FALSE
        ) +
        scale_fill_manual(values = node_color) +
        ggraph::scale_edge_color_manual(values = edge_color) +
        scale_size_continuous(range = c(1, 8)) +
        ggraph::theme_graph() +
        theme(
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent", color = NA),
          legend.position = "right",
          legend.background = element_rect(fill = "transparent", color = NA)
        )
      
    }
    
    
  }
  plot
}


arrange_coords <- function(coords, ratio = 0.95) {
  coords <-
    coords %>%
    plyr::dlply(.variables = .(class)) %>%
    purrr::map(function(x) {
      x <-
        x %>%
        dplyr::arrange(x)
      if (length(unique(coords$class)) == 4) {
        if (x$class[1] == "Molecule") {
          min_times <- 1
          max_times <- 1
        }
        
        if (x$class[1] == "Pathway") {
          min_times <- 1.1
          max_times <- ratio
        }
        
        if (x$class[1] == "Module") {
          min_times <- 1.1^2
          max_times <- ratio^2
        }
        
        if (x$class[1] == "Function_module") {
          min_times <- 1.1^3
          max_times <- ratio^3
        }
        
      }
      
      if (length(unique(coords$class)) == 3) {
        if (x$class[1] == "Pathway") {
          min_times <- 1
          max_times <- 1
        }
        
        if (x$class[1] == "Module") {
          min_times <- 1.1
          max_times <- ratio
        }
        
        if (x$class[1] == "Function_module") {
          min_times <- 1.1^2
          max_times <- ratio^2
        }
        
      }
      
      x$x <-
        seq(
          from = max(coords$x) - max_times * max(coords$x),
          to = max_times * max(coords$x),
          length.out = nrow(x)
        )
      x
    }) %>%
    do.call(rbind, .) %>%
    as.data.frame()
  
  coords <-
    coords %>%
    dplyr::arrange(index)
  coords
}
