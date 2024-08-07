library(plyr)
library(dplyr)
library(Seurat)
library(rjson)
library(shiny)
library(shinyjs)
library(shinydashboard)
library(tidyverse)
library(devtools)
library(DT)
library(varhandle)
library(rlist)
library(shinythemes)
library(viridisLite)
library(svglite)
source("utils.R")

#Start to read in the config file.
json_file <- rjson::fromJSON(file = './data/config.json');
json_data <- json_file$data;
#Read the config data
config <- json_file$config;
data_list <- list();
datasets <- 1:length(json_data);
dataset_groups <- sapply(json_data, function(x) x$group);
dataset_names <- sapply(json_data, function(x) x$name);
dataset_selector <- as.list(c(datasets));
names(dataset_selector) <- paste0(dataset_groups,"/",dataset_names);

logMessage <- function(...){
  outFileName <- sprintf("logs/appLog_%s.txt", Sys.Date());
  cat(format(Sys.time()), " ", sprintf(...), "\n", sep="", file=outFileName, append=TRUE);
  cat(format(Sys.time()), " ", sprintf(...), "\n", sep="");
}

#Use only the first dataset in the config file
dataset_group = dataset_groups[[1]]
dataset_name = dataset_names[[1]]
dataset = datasets[[1]]


calc_pt_size <- function(n) { 25 / n ^ 0.33 }

SetAllIdent <- function(object, ids) {
  Idents(object) <- ids
  return(object)
}

GetClusters <- function(object) {
  clusters <- data.frame(cell.name = names(object@active.ident), cluster = as.character(object@active.ident))
  rownames(clusters) <- NULL
  clusters$cell.name <- as.character(clusters$cell.name)
  return(clusters)
}

reloadConfig <- function(session=NULL){
  json_file <<- rjson::fromJSON(file = './data/config.json');
  json_data <<- json_file$data;
  #Read the config data
  config <<- json_file$config;
  logMessage("loading metadata...");
  datasets <<- 1:length(json_data);
  dataset_groups <<- sapply(json_data, function(x) x$group);
  dataset_names <<- sapply(json_data, function(x) x$name);
  dataset_selector <<- as.list(c(datasets));
  names(dataset_selector) <<- paste0(dataset_groups,"/",dataset_names);
  data_list <<- lapply(json_data, read_metadata);
  names(data_list) <<- paste0(sapply(data_list, function(x){x$group}),"/",sapply(data_list, function(x){x$name}));
  logMessage("all metadata loaded.");
  if(!is.null(session)){
    updateSelectInput(session, "selected_group", choices = dataset_groups, selected = dataset_groups[[1]]);
    updateSelectInput(session, "selected_dataset",
                      choices = dataset_names[dataset_groups == dataset_groups[[1]]],
                      selected = dataset_names[[1]]);
  }
}

read_metadata <- function(x){
  list(
    name = x$name,
    group = x$group,
    config = x,
    meta = if(!is.null(x$meta)){ x$meta } else { NULL },
    loaded = FALSE
  )  
}

read_data <- function(x) {
  # load data and metadata specified by the JSON string.
  # x: individual json string, with [name, file, clusters embedding]
  seurat_data <- readRDS(x$file)
  seurat_data <- SetAllIdent(seurat_data, x$cluster)
  ncells <- length(colnames(seurat_data))
  pt_size <- calc_pt_size(ncells)
  if (!is.null(x$pt_size)) {
    pt_size <- x$pt_size
  }
  font_scale <- 1
  if (!is.null(x$font_scale)) {
    font_scale <- x$font_scale;
  }
  colors <- seurat_data@misc[[sprintf("%s_colors", x$cluster)]]
  if (is.null(colors)) {
    set.seed(2);
    colors <- turbo(n_distinct(seurat_data@active.ident));
  }
  condition <- x$condition;
  geneCounts <- LayerData(seurat_data, layer="counts");
  genes <- sort(rownames(geneCounts));
  
  ## Identify potentially useful condition names
  condNames <- sapply(names(seurat_data@meta.data), 
                      function(x){length(unique(seurat_data@meta.data[[x]]))});
  condNames <- unique(c(condition, sort(names(which(condNames > 1 & condNames < 100)))));

  dimEmbedding <- x$embedding;
  if(is.null(dimEmbedding)){
    dimEmbedding <- "umap";
  }

  #Parser additions
  full_embedding <- as.data.frame(Embeddings(seurat_data, reduction=dimEmbedding))
  assign_clust <- as.data.frame(GetClusters(seurat_data))
  colorVec = mapvalues(assign_clust[, 2], from = unique(assign_clust[, 2]), to = toupper(colors))
  df_plot = cbind(full_embedding[,1:2], assign_clust[, 2], colorVec)
  colnames(df_plot) = c("dim1", "dim2", "cluster", "colorVec")
  y_range = max(full_embedding[, 2]) - min(full_embedding[, 2])
  x_domain = max(full_embedding[, 1]) - min(full_embedding[, 1])
  xScaleRatio_clusterPlot = y_range / x_domain
  yScaleRatio_clusterPlot = x_domain / y_range
  coords_title = group_by(df_plot, cluster) %>% dplyr::summarize(x_center = mean(dim1), y_center = mean(dim2))
  if (!is.null(x$label_coordinates)) {
    coords_title <- dplyr::bind_rows(x$label_coordinates)
    colnames(coords_title) <- c("cluster", "x_center", "y_center")
  }
  pt_size <- as.numeric(x$pt_size)
  
  #Add the full description name on mouse over
  if (is.null(x$cluster_name_mapping)) {
    cluster_names <- seurat_data@active.ident %>% levels()
    names(cluster_names) <- cluster_names
    x$cluster_name_mapping <- as.list(cluster_names)
  }
  desc_df = list.flatten(x$cluster_name_mapping)
  source_abbv = names(desc_df)
  dest_desc = as.character(list.flatten(x$cluster_name_mapping))
  df_plot$cluster_description = as.character(mapvalues(df_plot$cluster, from = source_abbv, to = dest_desc))

  #Differential expression data
  if(!is.null(x$diff_ex_file)){
    DE_condition <- read.csv(file = x$diff_ex_file, header = TRUE, sep = ",")
  } else {
    DE_condition <- NULL;
  }
  if(!is.null(x$diff_ex_cluster_file)){
    DE_cluster <- read.csv(file = x$diff_ex_cluster_file, header = TRUE, sep = ",")
  } else {
    DE_cluster <- NULL;
  }
  meta <- NULL;
  if(!is.null(x$meta)){
    meta <- x$meta;
  }
  
  dataConds <- as.character(unlist(seurat_data[[condition]]))
  
  #print(unique(dataConds));
  cellTypes <- unique(unlist(seurat_data[[x$cluster]]));
  clusterConds <- paste(Idents(seurat_data), dataConds, sep="_")
  seurat_data[["X__cluster"]] = as.character(Idents(seurat_data))
  seurat_data[["X__cond"]] = dataConds
  seurat_data[["X__clusterREV"]] = factor(Idents(seurat_data), 
                                          levels=sort(unique(as.character(Idents(seurat_data))), 
                                                      decreasing = TRUE))
  seurat_data[["X__clusterCondREV"]] = factor(clusterConds,
                                              levels=sort(unique(clusterConds), 
                                                          decreasing=TRUE))
  
  return(
    list(
      name = x$name,
      group = x$group,
      config = x,
      loaded = TRUE,
      seurat_data = seurat_data,
      ncells = ncells,
      pt_size = pt_size,
      font_scale = font_scale,
      embedding = dimEmbedding,
      condition = condition,
      colors = colors,
      geneCounts = geneCounts,
      genes = genes,
      meta = meta,
      condNames = condNames,
      clusters = as.character(Idents(seurat_data)),
      cellTypes = dataConds,
      
      #Parser additions
      plot_df = df_plot,
      x_scale_ratio_clusterPlot = xScaleRatio_clusterPlot,
      y_scale_ratio_clusterPlot = yScaleRatio_clusterPlot,
      title_coords = coords_title,
      DE_cluster_data = DE_cluster,
      DE_condition_data = DE_condition,
      cluster_name_mapping = x$cluster_name_mapping

    ))
}

reloadConfig();

server <- function(input, output, session) {
  ## Save system message outputs (i.e. STDERR) to a log file based on the App start time
  
  values <- reactiveValues(selectedGenes = "", selectedCluster = "",
                           conditionVariable = "", dataConds = "",
                           clusterConds = "")
  updateSelectInput(session, "selected_group", choices = dataset_groups, selected = dataset_groups[[1]])
  updateSelectInput(session, "selected_dataset",
                    choices = dataset_names[dataset_groups == dataset_groups[[1]]],
                    selected = dataset_names[[1]])
  
  winDims_debounced <- debounce(reactive(input$winDims), 300);
  
  genes_debounced <- debounce(reactive(input$selected_gene), 4000);

  observeEvent({ input$selected_group }, {
    groupNameChoices <- dataset_names[dataset_groups == input$selected_group];
    updateSelectInput(session, "selected_dataset",
                      choices = groupNameChoices,
                      selected = groupNameChoices[1])
  })

  #Updates dataset index on selection and updates gene list
  current_dataset_index <- eventReactive({ input$selected_dataset }, {
    mergedName <- paste0(input$selected_group, "/", input$selected_dataset);
    logMessage("Changing dataset to '%s'", mergedName);
    current_index <- dataset_selector[[mergedName]];
    if(data_list[[current_index]]$loaded == FALSE){
      logMessage("Loading dataset '%s' from RDS file...", mergedName);
      data_list[[current_index]] <<- read_data(data_list[[current_index]]$config);
      logMessage("Finished loading dataset '%s' from RDS file", mergedName);
    }
    return(current_index)
  }, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #Return current organoid and update values
  organoid <- eventReactive({ current_dataset_index() }, {
    return(data_list[[current_dataset_index()]])
  })
  
  #Update the UI elements on change
  observeEvent({ organoid() }, {
    updateSelectizeInput(session, 'selected_gene', choices = organoid()$genes, server = TRUE)
    updateSelectizeInput(session, 'selected_cluster', choices = sort(unique(organoid()$clusters)))
    updateSelectInput(session, 'conditionVariable', choices = organoid()$condNames)
    ## TODO: consider triggering the condition variable update
    updateSelectizeInput(session, 'selected_ctype', choices = unique(sort(organoid()$cellTypes)));
  })

  #Logging
  observeEvent({ input$client }, {
    logMessage("New client with ip: %s", input$client$ip)
  },
               ignoreNULL = TRUE)

  #Update expression plot from selectize input
  observeEvent({ genes_debounced() }, {
    values$selectedGene <- input$selected_gene
    logMessage("Gene selection from text input: %s", input$selected_gene)
  },
               ignoreNULL = TRUE, ignoreInit = TRUE)
  
  #Update condition names
  observeEvent(input$conditionVariable, {
    seurat_data <- organoid()$seurat_data;
    if(input$conditionVariable %in% names(seurat_data@meta.data)){
      dataConds <- as.character(seurat_data@meta.data[[input$conditionVariable]]);
      clusterConds <- paste(Idents(seurat_data), dataConds, sep="_");
      values$conditionVariable <- input$conditionVariable;
      values$dataConds <- dataConds;
      values$clusterConds <- clusterConds;
      updateSelectizeInput(session, 'selected_ctype',
                           choices = unique(sort(dataConds)));
    }
  }, ignoreNULL=TRUE, ignoreInit=TRUE)

  #Update expression plot on table row click 
  observeEvent({ input$clusterDE_gene_table_rows_selected }, {
    rowid <- input$clusterDE_gene_table_rows_selected
    gene_selected <- current_table()[rowid, 'gene']
    updateSelectizeInput(session, "selected_gene", choices = organoid()$genes, 
                         selected=gene_selected, server=TRUE)
    values$selectedGene <- gene_selected
  },
              ignoreNULL = TRUE, ignoreInit = TRUE)

  #Get plot window width using the cluster plot as a reference
  plot_window_width = eventReactive({ list(winDims_debounced(), input$pairEmbedding) }, {
    if(input$pairEmbedding){
      return(input$winDims[1]/2 - 25)
    } else {
      return(input$winDims[1] - 25)
    }
  })

  #Get plot window height using the cluster plot as a reference (force height = width)
  plot_window_height = eventReactive({ winDims_debounced() }, {
    return(input$winDims[2] - 145)
  })
  
  #Monitor hover information for biplot
  ##TODO: possibly not needed

  #Monitor cluster plot for changes and update selectedCluster field

  #Set the selectedCluster field to nothing when the reset button is clicked
  observeEvent(eventExpr = { input$reset_table }, handlerExpr = {
    values$selectedCluster <- ""
  })

  #Clear memory and reload metadata from configuration file
  observeEvent(eventExpr = { input$reload_config }, handlerExpr = {
    reloadConfig(session=session);
  })
  
  #Set the selectedCluster field to nothing when the the dataset is changed 
  observeEvent(eventExpr = { current_dataset_index() }, handlerExpr = {
    values$selectedCluster <- ""
  })

  current_gene_list <- eventReactive(c({ values$selectedGene }, { current_dataset_index() }), {
    gene_listy = values$selectedGene
    return(gene_listy)
  })
  
  ##GRAPHIC OUTPUTS (ideally here in the same order as the tabs in ui.R)
  output$cluster_plot <- renderPlot({
    GetClusterPlot(data_list, current_dataset_index(), input, values)
  }, width=plot_window_width, height=plot_window_height)
  output$expression_plot <- renderPlot({
    GetExpressionPlot(data_list, current_dataset_index(), genes_debounced(), input, values)
  }, width=plot_window_width, height=plot_window_height)
  output$bi_plot <- renderPlot({
    GetBiPlot(data_list, current_dataset_index(), genes_debounced(), input, values)
  }, width=plot_window_width, height=plot_window_height)
  output$pairVis <- renderPlot({
    GetPairVis(data_list, current_dataset_index(), genes_debounced(), input, values)
  }, width=plot_window_width, height=plot_window_height)
  output$heatmap_plot <- renderPlot({
    GetHeatmapPlot(data_list, current_dataset_index(), genes_debounced(), input, values)
  }, width=plot_window_width, height=plot_window_height)
  output$dot_plot <- renderPlot({
    GetDotPlot(data_list, current_dataset_index(), genes_debounced(), input, values)
  }, width=plot_window_width, height=plot_window_height)
  output$feature_vs_count_plot <- renderPlot({
    GetFeaturesVsCountsPlot(data_list, current_dataset_index(), input, values)
  }, width=plot_window_width, height=plot_window_height)

  ## Metadata description
  output$metadata_text <- renderUI({
    dataMeta <- (data_list[[current_dataset_index()]])$meta;
    if(!is.null(dataMeta)){
      res <- lapply(names(dataMeta), function(x){
        val <- dataMeta[[x]];
        if(length(val) == 1){
          if(startsWith(val, "http")){
            val <- tags$a(href=val, val);
          }
          meta.def <- tags$dd(val);
        } else {
          meta.def <- tags$dd(tags$ul(lapply(val, function(l){
            if(startsWith(l, "http")){
              l <- tags$a(href=l, target="_blank", l);
            }
            tags$li(l);
          })));
        }
        list(tags$dt(x), meta.def)
      });
    }
  });

  clusterString <- eventReactive({ values$selectedCluster }, {
    baseString = "all clusters"
    if (values$selectedCluster != "") {
      baseString = organoid()$cluster_name_mapping[values$selectedCluster]
    }
    return(sprintf("Genes differentially expressed in %s", baseString))
  })
  
  output$DotPlot_table <- DT::renderDT({
    GetDotPlotData(data_list, current_dataset_index(),
                   genes_debounced(), input, values) -> out.data;
    if(is.null(out.data)){
      NULL;
    } else {
      datatable(out.data);
    }
  });
  
  output$cluster_count_table <- DT::renderDT({
    cluster <- factor(Idents(organoid()$seurat_data))
    cluster <- factor(cluster, levels = sort(levels(cluster)));
    condition <- factor(unlist(organoid()$seurat_data[[input$conditionVariable]]));
    condition <- factor(condition, levels=sort(levels(condition)));
    out <- as_tibble(as.matrix(table(cluster, condition))) %>%
      arrange(cluster, condition) %>%
      pivot_wider(names_from="condition", values_from="n")
    datatable(out)
  })
    
  output$DE_cluster_table <- DT::renderDT({
    datatable(organoid()$DE_cluster_data)
  });

  output$DE_condition_table <- DT::renderDT({
    datatable(organoid()$DE_condition_data)
  });

  #TABLE OUTPUT
  #Format the cluster gene table and add links to Addgene and ENSEMBL

  decimal_columns <- c('avg_logFC', 'p_val', 'p_val_adj', 'avg_diff')
  important_columns <- c('gene', 'cluster_name', 'p_val')

  output$clusterDE_gene_table_title <- renderText({ clusterString() })
  output$clusterDE_gene_table <-
    DT::renderDT({
    datatable(organoid()$diff_eq_table,
                rownames = FALSE,
                extensions = c('Responsive'),
                selection = 'single',
                options =
                  list(
                    columnDefs =
                      list(
                        list(responsivePriority = 1, targets = important_columns),
                        list(
                          render = DT::JS(
                            "function(data, type, row, meta) {",
                            "return type === 'display'?",
                            "'<a href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=' + data + '\">' + data + '</a>' : data;",
                            "}"), targets = c(0)) #,
                      )
                  )
                ) %>%
                formatSignif(decimal_columns[decimal_columns %in% colnames(organoid()$diff_eq_table)], 3)
  },
    server = TRUE
    )
  
  output$save_file <- downloadHandler(
    filename = function() {
      tabName <- gsub(" ", "_", input$tabPanel)
      if(input$tabPanel %in% c("Cluster DE", "Condition DE", "Cluster Counts", "Dot Plot Table")){
        paste0(format.Date(Sys.time(), "%Y-%m-%d_%H%M%S_"), tabName, ".csv")
      } else {
        paste0(format.Date(Sys.time(), "%Y-%m-%d_%H%M%S_"), tabName, ".", input$image_type)
      }
    },
    content = function(file) {
      logMessage("Saving %s to file", input$tabPanel);
      if(input$tabPanel == "Condition DE"){
        organoid()$DE_condition_data %>%
          write_csv(file)
      } else if(input$tabPanel == "Cluster DE"){
        organoid()$DE_cluster_data %>%
          write_csv(file)
      } else if(input$tabPanel == "Cluster Counts"){
        cluster <- Idents(organoid()$seurat_data)
        condition <- unlist(organoid()$seurat_data[[organoid()$condition]])
        as_tibble(as.matrix(table(cluster, condition))) %>%
          pivot_wider(names_from="condition", values_from="n") %>%
          write_csv(file)
      } else if(input$tabPanel == "Dot Plot Table"){
        GetDotPlotData(data_list, current_dataset_index(),
                       genes_debounced(), input, values) %>%
          mutate(across(where(is.numeric), signif, 4)) %>%
          write_csv(file)
      } else { ##FILE OUTPUTS (ideally here in the same order as the tabs in ui.R)
        if(input$tabPanel == "Cluster Plot"){
          GetClusterPlot(data_list, current_dataset_index(), input, values)
          ggsave(file, width = 11, height=plot_window_height() / plot_window_width() * 11, bg="white")
        } else if(input$tabPanel == "Expression Plot"){
          GetExpressionPlot(data_list, current_dataset_index(), current_gene_list(), input, values)
          ggsave(file, width = 11, height=plot_window_height() / plot_window_width() * 11, bg="white")
        } else if(input$tabPanel == "Bi Plot"){
          if(input$pairEmbedding){
            plot_grid(nrow = 1,
                      GetBiPlot(data_list, current_dataset_index(), current_gene_list(), input, values),
                      GetPairVis(data_list, current_dataset_index(), current_gene_list(), input, values))
            ggsave(file, width = 22, height=plot_window_height() / plot_window_width() * 11, bg="white")
          } else {
            GetBiPlot(data_list, current_dataset_index(), current_gene_list(), input, values)
            ggsave(file, width = 11, height=plot_window_height() / plot_window_width() * 11, bg="white")
          }
        } else if(input$tabPanel == "Heat Map"){
          GetHeatmapPlot(data_list, current_dataset_index(), current_gene_list(), input, values)
          ggsave(file, width = 11, height=plot_window_height() / plot_window_width() * 11, bg="white")
        } else if(input$tabPanel == "Dot Plot"){
          GetDotPlot(data_list, current_dataset_index(), current_gene_list(), input, values)
          ggsave(file, width = 11, height=plot_window_height() / plot_window_width() * 11, bg="white")
        } else if(input$tabPanel == "Feature/Count Plot"){
          GetFeaturesVsCountsPlot(values$data_list, current_dataset_index(), input, values)
          ggsave(file, width = 11, height=plot_window_height() / plot_window_width() * 11, bg="white")
        } else if(input$tabPanel == "Feature/Count Plot"){
          GetFeaturesVsCountsPlot(values$data_list, current_dataset_index(), input, values)
          ggsave(file, width = 11, height=plot_window_height() / plot_window_width() * 11, bg="white")
        } else {
          png(file, width=2200, height=1600, pointsize=20)
          plot.default(NA, xlim=c(0, 1), ylim=c(0,1))
          invisible(dev.off())
        }
      }
    }
  )
}
