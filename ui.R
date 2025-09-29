library(shiny)
library(shinyjs)
library(DT)
## extra packages needed to get this working
library(shinythemes)
library(varhandle)
library(rlist)
library(logging)
library(cowplot)

json_file <- rjson::fromJSON(file = './data/config.json')
json_config <- json_file$config
window_title <- json_config$ui_title
title_link_text <- json_config$title_link_text
title_link_url <- json_config$title_link_url

ui <- fluidPage(
  tags$head(tags$script('
                        var winDims = [0, 0];
                        var plotElt = document;
                        $(document).on("shiny:connected", function(e) {
                            plotElt = document.getElementById("plotContainer");
                            winDims[0] = plotElt.clientWidth;
                            winDims[1] = window.innerHeight;
                            Shiny.onInputChange("winDims", winDims);
                        });
                        $(window).resize(function(e) {
                            winDims[0] = plotElt.clientWidth;
                            winDims[1] = window.innerHeight;
                            Shiny.onInputChange("winDims", winDims);
                        });
                    ')),
  titlePanel(title = window_title,
             windowTitle = window_title),  # this is to display on the headers of web browser.
  shinyjs::useShinyjs(),
  sidebarLayout(fluid = TRUE,
                sidebarPanel(width = 3,
                             selectInput(inputId = "selected_group", label = "Group", choices = NULL),
                             selectInput(inputId = "selected_dataset", label = "Dataset", choices = NULL),
                             selectizeInput(inputId = "selected_gene", label = "Gene(s)", choices = NULL,
                                            multiple=TRUE,
                                            options = list(placeholder = 'Select a gene')),
                             selectizeInput(inputId = "selected_cluster", label = "Cluster(s)", choices = NULL,
                                            multiple=TRUE,
                                            options = list(placeholder = '[All]')),
                             selectInput("conditionVariable", label="Condition Variable", choices=c("condition"), selectize=TRUE),
                             selectizeInput(inputId = "selected_ctype", label = "Condition(s)", choices = NULL,
                                            multiple=TRUE,
                                            options = list(placeholder = '[All]')),
                             fluidRow(column(6,
                                             radioButtons(inputId="colour_scale", label="Colour Scale", 
                                                          choices = c("MIMR Red", "Viridis"))),
                                      column(6,
                                             radioButtons(inputId="image_type", label="Image Type", 
                                                          choices = c("PNG"="png", "PDF"="pdf", "SVG"="svg"),
                                                          selected="png"))),
                             checkboxInput("splitByCondition", label="Split by condition", value = FALSE),
                             checkboxInput("mergeCluster", label="Merge Clusters", value = FALSE),
                             checkboxInput("mergeGenes", label="Merge Gene Counts", value = FALSE),
                             checkboxInput("pairEmbedding", label="Pair with Visualisation", value = FALSE),
                             fluidRow(column(12, downloadButton(outputId = "save_file", label = "Save output"), align = "center"),
                             h5(strong("Controls")),
                             fluidRow(column(6, actionButton(inputId = "reload_config", label = "Reload config"), align = "center"),
                                      column(6, actionButton(inputId = "reset_table", label = "Reset table"), align = "center")))
                ),
                mainPanel(fluid = TRUE, width = 9,
                          fluidRow(width = 12,
                                   tags$div(id="plotContainer",
                                            tabsetPanel(id = "tabPanel",
                                                        tabPanel("Overview",
                                                                 uiOutput(outputId = "metadata_text", inline = TRUE)),
                                                        tabPanel("Cluster Plot",
                                                                 plotOutput(outputId = "cluster_plot", inline = TRUE)),
                                                        tabPanel("Expression Plot",
                                                                 plotOutput(outputId = "expression_plot", inline = TRUE)),
                                                        tabPanel("Bi Plot", tags$div(
                                                          plotOutput(outputId = "bi_plot", inline = TRUE,
                                                                     #hover=hoverOpts("biPlotHover", delay = 300,
                                                                     #                 delayType="debounce"),
                                                                     brush=brushOpts("biPlotBrush")
                                                          ),
                                                          tags$span(conditionalPanel("input.pairEmbedding", style = "display: inline-block;",
                                                                           plotOutput(outputId = "pairVis", inline = TRUE))))
                                                        ),
                                                        tabPanel("Heat Map",
                                                                 plotOutput(outputId = "heatmap_plot", inline = TRUE)),
                                                        tabPanel("Dot Plot",
                                                                 plotOutput(outputId = "dot_plot", inline = TRUE)),
                                                        tabPanel("Feature/Count Plot",
                                                                 plotOutput(outputId = "feature_vs_count_plot", inline = TRUE)),
                                                        tabPanel("Dot Plot Table",
                                                                 DTOutput("DotPlot_table")),
                                                        tabPanel("Condition DE",
                                                                 DTOutput('DE_condition_table')),
                                                        tabPanel("Cluster DE",
                                                                 DTOutput('DE_cluster_table')),
                                                        tabPanel("Cluster Counts",
                                                                 DTOutput('cluster_count_table'))
                                            )
                                   )
                          )
                )
  ),
  windowTitle = window_title
)
