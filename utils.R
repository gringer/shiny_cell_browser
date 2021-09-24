library(Seurat)
library(dplyr)
library(ggplot2)
library(plyr)
library(dplyr)
library(varhandle)
library(reshape2)
library(patchwork)
library(viridis)

##Helper calculation and data functions

#All plots were desigend around a width of 330 pixels, so we scale around that for different screen sizes
scaleRatio <- function(inputWidth) {
  return(inputWidth / 330)
}

PercentAbove <- function(x) {
  return(length(x = x[x > 0]) / length(x = x))
}

MaxMutate <- function(x) {
  return(x / max(x))
}

get_shared_genes <- function(inputGeneList1, inputGeneList2, topN) {
  gene_list1 = dplyr::distinct(as.data.frame(inputGeneList1) %>% mutate_if(is.factor, as.character))
  gene_list2 = dplyr::distinct(as.data.frame(inputGeneList2) %>% mutate_if(is.factor, as.character))
  colnames(gene_list1) = c("gene")
  colnames(gene_list2) = c("gene")
  gene_list1$geneMod <- toupper(sub("-ENS.*$", "", gene_list1$gene))
  gene_list2$geneMod <- toupper(sub("-ENS.*$", "", gene_list2$gene))
  sharedPoss <- na.omit(match(gene_list1$geneMod, gene_list2$geneMod))
  if(length(sharedPoss) > 0){
    geneList <- head(gene_list2$gene[sharedPoss], topN)
    cat("Found genes: ", geneList)
    return(geneList)
  } else {
    cat("No genes found\n")
    print(head(gene_list1))
    print(head(gene_list2))
    return(NULL)
  }
}

FetchGenes <- function(
  object,
  vars.all = NULL
) {
  cells.use <- colnames(object)
  gene.check <- vars.all %in% rownames(GetAssayData(object))
  if (!all(gene.check)) {
    for (my.var in vars.all) {
      if (!(my.var %in% rownames(GetAssayData(object)))) {
        stop(paste("Error:", my.var, "not found"))
      }
    }
  }
  data.expression <- GetAssayData(object)
  data.expression <- t(as.matrix(data.expression[vars.all[gene.check], cells.use, drop=FALSE]))
  return(as.matrix(x = data.expression))
}

##Helper plotting functions

GetClusterPlot <- function(inputDataList, inputDataIndex, inputOpts) {

  inputDataObj = inputDataList[[inputDataIndex]]
  DimPlot(inputDataObj$seurat_data, cols=inputDataObj$colors, reduction=inputDataObj$embedding,
          pt.size=2, split.by=if(inputOpts$splitByCondition){inputDataObj$condition} else {NULL})
}

GetPlotData <- function(inputDataObj, inputGene) {
  print(inputGene)
  if(length(inputGene) == 1){
    single_gene <- mutate(inputDataObj$plot_df[, 1:2], 
                          gene = as.numeric(FetchGenes(inputDataObj$seurat_data, inputGene))) %>% arrange(gene)
    colnames(single_gene) = c("dim1", "dim2", "gene")
    return(single_gene)
  } else {
    return(NULL)
  }
}

GetExpressionPlot <- function(inputDataList, inputDataIndex, inputGeneList, inputOpts) {

    inputDataObj <- inputDataList[[inputDataIndex]]
    seuratObj <- inputDataObj$seurat_data

    #On initialization, check if the inputGene is not defined
    if ((length(inputGeneList) == 0) || (inputGeneList == "")) {
      return(NULL)
    }
    ## Fetch expression data
    GetAssayData(seuratObj, slot="counts", assay="RNA")[inputGeneList,, drop=FALSE] %>%
        data.frame() %>%
        rownames_to_column("feature") %>%
        as_tibble() %>%
        pivot_longer(cols = -1, names_to="cell", values_to="expr") -> feature.tbl;
    colnames(feature.tbl) <- sub("-ENS.*$", "", colnames(feature.tbl));
    ## Fetch dimensional reduction
    Embeddings(seuratObj, reduction=inputDataObj$embedding) %>%
        data.frame() %>%
        rownames_to_column("cell") %>%
        as_tibble() %>%
        mutate(cell=gsub("-", ".", cell)) -> cell.tbl;
    ## identify reduction names
    cxName <- colnames(cell.tbl)[2];
    cyName <- colnames(cell.tbl)[3];
    rangeX <- range(cell.tbl[,2]);
    rangeY <- range(cell.tbl[,3]);
    ## create cell groups
    if((length(inputOpts$selected_ctype) > 1) | (length(inputOpts$selected_cluster) > 1)){
      cell.tbl$group <- as.character(unlist(
        seuratObj[[if(length(inputOpts$selected_ctype) > 1){
          if(length(inputOpts$selected_cluster) > 1){
            "X__clusterCondREV"
          } else {"X__cond"}
        } else if(length(inputOpts$selected_cluster) > 1){
          "X__cluster"
        }]]
      ));
    } else {
      cell.tbl$group = 1;
    }
    ## filter cluster / cell type at cell level
    cell.tbl$includeFilter <- TRUE;
    if(length(inputOpts$selected_cluster) > 0){
      cell.tbl$includeFilter <- as.character(unlist(seuratObj[["X__cluster"]])) %in% 
                 inputOpts$selected_cluster;
    }
    if(length(inputOpts$selected_ctype) > 0){
      cell.tbl$includeFilter <- cell.tbl$includeFilter & 
        as.character(unlist(seuratObj[["X__cond"]])) %in% 
                 inputOpts$selected_ctype;
    }
    cell.tbl <- cell.tbl %>% filter(includeFilter);
    ## merge expression + cell data
    cell.tbl %>%
        left_join(feature.tbl, by="cell") %>%
        arrange(expr) -> merged.tbl
    ## plot object
    merged.tbl %>% ggplot() +
      aes(x=!!sym(cxName), y=!!sym(cyName), col=expr) +
      xlim(rangeX[1], rangeX[2]) +
      ylim(rangeY[1], rangeY[2]) +
      geom_point() +
      theme_bw() +
      theme(strip.background = element_blank(), strip.text.x = element_text(face="bold")) -> res;
    if((length(inputOpts$selected_cluster) > 1) | (length(inputOpts$selected_ctype) > 1)){
      res <- res + facet_wrap(~ group + feature);
    } else if(length(inputGeneList) > 1){
      res <- res + facet_wrap(~ feature);
    }
    ## update dot colours
    if(inputOpts$colour_scale == "Viridis"){
        res <- res + scale_colour_viridis();
    } else {
        res <- res + scale_colour_gradient(low = "lightgrey", high="#e31837");
    }
    return(res);
}

GetDotPlot <- function(inputDataList, inputDataIndex, inputGeneList, inputOpts) {

  inputDataObj <- inputDataList[[inputDataIndex]]

  if (length(inputGeneList) == 0) {
    return(NULL)
  }
  else {
    res <- DotPlot(inputDataObj$seurat_data, features=inputGeneList,
            dot.min=0.0001, scale.by="size", scale=TRUE,
            col.min=0, col.max=1,
            group.by=if(inputOpts$splitByCondition){"X__clusterCondREV"} else {"X__clusterREV"},
            cols = c("lightgrey", "#e31837")) +
      ylab("Cluster") +
      scale_x_discrete(labels=function(x){sub("-ENSM.*", "", x)}) +
      theme(axis.text.x=element_text(angle=45, hjust=1))
    if(inputOpts$colour_scale == "Viridis"){
      suppressWarnings(res <- res + scale_color_viridis_c())
    }
    return(res)
  }
}

GetHeatmapPlot <- function(inputDataList, inputDataIndex, inputGeneList, inputOpts) {

  inputDataObj = inputDataList[[inputDataIndex]]

  if (length(inputGeneList) == 0) {
    return(NULL)
  }
  else {
    print(table(unlist(inputDataObj$seurat_data[["X__clusterCondREV"]])))
    res <- DoHeatmap(inputDataObj$seurat_data, features=inputGeneList,
            slot="counts", assay="RNA", raster=FALSE,
            group.by=if(inputOpts$splitByCondition){"X__clusterCondREV"} else {"X__clusterREV"}) +
      scale_y_discrete(labels=function(x){sub("-ENSM.*", "", x)})
    if(inputOpts$colour_scale == "Viridis"){
      suppressWarnings(res <- res + scale_fill_viridis())
    }
    return(res)
  }
}
