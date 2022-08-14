library(plyr)
library(dplyr)
library(Seurat)
library(ggplot2)
library(varhandle)
library(reshape2)
library(patchwork)
library(viridis)
library(cowplot)

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
    if ((length(inputGeneList) == 0) || ((length(inputGeneList) == 1) && (inputGeneList == ""))) {
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
      aes(x=!!sym(cxName), y=!!sym(cyName), col=log1p(expr)) +
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
  seuratObj <- inputDataObj$seurat_data

  if (length(inputGeneList) == 0) {
    return(NULL)
  }
  
  chooseCells <- TRUE;
  if(length(inputOpts$selected_cluster) > 0){
    chooseCells <- as.character(unlist(seuratObj[["X__cluster"]])) %in% 
      inputOpts$selected_cluster;
  }
  if(length(inputOpts$selected_ctype) > 0){
    chooseCells <- chooseCells & 
      as.character(unlist(seuratObj[["X__cond"]])) %in% 
      inputOpts$selected_ctype;
  }
  if(!all(chooseCells)){
    seuratObj <- subset(seuratObj, cells = which(chooseCells));
    ## refactor groups to stop empties from showing
    seuratObj[["X__clusterCondREV"]] <- factor(unlist(seuratObj[["X__clusterCondREV"]]));
    seuratObj[["X__clusterREV"]] <- factor(unlist(seuratObj[["X__clusterREV"]]));
  }
  ## Create table linking cells to conditions
  tibble(cellID = colnames(seuratObj),
         cell.identity = as.character(seuratObj@meta.data[[
           if(inputOpts$splitByCondition){"X__clusterCondREV"}
           else {"X__clusterREV"}]])) -> cell.type.tbl
  ## Create group summary table for DotPlot graph
  seuratObj[["RNA"]]@counts[inputGeneList, , drop=FALSE] %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    as_tibble() %>%
    pivot_longer(cols = !starts_with("gene"), names_to = "cellID", values_to = "count") %>%
    pivot_wider(names_from = "gene", values_from = "count") %>%
    left_join(cell.type.tbl, by="cellID") %>%
    pivot_longer(cols = !starts_with(c("cellID", "cell.identity")), 
                 names_to = "gene", values_to = "count") %>%
    group_by(cell.identity, gene) %>%
    dplyr::summarise(pctExpressed = round(100*sum(count != 0) / n(), 1),
                     logMeanExpr = ifelse(pctExpressed == 0, 0, 
                                          log1p(mean(count[count > 0])) / log1p(2)),
                     .groups = "keep") -> dotplot.data
  # Draw dotplot graph
  dotplot.data %>%
    ggplot() +
    aes(x=gene, y=cell.identity, size=pctExpressed, col=logMeanExpr) +
    geom_point() +
    scale_size(limits=c(0,100)) +
    scale_colour_gradient(low="lightgrey", high="#e31837") +
    scale_x_discrete(labels=function(x){sub("-ENSM.*", "", x)}) +
    xlab("Features") +
    ylab("Cluster") +
    theme_cowplot() +
    guides(size=guide_legend(title = "% Expressed"),
           colour=guide_colourbar(title = expression(atop(log[2]~Mean,Expression)),
                                  draw.ulim = TRUE, draw.llim = TRUE)) +
    theme(axis.text.x=element_text(angle = 45, hjust=1)) -> res
  if(inputOpts$colour_scale == "Viridis"){
    suppressWarnings(res <- res + scale_color_viridis_c());
  }
  return(res);
}

GetHeatmapPlot <- function(inputDataList, inputDataIndex, inputGeneList, inputOpts) {

  inputDataObj <- inputDataList[[inputDataIndex]]
  seuratObj <- inputDataObj$seurat_data
  doViridis <- (inputOpts$colour_scale == "Viridis");
  
  if (length(inputGeneList) == 0) {
    return(NULL)
  }
  
  chooseCells <- TRUE;
  if(length(inputOpts$selected_cluster) > 0){
    chooseCells <- as.character(unlist(seuratObj[["X__cluster"]])) %in% 
      inputOpts$selected_cluster;
  }
  if(length(inputOpts$selected_ctype) > 0){
    chooseCells <- chooseCells & 
      as.character(unlist(seuratObj[["X__cond"]])) %in% 
      inputOpts$selected_ctype;
  }
  if(!all(chooseCells)){
    seuratObj <- subset(seuratObj, cells = which(chooseCells));
    ## refactor groups to stop empties from showing
    seuratObj[["X__clusterCondREV"]] <- factor(unlist(seuratObj[["X__clusterCondREV"]]));
    seuratObj[["X__clusterREV"]] <- factor(unlist(seuratObj[["X__clusterREV"]]));
  }
  
  ## Create table linking cells to conditions
  tibble(cellID = colnames(seuratObj),
         cell.identity = as.character(seuratObj@meta.data[[
           if(inputOpts$splitByCondition){"X__clusterCondREV"}
           else {"X__clusterREV"}]])) %>%
    arrange(cell.identity, cellID) -> cell.type.tbl
  
  ## determine scale upper limit based on maximum average gene expression
  maxExpr <- ceiling(log1p(max(rowMeans(seuratObj[["RNA"]]@counts))) / log1p(2));
  ## create colour palette
  colPal <- if(doViridis) {viridis(n=256)} else {
    colorRampPalette(colors=c("lightgrey", "#e31837"), space="Lab")(256)};
  ## get group position centre points
  nCells <- ncol(seuratObj);
  nGenes <- length(inputGeneList);
  groupPos.rle <- rle(cell.type.tbl$cell.identity)
  groupPos <- tibble(cell.identity = factor(groupPos.rle$values),
                     gLength = groupPos.rle$lengths,
                     gStart = head(c(1,cumsum(groupPos.rle$lengths)), -1),
                     gEnd = cumsum(groupPos.rle$lengths),
                     gMid = (gStart + gEnd) / 2);

  seuratObj[["RNA"]]@counts[inputGeneList,,drop=FALSE] %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    as_tibble() %>%
    mutate(gene = factor(gene, levels=rev(inputGeneList))) %>%
    ## convert to long format for easier data manipulation
    pivot_longer(cols = !starts_with("gene"), names_to = "cellID", values_to = "count") %>%
    ## change counts to log2-ish data values
    mutate(logCount = ifelse(count == 0, 0, log1p(count) / log1p(2))) %>%
    ## clip to scale limit
    mutate(logCountBreak = 
             round(255 * ifelse(logCount > maxExpr,
                                maxExpr, logCount) / maxExpr) + 1) %>%
    mutate(logCountCol = colPal[logCountBreak]) %>%
    ## add cell categories and sort cells
    left_join(cell.type.tbl, by="cellID") %>%
    arrange(cell.identity, cellID, gene) -> int.tbl
  
  int.tbl %>%
    group_by(gene, cell.identity) %>%
    dplyr::summarise(logCount = max(logCount), .groups="keep") %>%
    left_join(groupPos, by="cell.identity") -> maxSummary.tbl
  
  int.tbl %>%
    ## convert to matrix
    select(gene, cellID, logCountBreak) %>%
    pivot_wider(names_from=cellID, values_from=logCountBreak) %>%
    column_to_rownames(var="gene") %>%
    as.matrix() -> heatmap.data.int
  
  int.tbl %>%
    ## convert to matrix
    select(gene, cellID, logCountCol) %>%
    pivot_wider(names_from=cellID, values_from=logCountCol) %>%
    column_to_rownames(var="gene") %>%
    as.matrix() -> heatmap.data

  ## flip data so that order matches the displayed order on the plot 
  heatmap.data <- heatmap.data[nGenes:1,,drop=FALSE];

  ## Create base visualisation as dotplot
  maxSummary.tbl %>%
    ggplot() +
    aes(x=gMid, y=gene, col=logCount) +
    ## set up extents
    #geom_rect(aes(xmin=0, xmax=nCells, ymin=0, ymax=nGenes), inherit.aes=FALSE) +
    geom_point() +
    ## add heatmap
    #annotation_raster(heatmap.data, 0.5, nCells+0.5, 0.5, nGenes+0.5) +
    suppressWarnings(scale_x_discrete(name="Identity", limits=groupPos$gMid,
                                      labels=groupPos$cell.identity, position="top",
                                      expand=expansion(mult=0))) +
    theme_classic() +
    theme(axis.text.x=element_text(angle = -45, hjust=1),
          axis.line = element_blank()) +
    scale_colour_gradient(low="lightgrey", high="#e31837") +
    guides(col=guide_colourbar(title = expression(atop(log[2]~Mean,Expression)))) +
    geom_rect(aes(xmin=gStart, xmax=gEnd, ymin=nGenes+0.6, ymax=nGenes+0.7), 
              data = groupPos,
              fill = viridis(nrow(groupPos), option = "H", 
                             begin=0.1)[as.integer(groupPos$cell.identity)],
              inherit.aes=FALSE) +
    geom_segment(data=groupPos,
                 x=groupPos$gEnd + 0.5, xend=groupPos$gEnd + 0.5, y=0.5, yend=nGenes+0.5, 
                 inherit.aes=FALSE, lwd=3, col="grey") -> res
  if(doViridis){
    res <- suppressWarnings(res + scale_colour_viridis_c());
  }

  ## Work out amount to add between groups for a consistent gap
  gapAdd <- (nCells * 0.02) / nrow(groupPos);
  ## Add cell-level heatmap data per-group as a raster annotation
  for(gLine in (1:nrow(groupPos))){
    gStart <- unlist(groupPos[gLine, "gStart"]);
    gEnd <- unlist(groupPos[gLine, "gEnd"]);
    gLength <- unlist(groupPos[gLine, "gLength"]);
    subData <- heatmap.data.int[,gStart:gEnd, drop=FALSE];
    subText <- heatmap.data[,gStart:gEnd, drop=FALSE];
    subText <- subText[,hclust(dist(t(subData)))$order, drop=FALSE];
    res +
      annotation_raster(subText,
                        gStart - 0.5 + gapAdd, gEnd + 0.5 - gapAdd,
                        0.5, nGenes+0.5) -> res
  }
  return(res)
}
