#!/usr/bin/Rscript

cat("Loading R libraries...");
suppressMessages({
    library(tidyverse);
    library(Seurat);
    library(rjson);
});
cat(" done\n");

configFileName <-
    if(length(commandArgs(TRUE)) > 0){
        commandArgs(TRUE)[1];
    } else {
        "./config.json";
    };

dataSetsToCheck <- if(length(commandArgs(TRUE)) > 1){
        commandArgs(TRUE)[2];
    } else {
        "";
    };

#Start to read in the config file.
json_file <- rjson::fromJSON(file = configFileName);
json_data <- json_file$data;
datasets <- 1:length(json_data);
dataset_names <- sapply(json_data, function(x) x$name);
dataset_selector <- as.list(c(datasets));
names(dataset_selector) <- c(dataset_names);

if(dataSetsToCheck == ""){
    dataSetsToCheck <- dataset_names;
}

errorStates <- list();

for(dPos in datasets){
    dsObj <- json_data[[dPos]];
    dName <- dsObj$name;
    if(!(dName %in% dataSetsToCheck)){
        next;
    }
    cat("\n",dName,":\n", sep="");
    errorStates[[dName]] <- NULL;
    #print(str(dsObj));
    if(!file.exists(sub("data/", "", sub("^/app/", "", dsObj$file)))){
        cat(sprintf(" - Seurat file existence... FAILED - '%s' not found!\n", dsObj$file));
        errorStates[[dName]] <- c(errorStates[[dName]], "FileCheck");
        next;
    } else {
        cat(" - Seurat file existence... PASSED!\n");
    }
    seuratObj <- NULL;
    cat(" - Seurat object loading... ");
    try({
        seuratObj <- readRDS(sub("data/", "", sub("^/app/", "", dsObj$file)));
    });
    if(class(seuratObj) == "Seurat"){
        cat("PASSED!\n");
    } else {
        cat(sprintf("FAILED - loaded class is '%s'!\n", class(seuratObj)));
        errorStates[[dName]] <- c(errorStates[[dName]], "LoadCheck");
        next;
    }
    if(prod(dim(LayerData(nina.sc, "counts"))) > 0){
        cat(" - Gene counts object existence... PASSED!\n");
    } else {
        cat(sprintf(" - Gene counts... FAILED - no count data found!\n"));
        errorStates[[dName]] <- c(errorStates[[dName]], "GeneCounts");
    }
    if(dsObj$embedding %in% names(seuratObj)){
        cat(" - UMAP embedding... PASSED!\n");
    } else {
        cat(sprintf(" - UMAP embedding... FAILED - embedding '%s' not found in Seurat object!\n", dsObj$embedding));
        errorStates[[dName]] <- c(errorStates[[dName]], "UMAPEmbedding");
    }
    seuratMeta <- seuratObj@meta.data;
    if(dsObj$cluster %in% names(seuratMeta)){
        cat(" - Cluster definition... PASSED!\n");
    } else {
        cat(sprintf(" - Cluster definition... FAILED - cluster classification '%s' not found in Seurat object metadata!\n", dsObj$cluster));
        errorStates[[dName]] <- c(errorStates[[dName]], "ClusterDef");
    }
    if(dsObj$condition %in% names(seuratMeta)){
        cat(" - Condition definition... PASSED!\n");
    } else {
        cat(sprintf(" - Condition definition... FAILED - condition '%s' not found in Seurat object metadata!\n", dsObj$condition));
        errorStates[[dName]] <- c(errorStates[[dName]], "ConditionDef");
    }
    if(("diff_ex_file" %in% names(dsObj))){
        deFile <- sub("data/", "", sub("^/app/", "", dsObj$diff_ex_file));
        if(!file.exists(deFile)){
            cat(sprintf(" - DE file... FAILED - file '%s' not found!\n", deFile));
            errorStates[[dName]] <- c(errorStates[[dName]], "DEFileNotFound");
        } else {
            cat(" - DE file... PASSED!\n");
        }
    }
    if(("diff_ex_cluster_file" %in% names(dsObj))){
        deFile <- sub("data/", "", sub("^/app/", "", dsObj$diff_ex_cluster_file));
        if(!file.exists(deFile)){
            cat(sprintf(" - DE cluster file... FAILED - file '%s' not found!\n", deFile));
            errorStates[[dName]] <- c(errorStates[[dName]], "DEClusterFileNotFound");
        } else {
            cat(" - DE cluster file... PASSED!\n");
        }
    }
}

if((length(errorStates) > 0) && (sum(sapply(errorStates, length)) > 0)){
    cat("\n\nSummarising errors:\n");
} else {
    cat("\n\nAll done, no errors found! [please complain to David Eccles <bioinformatics@gringene.org> if there are actually errors]\n");
}
for(dName in names(errorStates)){
    if(length(errorStates[[dName]]) > 0){
        cat("  ", dName, ":", paste(errorStates[[dName]], collapse=";"), sep="");
    }
}

if((length(errorStates) > 0) && (sum(sapply(errorStates, length)) > 0)){
    cat("\n\nTo rerun with only a single dataset, use the following syntax:\n");
    cat(sprintf("       %s <config_file> '<dataset>'\n", "./config_verifier.r"));
    cat(sprintf("  i.e. %s %s '%s'\n", "./config_verifier.r", configFileName, head(names(errorStates), 1)));
}
