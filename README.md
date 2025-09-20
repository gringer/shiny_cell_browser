
## Shiny Single Cell Browser

Interactive visualization of single cell RNAseq datasets. 

  - Visualize cluster distribution, marker gene expression, and cluster-averaged expression of lists of genes. 
  - Select or click on a gene to show its expression on t-SNE/UMAP plots, select a cluster to show its marker genes.
  - Specify pre-analyzed datasets ([Seurat](https://github.com/satijalab/seurat) format) in the JSON config file as data source. Easily switch between different datasets.

## Setting up and launch the App
  
  - Download the App, `git clone https://github.com/gringer/shiny_cell_browser.git`.
  - This application supports Seurat v3, v4, and v5 objects.
  - Install dependencies as listed [below](#dependencies).
  - Prepare data
    - Store Seurat data objects as `.rds` files. Application loading time can be improved by creating a temporary copy and removing unnecessary data prior to saving the data object, e.g. via [`DietSeurat`](https://github.com/satijalab/seurat/issues/3892#issuecomment-756923135):
    ```
    sc.forApp <- sc
    DefaultAssay(sc.forApp) <- "RNA"
    sc.forApp <- DietSeurat(sc.forApp, dimreducs = "umap", assays="RNA")
    saveRDS(sc.forApp, "fileToSave.rds")
    ```
    - Store marker gene differential expression table in `.csv` files (column names must contain `gene` and `cluster`).
    - Optionally, store cluster colors as a vector in `seurat_data@misc[[sprintf("%s_colors", cluster_name)]]`.
  - Specify the file paths and parameters by creating a `data/config.json` file. Follow the example in [`data/example_config.json`](data/example_config.json). The App will load the files on startup. 
  - Specify the visualization config and data file paths by creating a `data/config.json` file and following the example in [`data/example_config.json`](data/example_config.json). 
    - Multiple datasets can be configured in the same browser.
    - The browser-level config includes the browser title and url link
    - The dataset-level config options are listed below:
      - `group`: the dataset group (for categorising different datasets in the same application).
      - `name`: the dataset name.
      - `file`: the `.rds` file path.
      - `cluster`: the name of the column to use for labelling cell clusters, usually used to describe cell types.
      - `condition`: the name of the column to split cells by (usually the same as cluster, but could be a different variable such as treatment).
      - `embedding`: the type of 2D embedding (e.g. tsne or umap).
      - `diff_ex_cluster`: the name of the `@meta.data` cluster id column that corresponds to the cluster ids in the differential expression `csv` file. In most cases, this is the same as `cluster`.
      - `diff_ex_file`: the marker gene differential expression `csv` file.
      - `diff_ex_cluster_file`: alternative differential expression `csv` file, typically used for cluster-specific differential expression.
      - `cluster_name_mapping` (optional): a mapping from the Seurat cluster ids to more readable cluster names.
      - `pt_size` (optional): if set, overrides the automatically computed point size in embedding plots.
      - `font_scale` (optional): if set, scales the font size of cluster labels by this factor.
      - `label_coordinates` (optional): if set, the cluster labels will be placed at these coordinates rather than at the center of each cluster.
  - To launch Single Cell Browser locally, run the following code.  
  ```
  cd shiny_cell_browser
  R -e "shiny::runApp('./', port=4242)
  ## or store the lunch script in run_app.sh and run the following
  ./run_app.sh 
  ```
  - This should launch the web browser at `http://127.0.0.1:4242/`. For other computers in the local network to access the web app, run `R -e "shiny::runApp('./', host='0.0.0.0' port=4242)`. Then visit `your-ip-address:1234`.
  - The App can be deployed on a web server using Docker, [shinyapps.io](https://www.shinyapps.io), or Singularity (see the recipe files provided for build/running information).
  
Example `config.json` file: 

```
{
    "data": [
        {
            "name": "My 1st sample",
            "file": "path/to/seurat/data.rds",
            "cluster": "res.1",
            "embedding": "umap",
            "diff_ex_cluster": "res.1", 
	    "diff_ex_file": "path/to/differential_expression/markers.csv"
        },
        {
            "name": "My 2nd sample",
            "file": "path/to/seurat/data2.rds",
            "cluster": "res.1_rename1",
            "embedding": "tsne",
            "diff_ex_cluster": "res.1", 
            "diff_ex_file": "path/to/differential_expression/markers.csv",

            "cluster_name_mapping": {
                "C1": "Neurons",
                "C2": "Astrocytes",
                "C3": "Neural Progenitors",
                "Note": "cluster_name_mapping is optional"
            },
            "pt_size": 2,
            "font_scale": 0.75
        }
    ],
    "config": {
        "ui_title": "Single Cell Browser",
        "title_link_text": "Optional subtitle (e.g. your lab)",
        "title_link_url": "http://optional-link-to-your-lab.com"
    }

}

```
A full example Seurat v4 dataset can be found [here](https://doi.org/10.5281/zenodo.5540103), which is a single-cell dataset associated with [this paper](https://doi.org/10.1038/s41590-021-01067-0).

## Troubleshooting

* Cluster names must be numeric (i.e. the variable linked to "cluster" in the config file)
* The variable linked to "cluster" should be defined in the Seurat object. The `StashIdent` function can be used to store cluster identities as a new named variable within the Seurat object.

### Dependencies

Check the Dockerfile. Here are the R packages that are needed as of 2025-Sep-19:

    install.packages(c("xml2", "openssl", "rjson"))
    install.packages(c("SeuratObject", "Seurat"))
    install.packages(c("shinydashboard"))
    install.packages(c("DT"))
    install.packages(c("varhandle"))
    install.packages(c("shinyjs"))
    install.packages(c("rlist"))
    install.packages(c("shinythemes"))
    install.packages(c("viridis"))
    install.packages(c("logging"))
    install.packages(c("svglite"))
    install.packages(c("systemfonts"))

If these packages don't install properly, it is likely that additional libraries need to be installed on the system that the shiny cell browser is running on. The Dockerfile may help in working out what additional system libraries need to be installed.
  
## Updates

see [updates.md](UPDATES.md)



