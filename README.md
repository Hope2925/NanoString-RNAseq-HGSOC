## Addressing transcriptomic assay heterogeneity for predictive modeling in cancer

## Data
* Expression and log2FC data can be found in `hgsoc-gene-explorer/data`:
  * `expression_data.rds` (Expression of genes in Pre and Post NACT (NanoString, RNA-seq (RPKM), Microarray))
  * `l2fc_data.rds` (Log2FC of genes in longitudinal data (NanoString, RNA-seq, Microarray))
  * `matched_data.rds` (Expression of genes in matched NanoString and RNA-seq (RPKM))
  * `SC_pseudo.rds` (Cell-type specific CPM of genes)
* Metadata can be found either in our original publication or as metadata.rds 

## Analyses
Additional information is available in the READMEs within each folder
* `RNA-Nano-Matched`: Assessing the harmonization of RNA-seq and NanoString based on FFPE samples sequenced with both assays.
* `RNA-PrePost`: 
    * Assessing initial harmonization of RNA-seq, NanoString, & Microarray. 
    * Identifying biomarkers of HGSOC survival with RNA-seq and NanoString models.
* `RNA-scRNA-Microarray`: 
    * scRNA-seq and bulk RNA-seq comparison
    * Cell-type specific trends (expression and clinical correlation)
    * Non-longitudinal Microarray analyses
    * Gene networks (hdWGCNA & WGCNA)

## Gene-explorer Application
* You can access our portal most easily at [https://hopetownsend.shinyapps.io/hgsoc-gene-explorer/](https://hopetownsend.shinyapps.io/hgsoc-gene-explorer/).
* Code and details for running the application locally are found in `hgsoc-gene-explorer`

## Citation
Please cite X
