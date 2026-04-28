## Addressing transcriptomic assay heterogeneity for predictive modeling in cancer

## Data
### Most easily accessible and readable:
* Expression and log2FC data can be found in `hgsoc-gene-explorer/data`:
  * `expression_data.rds` (Expression of genes in Pre and Post NACT (NanoString, RNA-seq (RPKM), Microarray))
  * `l2fc_data.rds` (Log2FC of genes in longitudinal data (NanoString, RNA-seq, Microarray))
  * `matched_data.rds` (Expression of genes in matched NanoString and RNA-seq (RPKM))
  * `SC_pseudo.rds` (Cell-type specific CPM of genes)
* Metadata can be found either in our original publication or as metadata.rds 
### All:
* All data used for analyses (as detailed in the subsequent analysis folders) are found in the `data/` folder. 
* TOMs from WGCNA and hdWGCNA exceeded the file size limit for Github and are therefore instead stored on [Zenodo](https://zenodo.org/records/19196835).
* BAMS are not provided also due to size limit and because some samples are protected by European standards that prevent sequencing data release. Code used and all counts from BAMS are provided instead.
* The following files were too large for Github and can be instead found on [Zenodo](https://zenodo.org/records/19865826):
    * `data/counts/GSE165897_UMIcounts_HGSOC.tsv.gz`
    * `data/final_counts/GSE165897_HGSOC_SCEobject.rds`
    * `data/Precog/HGSOC_Crijins_GSE13876_GPL7759.rds`

## Results
* Results for analyses can be found within each Analysis folder (detailed below). Some results may be used as data for other sections and therefore may be found in the data folder above.

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
