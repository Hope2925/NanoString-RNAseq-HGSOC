load_all_data <- function(testing=FALSE) {
  if (testing) {
    list(
    exp_data = readRDS("hgsoc-gene-explorer/data/expression_data.rds"),
    meta = readRDS("hgsoc-gene-explorer/data/metadata.rds"),
    l2fc_data = readRDS("hgsoc-gene-explorer/data/l2fc_data.rds"), 
    matched_data = readRDS("hgsoc-gene-explorer/data/matched_data.rds"), 
    sc_data = readRDS("hgsoc-gene-explorer/data/SC_pseudo.rds")
  )

  } else {
    list(
    exp_data = readRDS("data/expression_data.rds"),
    meta = readRDS("data/metadata.rds"),
    l2fc_data = readRDS("data/l2fc_data.rds"), 
    matched_data = readRDS("data/matched_data.rds"), 
    sc_data = readRDS("data/SC_pseudo.rds")
  )
  }
  
}

#pre_median_df, post_median_df, full_median_df