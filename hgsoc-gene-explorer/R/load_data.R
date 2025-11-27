load_all_data <- function() {
  list(
    exp_data = readRDS("data/expression_data.rds"),
    meta = readRDS("data/metadata.rds"),
    l2fc_data = readRDS("data/l2fc_data.rds"), 
    matched_data = readRDS("data/matched_data.rds")
  )
}