#' Parse sample metadata files
#'
#'
#' @param path
#' @param group
#'
#' @return
#'
#' @export

parse_sample_meta <- function(path, group) {
  all_samp_data <- read_csv(path, show_col_types = FALSE)
  all_samp_data$modified_id <- gsub(all_samp_data$sample, pattern = "-", replacement = "_", fixed = TRUE) #TODO: check that this is needed with new data format
  group_data <- strsplit(all_samp_data$report_group, split = ";")
  samp_data <- all_samp_data[map_lgl(group_data, function(x) group %in% x), ]
  return(samp_data)
}


#' Parse reference metadata files
#'
#'
#' @param
#' @param
#' @param
#'
#' @return
#'
#' @export

#parse_reference_metadata
