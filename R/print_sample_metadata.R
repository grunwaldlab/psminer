
#' Print a table of sample metadata
#'
#' Selects columns in the sample metadata to print and format the result for
#' use as a static table in PDF or an interactive table in HTML.
#'
#' @param sample_data A tibble/data.frame with the sample metadata
#' @param interactive Whether to use an HTML-based interactive format or not (default: TRUE)
#' @param ... Passed to `DT::datatable`.
#'
#' @export
print_sample_metadata <- function(sample_data, interactive = knitr::is_html_output(), ...) {

  # Subset and reformat data for printing
  # column_key <- c(
  #   sample_id = 'Sample ID',
  #   path = 'Read paths',
  #   path_2 = 'Reverse Reads',
  #   reference_id = 'Reference ID',
  #   reference_name =  'Reference'
  # )
  # formatted_data <- sample_data[, names(column_key)]
  # colnames(formatted_data) <- column_key
  formatted_data <- sample_data

  # Print table
  if (interactive) {
    DT::datatable(formatted_data, class = "display nowrap", ...) %>%
      formatStyle(colnames(formatted_data), "white-space" = "nowrap")

  } else {
    print_static_table(formatted_data, compressed_cols = c('Forward Reads', 'Reverse Reads', 'Reference'))
  }

}
