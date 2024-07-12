
#' Print a table of sample metadata
#'
#' Selects columns in the sample metadata to print and format the result for
#' use as a static table in PDF or an interactive table in HTML.
#'
#' @param input The path to one or more folders that contain
#'   pathogensurveillance output or a table in the format of the
#'   [sample_meta_parsed()] output.
#' @param interactive Whether to use an HTML-based interactive format or not (default: TRUE)
#' @param ... Passed to `DT::datatable`.
#'
#' @export
sample_meta_table <- function(input, interactive = knitr::is_html_output(), ...) {

  # Parse the input if it is a file/folder path
  if (is.data.frame(input)) {
    sample_data <- input
  } else {
    sample_data <- sample_meta_parsed(input)
  }

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
