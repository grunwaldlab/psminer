
#' Print a table of sample metadata
#'
#' Selects columns in the sample metadata to print and format the result for
#' use as a static table in PDF or an interactive table in HTML.
#'
#' @param sample_data A tibble/data.frame with the sample metadata
#' @param interactive Whether to use an HTML-based interactive format or not (default: TRUE)
#' @export
make_ani_nj_tree <- function(sample_data, interactive = TRUE) {

  # Subset and reformat data for printing
  column_key <- c(
    sample = 'Sample ID',
    fastq_1 = 'Forward Reads',
    fastq_2 = 'Reverse Reads',
    reference_id = 'Reference ID',
    reference =  'Reference'
  )
  formatted_data <- sample_data[, names(column_key)]
  colnames(formatted_data) <- column_key

  # Print table
  if (interactive) {
    DT::datatable(formatted_data)
  } else {
    kableExtra::kbl(formatted_data, booktabs = TRUE)
  }

}
