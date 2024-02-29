#' Print best match table
#'
#' Prints a table with the best matches for each sample given a pairwise distance matrix.
#'
#' @param pairwise_matrix A triangular `matrix` with all pairwise comparisons between samples and references
#' @param sample_data The sample metadata
#' @param ref_data The reference metadata
#' @param metric_name The name of the metric used in the matrix. This is only used to name columns in the output.
#' @param interactive Whether to print an interactive HTML-based table or one for use with PDF
#' @param ... Passed to `DT::datatable`.
#'
#' @return Returns the table to print
#' @export
#'
#' @examples
print_best_match_table <- function(pairwise_matrix, sample_data, ref_data, metric_name, interactive = knitr::is_html_output(), ...) {
  # Create table to print
  output <- do.call(rbind, lapply(sample_data$sample_id, function(id) { # loop over sample IDs and combine results into a table
    ref_samp_comp <- pairwise_matrix[id, colnames(pairwise_matrix) %in% ref_data$reference_id]
    best_match <- which.max(ref_samp_comp)
    next_best_match <-  which.max(ref_samp_comp[, -best_match])
    data.frame(
      check.names = FALSE,
      sample = sample_data$sample_name[sample_data$sample_id == id],
      best_match =  ref_data$reference_name[ref_data$reference_id == names(best_match)],
      value = unname(ref_samp_comp[best_match]),
      next_match = ref_data$reference_name[ref_data$reference_id == names(next_best_match)],
      next_value = unname(ref_samp_comp[next_best_match])
    )
  }))
  names(output) <- c('Sample', 'Best match', metric_name, '2nd Best match', paste('2nd', metric_name))
  row.names(output) <- NULL

  # Print table
  if (interactive) {
    DT::datatable(output, class = "display nowrap", ...) %>%
      formatStyle(colnames(output), "white-space" = "nowrap")

  } else {
    print_static_table(output)
  }
}
