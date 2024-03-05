#' Print best ANI match table
#'
#' Prints a table with the best matches for each sample given a pairwise distance matrix.
#'
#' @param pairwise_matrix A triangular `matrix` with all pairwise comparisons between samples and references
#' @param sample_data The sample metadata
#' @param ref_data The reference metadata
#' @param interactive Whether to print an interactive HTML-based table or one for use with PDF
#' @param ... Passed to `DT::datatable`.
#'
#' @return Returns the table to print
#' @export
#'
#' @examples
print_ani_table <- function(pairwise_matrix, sample_data, ref_data, interactive = knitr::is_html_output(), ...) {
  # Create table to print
  output <- do.call(rbind, lapply(sample_data$sample_id, function(id) { # loop over sample IDs and combine results into a table
    ref_samp_comp <- pairwise_matrix[id, colnames(pairwise_matrix) %in% ref_data$reference_id]
    best_match <- which.max(ref_samp_comp)
    next_best_match <-  which.max(ref_samp_comp[, -best_match])
    data.frame(
      check.names = FALSE,
      'Sample' = sample_data$sample_name[sample_data$sample_id == id],
      'Best match' =  ref_data$reference_name[ref_data$reference_id == names(best_match)],
      'ANI (%)' = format_number(unname(unlist(ref_samp_comp[best_match]))),
      '2nd Best match' = ref_data$reference_name[ref_data$reference_id == names(next_best_match)],
      '2nd ANI (%)' = format_number(unname(unlist(ref_samp_comp[next_best_match])))
    )
  }))
  row.names(output) <- NULL

  # Print table
  if (interactive) {
    DT::datatable(output, class = "display nowrap", ...) %>%
      formatStyle(colnames(output), "white-space" = "nowrap")

  } else {
    print_static_table(output)
  }
}



#' Print best POCP table
#'
#' Prints a table with the highest  for each sample given a pairwise distance matrix.
#'
#' @param pairwise_matrix A `matrix` with all pairwise comparisons between samples and references
#' @param sample_data The sample metadata
#' @param ref_data The reference metadata
#' @param interactive Whether to print an interactive HTML-based table or one for use with PDF
#' @param ... Passed to `DT::datatable`.
#'
#' @return Returns the table to print
#' @export
#'
#' @examples
print_pocp_table <- function(pairwise_matrix, sample_data, ref_data, interactive = knitr::is_html_output(), ...) {
  # Create table to print
  output <- do.call(rbind, lapply(sample_data$sample_id, function(id) { # loop over sample IDs and combine results into a table
    ref_samp_comp_1 <- pairwise_matrix[id, colnames(pairwise_matrix) %in% ref_data$reference_id, drop = FALSE]
    ref_samp_comp_2 <- t(pairwise_matrix[rownames(pairwise_matrix) %in% ref_data$reference_id, id, drop = FALSE])
    ref_samp_comp_max <- mapply(max, unlist(ref_samp_comp_1), unlist(ref_samp_comp_2))
    best_match <- which.max(ref_samp_comp_max)
    data.frame(
      check.names = FALSE,
      'Sample' = sample_data$sample_name[sample_data$sample_id == id],
      'Best match' =  ref_data$reference_name[ref_data$reference_id == names(best_match)],
      'POCP 1 (%)' = format_number(unname(unlist(ref_samp_comp_1[best_match]))),
      'POCP 2 (%)' = format_number(unname(unlist(ref_samp_comp_2[best_match])))
    )
  }))
  row.names(output) <- NULL

  # Print table
  if (interactive) {
    DT::datatable(output, class = "display nowrap", ...) %>%
      formatStyle(colnames(output), "white-space" = "nowrap")

  } else {
    print_static_table(output)
  }
}


