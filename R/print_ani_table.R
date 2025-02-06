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
    ref_samp_comp <- pairwise_matrix[id, colnames(pairwise_matrix) %in% ref_data$ref_id, drop = FALSE]
    best_match <- which.max(ref_samp_comp)
    if (ncol(ref_samp_comp) >= 2) {
      next_best_match <- which.max(ref_samp_comp[, -best_match])
      next_best_match_name <- ref_data$ref_name[ref_data$ref_id == names(next_best_match)]
      next_best_match_ani <- format_number(unname(unlist(ref_samp_comp[next_best_match])))
    } else {
      next_best_match_name <- "NA (Too few references)"
      next_best_match_ani <- "NA (Too few references)"
    }
    data.frame(
      check.names = FALSE,
      'Sample' = sample_data$name[sample_data$sample_id == id],
      'Best match' =  ref_data$ref_name[ref_data$ref_id == names(best_match)],
      'ANI (%)' = format_number(unname(unlist(ref_samp_comp[best_match]))),
      '2nd Best match' = next_best_match_name,
      '2nd ANI (%)' = next_best_match_ani
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
  sample_ids_in_matrix <- sample_data$sample_id[sample_data$sample_id %in% row.names(pairwise_matrix)]
  output <- do.call(rbind, lapply(sample_ids_in_matrix, function(id) { # loop over sample IDs and combine results into a table
    if (any(colnames(pairwise_matrix) %in% ref_data$ref_id)) {
      ref_samp_comp <- pairwise_matrix[id, colnames(pairwise_matrix) %in% ref_data$ref_id, drop = FALSE]
      best_match <- which.max(ref_samp_comp)
      best_ref_match_name <- ref_data$ref_name[ref_data$ref_id == names(best_match)]
      best_ref_match_pocp <- format_number(unname(unlist(ref_samp_comp[best_match])))
    } else {
      best_ref_match_name <- "NA (No reference used)"
      best_ref_match_pocp <- "NA (No reference used)"
    }
    if (nrow(sample_data) > 1) {
      samp_samp_comp <- pairwise_matrix[id, ! colnames(pairwise_matrix) %in% ref_data$ref_id & colnames(pairwise_matrix) != id, drop = FALSE]
      best_match <- which.max(samp_samp_comp)
      best_samp_match_name <- sample_data$name[sample_data$sample_id == names(best_match)]
      best_samp_match_pocp <- format_number(unname(unlist(samp_samp_comp[best_match])))
    } else {
      best_samp_match_name <- "NA (Too few samples)"
      best_samp_match_pocp <- "NA (Too few samples)"
    }
    data.frame(
      check.names = FALSE,
      'Sample' = sample_data$name[sample_data$sample_id == id],
      'Best reference match' =  best_ref_match_name,
      'Reference POCP (%)' = best_ref_match_pocp,
      'Best sample match' =  best_samp_match_name,
      'Sample POCP (%)' = best_samp_match_pocp
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


