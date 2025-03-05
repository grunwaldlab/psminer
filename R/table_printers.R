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


#' Get table of pipeline status data
#'
#' Return a formatted interactive table with the data on the issues encountered
#' by the pipeline, one row for each issue. The contents of all status message
#' files found in the given paths will be combined.
#'
#' @param paths The path to one or more folders that contain
#'   pathogensurveillance output.
#' @param interactive Whether to produce interactive tables
#'   (TRUE) or static tables (FALSE). Defaults to TRUE if the environment
#'   supports HTML output, otherwise FALSE. Interactive tables offer enhanced
#'   browsing capabilities, while static tables are best for printed pdf
#'   reports.
#' @param ... Passed to `DT::datatable`.
#'
#' @return A table with details for errors, warnings, notes.
#'
#' @export
status_message_table <- function(paths, interactive = knitr::is_html_output(), ...) {
  message_data <- status_message_parsed(paths)

  # Make step status column and sort table by it
  symbol_key <- c(
    NOTE = '\u2705',
    WARNING = '\uD83D\uDFE1',
    ERROR = '\u274C'
  )
  message_data <- message_data[order(match(message_data$level, names(symbol_key))), ]
  message_data$status <- symbol_key[message_data$level]

  # Format the analysis step to look nicer
  message_data$workflow <- tools::toTitleCase(gsub(tolower(message_data$workflow), pattern = '_', replacement = ' '))
  message_data$level <- tools::toTitleCase(gsub(tolower(message_data$level), pattern = '_', replacement = ' '))

  # Rename and reorder columns
  col_name_key <- c(
    level = 'Type',
    workflow = 'Pipeline Step',
    message = 'Message',
    sample_id = 'Sample ID',
    reference_id = 'Reference ID'
  )
  print_data <- message_data[, names(col_name_key)]
  colnames(print_data) <- col_name_key

  # Create table for output
  if (interactive) {
    output <- DT::datatable(print_data, options = list(pageLength = 5, autoWidth = TRUE), escape = FALSE, rownames = message_data$status, ...) %>%
      DT::formatStyle(colnames(print_data), "white-space" = "nowrap")
  } else {
    output <- print_static_table(print_data)
  }

  return(output)
}


#' Get summary table of pipeline status data
#'
#' Return a formatted interactive table with the numbers of issues encountered
#' by the pipeline, one row for each pipeline step. The contents of all status
#' message files found in the given paths will be combined.
#'
#' @param paths The path to one or more folders that contain
#'   pathogensurveillance output.
#' @param interactive Whether to produce interactive tables
#'   (TRUE) or static tables (FALSE). Defaults to TRUE if the environment
#'   supports HTML output, otherwise FALSE. Interactive tables offer enhanced
#'   browsing capabilities, while static tables are best for printed pdf
#'   reports.
#' @param ... Passed to `DT::datatable`.
#'
#' @return A table with counts of errors, warnings, notes.
#'
#' @export
status_message_table_summary <- function(paths, interactive = knitr::is_html_output(), ...) {

  col_name_key <- c(
    workflow = 'Pipeline Step',
    errors = 'Errors',
    warnings = 'Warnings',
    notes = 'Notes'
  )

  message_data <- status_message_parsed_summary(paths)

  if (nrow(message_data) == 0) {
    print_data <- list(character(0))[rep(1, length(col_name_key))]
    names(print_data) <- col_name_key
    print_data <- tibble::as_tibble(print_data)
  } else {
    # Make step status column and sort table by it
    message_data$status <- ifelse(message_data$warnings > 0, 2, 1)
    message_data$status <- ifelse(message_data$errors > 0, 3, message_data$status)
    message_data <- message_data[order(message_data$status), ]
    symbol_key <- c(
      '\u2705',
      '\uD83D\uDFE1',
      '\u274C'
    )
    message_data$status <- symbol_key[message_data$status]

    # Format the analysis step to look nicer
    message_data$workflow <- tools::toTitleCase(gsub(tolower(message_data$workflow), pattern = '_', replacement = ' '))

    # Rename and reorder columns
    print_data <- message_data[, names(col_name_key)]
    colnames(print_data) <- col_name_key
  }

  # Create table for output
  if (nrow(print_data) == 0) {
    status_symbols <- character(0)
  } else {
    status_symbols <- message_data$status
  }

  if (interactive) {
    output <- DT::datatable(print_data, options = list(pageLength = 5, autoWidth = TRUE), escape = FALSE, rownames = status_symbols, ...) %>%
      DT::formatStyle(colnames(print_data), "white-space" = "nowrap")
  } else {
    output <- print_static_table(print_data)
  }

  return(output)
}
