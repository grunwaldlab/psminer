#' Get table of pipeline status data
#'
#' Return a formatted interactive table with the data on the issues encountered
#' by the pipeline, one row for each issue. The contents of all status message
#' files found in the given paths will be combined.
#'
#' @param paths The path to one or more folders that contain
#'   pathogensurveillance output.
#'
#' @return A [DT::datatable] with details for errors, warnings, notes.
#'
#' @export
status_message_table <- function(paths) {
  message_data <- status_message_parsed(paths)

  # Make step status column and sort table by it
  symbol_key <- c(
    NOTE = '✅' ,
    WARNING = '🟡',
    ERROR = '❌'
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
  DT::datatable(print_data, options = list(pageLength = 5, autoWidth = TRUE), escape = FALSE, rownames = message_data$status) %>%
    DT::formatStyle(colnames(print_data), "white-space" = "nowrap")
}


#' Get summary table of pipeline status data
#'
#' Return a formatted interactive table with the numbers of issues encountered
#' by the pipeline, one row for each pipeline step. The contents of all status
#' message files found in the given paths will be combined.
#'
#' @param paths The path to one or more folders that contain
#'   pathogensurveillance output.
#'
#' @return A [DT::datatable] with counts of errors, warnings, notes.
#'
#' @export
status_message_table_summary <- function(paths) {

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
      '✅' ,
      '🟡',
      '❌'
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
  DT::datatable(print_data, options = list(pageLength = 5, autoWidth = TRUE), escape = FALSE, rownames = status_symbols) %>%
    DT::formatStyle(colnames(print_data), "white-space" = "nowrap")
}
