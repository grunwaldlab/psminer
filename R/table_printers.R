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
sample_meta_table <- function(input, interactive = FALSE, ...) {

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
#' @param summarize_by What variable to summarize results by, if any. Can be one
#'   of `'sample'`, `'message'`, `'workflow'`, or `NULL`. By default, all values
#'   in the message data is shown on its own row.
#' @param interactive Whether to produce interactive tables
#'   (TRUE) or static tables (FALSE). Defaults to TRUE if the environment
#'   supports HTML output, otherwise FALSE. Interactive tables offer enhanced
#'   browsing capabilities, while static tables are best for printed pdf
#'   reports.
#' @param ... Passed to `DT::datatable` for interactive output.
#'
#' @return A table with details for errors, warnings, notes.
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'psminer')
#' status_message_table(path)
#' status_message_table(path, interactive = TRUE)
#' status_message_table(path, summarize_by = 'sample')
#' status_message_table(path, summarize_by = 'workflow')
#' status_message_table(path, summarize_by = 'message')
#'
#' @export
status_message_table <- function(paths, summarize_by = NULL, interactive = FALSE, ...) {
  # Check parameters
  valid_summarize_by <- c('sample', 'message', 'workflow')
  if (! is.null(summarize_by) && ! summarize_by %in% valid_summarize_by) {
    stop('The `summarize_by` parameter must be NULL or one of: ',
         paste0(valid_summarize_by, collapse = ', '))
  }

  # Find and parse message data
  message_data <- status_message_parsed(paths)

  # Make step status column
  symbol_key <- c(
    NOTE = '\u2705',
    WARNING = '\uD83D\uDFE1',
    ERROR = '\u274C'
  )
  symbol_key <- factor(symbol_key, levels = symbol_key, ordered = TRUE)
  message_data$status <- symbol_key[message_data$level]

  # Summarize by sample
  summarize_by_factor <- function(table, by) {
    table <- do.call(rbind, lapply(split(table, table[[by]]), function(part) {
      data.frame(
        x = unique(part[[by]]),
        notes = sum(part$level == 'NOTE'),
        warnings = sum(part$level == 'WARNING'),
        errors = sum(part$level == 'ERROR'),
        status = max(part$status)
      )
    }))
    colnames(table)[1] <- by
    rownames(table) <- NULL
    return(table)
  }
  if (! is.null(summarize_by)) {
    if (summarize_by == 'sample') {
      message_data <- summarize_by_factor(message_data, 'sample_id')
    }
    if (summarize_by == 'workflow') {
      message_data <- summarize_by_factor(message_data, 'workflow')
    }
    if (summarize_by == 'message') {
      message_data <- summarize_by_factor(message_data, 'message')
      message_data$count <- rowSums(message_data[, c('notes', 'warnings', 'errors')])
      message_data <- message_data[, c('message', 'count', 'status')]
    }
  }

  # Format the analysis step to look nicer
  make_text_nice <- function(x) {
    x <- tools::toTitleCase(gsub(tolower(x), pattern = '_+', replacement = ' '))
    first_letter <- substr(x, 1, 1)
    vapply(seq_len(length(x)), FUN.VALUE = character(1), function(i) {
      sub(x[i], pattern = '^.?', replacement = toupper(first_letter[i]))
    })
  }
  columns_to_prettify <- c('workflow', 'level')
  columns_to_prettify <- columns_to_prettify[columns_to_prettify %in% colnames(message_data)]
  message_data[columns_to_prettify] <- lapply(message_data[columns_to_prettify], make_text_nice)

  # Sort by how bad the status is
  message_data <- message_data[order(message_data$status, decreasing = TRUE), ]

  # Reorder columns
  col_name_key <- c(
    status = 'Status',
    level = 'Type',
    sample_id = 'Sample ID',
    reference_id = 'Reference ID',
    workflow = 'Pipeline Step',
    count = 'Count',
    message = 'Message',
    notes = 'Notes',
    warnings = 'Warnings',
    errors = 'Errors'
  )
  col_name_key <- col_name_key[names(col_name_key) %in% colnames(message_data)]
  message_data <- message_data[, names(col_name_key)]

  # Rename columns
  print_data <- message_data
  colnames(print_data) <- col_name_key

  # Create pretty table to print
  if (interactive) {
    print_data$Status <- NULL
    output <- DT::datatable(print_data, options = list(pageLength = 5, autoWidth = TRUE),
                            escape = FALSE, rownames = as.character(message_data$status), ...)
    output <- DT::formatStyle(output, colnames(print_data), "white-space" = "nowrap")
  } else {
    rownames(print_data) <- NULL
    colnames(print_data)[colnames(print_data) == 'Status'] <- ''
    output <- print_static_table(print_data)
  }
  print(output)

  # Return an invisible less pretty table for downstream analysis
  return(invisible(tibble::as_tibble(message_data)))
}
