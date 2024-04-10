# This file contains internal functions that are generally useful in multiple contexts



#' Wrapper function to print static tables
#'
#' This wrapper function is meant to allow the method of printing tables used by
#' many functions to be changed in the future. This is primarily used for PDF
#' output.
#'
#' @param data The table to print.
#' @param compressed_cols The named of columns to attempt to make shorter by
#'   moving shared text at the beginning and end to footnotes and replacing with
#'   `...`
#' @param max_nchar The number of characters that the longest entry in a column
#'   can be before the `compressed_cols` option takes effect.
#'
#' @keywords internal
print_static_table <- function(data, compressed_cols = NULL, max_nchar = 20) {
  # Compress column contents
  if (!is.null(compressed_cols)) {
    # Dont compress column that already have small values
    compress_needed <- unlist(lapply(compressed_cols, function(col) {
      ! all(is.na(data[[col]])) && max(nchar(data[[col]])) > max_nchar
    }))
    compressed_cols <- compressed_cols[compress_needed]
    if (length(compressed_cols) > 0 ) {
      # Remove starts
      starts <- lapply(data[compressed_cols], shared_char)
      data[compressed_cols] <- lapply(compressed_cols, function(col_name) {
        column <- data[[col_name]]
        start <- starts[[col_name]]
        if (start != "") {
          column <- sub(column, pattern = paste0("^", start), replacement = "…")
        }
        return(column)
      })
      # Remove ends
      ends <- lapply(data[compressed_cols], shared_char, end = TRUE)
      data[compressed_cols] <- lapply(compressed_cols, function(col_name) {
        column <- data[[col_name]]
        end <- ends[[col_name]]
        if (end != "") {
          column <- sub(column, pattern = paste0(end, "$"), replacement = "…")
        }
        return(column)
      })
      # Create footnotes
      footnotes <- unlist(lapply(compressed_cols, function(col_name) {
        start <- starts[[col_name]]
        end <- ends[[col_name]]
        if (start != "" && end != "") {
          note <- paste0('All values in column "', col_name, '" start with "', start, '" and end with "', end, '".')
        } else if (start != "" ) {
          note <- paste0('All values in column "', col_name, '" start with "', start, '".')
        } else if (end != "" ) {
          note <- paste0('All values in column "', col_name, '" end with "', end, '".')
        } else {
          note <- NA_character_
        }
      }))
      names(footnotes) <- compressed_cols
      footnotes <- footnotes[! is.na(footnotes)]
      # Modify column names ( cant figure out how to get superscripts to render correctly )
      colnames(data)[colnames(data) %in% names(footnotes)] <- paste0(colnames(data)[colnames(data) %in% names(footnotes)], ' (', seq_along(footnotes), ')')
    } else {
      footnotes <- character(0)
    }
  } else {
    footnotes <- character(0)
  }

  # Print table
  is_multi_page_table <- nrow(data) > 45
  table <- kableExtra::kbl(data, booktabs = TRUE, longtable = is_multi_page_table) %>%
    kable_styling(full_width = FALSE, latex_options = c("hold_position", "repeat_header", "scale_down"))
  if (length(footnotes) > 0) {
    table <- footnote(table, number = footnotes)
  }
  return(table)
}


#' Identify conserved start/ends of characters
#'
#' Find any starts or ends of strings that are the same in all elements of a character vector
#'
#' @param col The character vector
#' @param end If `TRUE` then return shared end instead of shared start
#'
#' @keywords internal
shared_char <- function(col, end = FALSE) {

  reverse_char <- function(x) {
    split <- strsplit(x, split = "")
    reversed <- lapply(split, rev)
    return(unlist(lapply(reversed, paste, collapse = "")))
  }

  if (all(is.na(col))) {
    return("")
  }
  shorest_length <- min(unlist(lapply(col, nchar)), na.rm = TRUE)
  result <- ""
  if (end) {
    col <- reverse_char(col)
  }
  indexes <- 1:shorest_length
  for (index in indexes) {
    unique_starts <- unique(substr(col, start = 1, stop = index))
    if (length(unique_starts) == 1) {
      result <- unique_starts
    } else {
      break
    }
  }
  if (end) {
    result <- reverse_char(result)
  }
  return(result)
}



#' Format a number for printing
#'
#' Shortens a number to a specified number of significant digits.
#' Taken from https://stackoverflow.com/questions/3245862/format-numbers-to-significant-figures-nicely-in-r
#'
#' @keywords internal
format_number <- function(nums, sig_fig = 4) {
  formatC(signif(nums, digits = sig_fig), digits = sig_fig, format="fg", flag="#")
}
