#' Get parsed pipeline status data
#'
#' Return a [tibble::tibble()] (table) with the status messages produced by the pathogensuriveillance
#' pipeline. The contents of all status message files found in the given paths
#' will be combined.
#'
#' @param paths The path to one or more folders that contain
#'   pathogensurveillance output.
#'
#' @return A [tibble::tibble()] with the messages from all input paths
#'
#' @export
status_message_parsed <- function(paths) {
  path_data <- status_message_path_table(paths)
  if (nrow(path_data) > 0) {
    output <- do.call(rbind, lapply(1:nrow(path_data), function(index) {
      table <- readr::read_csv(path_data$path[index], col_types = 'ccccc')
      table$report_group_id <- path_data$report_group_id[index]
      return(table)
    }))
    output <- output[, c("sample_id", "reference_id", "report_group_id", "workflow", "level", "message")]
  } else {
    output <- tibble(
      sample_id = character(0),
      reference_id = character(0),
      report_group_id = character(0),
      workflow = character(0),
      level = character(0),
      message = character(0)
    )
  }
  return(output)
}

#' Get parsed sample metadata
#'
#' Return a [tibble::tibble()] (table) with sample metadata. The contents of all sample
#' metadata files found in the given paths will be combined.
#'
#' @param paths The path to one or more folders that contain
#'   pathogensurveillance output.
#'
#' @return A [tibble::tibble()] with the messages from all input paths
#'
#' @export
sample_meta_parsed <- function(paths) {
  dplyr::bind_rows(lapply(sample_meta_path(paths), readr::read_csv, show_col_types = FALSE))
}

#' Get parsed report groups
#'
#' Return a [tibble::tibble()] (table) with report group ID. The groups found in the given
#' paths will be combined.
#'
#' @param paths The path to one or more folders that contain
#'   pathogensurveillance output.
#'
#' @return A [base::character()] vector with the groups from all input paths
#'
#' @export
report_group_parsed <- function(paths) {
  report_group_path_data(paths)$report_group_id
}

#' Get parsed sendsketch results
#'
#' Return a [tibble::tibble()] (table) with the results from sendsketch. The
#' contents of all sendsketch files found in the given paths will be combined.
#'
#' @param paths The path to one or more folders that contain
#'   pathogensurveillance output.
#' @param only_best Only return the best hit for each combination of report
#'   group and sample. For more control/details on how the top hit is selected,
#'   see [sendsketch_best_hits()].
#'
#' @return A [tibble::tibble()] with the sendsketch output combined
#'
#' @export
sendsketch_parsed <- function(paths, only_best = FALSE) {
  path_data <- sendsketch_path_data(paths)

  # Internal function to parse a single file
  parse_one_file <- function(path, report_group_id, sample_id) {
    data <- readr::read_tsv(path,
                            skip = 2,
                            show_col_types = FALSE,
                            col_types = 'ccccdccccddddddddddddddcccdddddddcccc')
    return(dplyr::bind_cols(
      sample_id = rep(sample_id, nrow(data)),
      report_group_id = rep(report_group_id, nrow(data)),
      data
    ))
  }

  # Use purrr to parse all files and combine them into one data frame
  sketch_data <- purrr::map_dfr(1:nrow(path_data), function(i) {
    parse_one_file(path_data$path[i], path_data$report_group_id[i], path_data$sample_id[i])
  })

  # Convert percentage fields from character to numeric
  sketch_data$WKID <- as.numeric(gsub("%", "", sketch_data$WKID))
  sketch_data$ANI <- as.numeric(gsub("%", "", sketch_data$ANI))
  sketch_data$Complt <- as.numeric(gsub("%", "", sketch_data$Complt))

  # Filter for top hits
  if (only_best) {
    sketch_data <- sendsketch_best_hits(sketch_data)
  }

  return(sketch_data)
}
