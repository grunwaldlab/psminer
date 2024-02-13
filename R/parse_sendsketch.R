#' Parse SendSketch Output Files
#'
#' Parses multiple SendSketch output files into a single data frame.
#' Each file is read, and a new column with the sample ID derived from the file name is added.
#'
#' @param paths A character vector of file paths to SendSketch output files like those produced by the nextflow pipeline pathogensurveillance.
#' @return A data frame combining all parsed files with an additional `sample_id` column.
#' @importFrom readr read_tsv
#' @importFrom dplyr bind_cols
#' @importFrom purrr map_dfr
#' @export
#' @examples
#' sendsketch_paths <- list.files(file.path(params$inputs, "inputs", "sendsketch"), full.names = TRUE)
#' sketch_data <- parse_sendsketch(sendsketch_paths)
parse_sendsketch <- function(paths) {

  # Check if 'paths' is a character vector
  if (!is.character(paths)) {
    stop("'paths' must be a character vector of file paths.")
  }

  # Check if 'paths' is not empty
  if (length(paths) == 0) {
    stop("'paths' is empty. Please provide valid file paths.")
  }

  # Check if all files exist
  if (any(!file.exists(paths))) {
    missing_files <- paths[!file.exists(paths)]
    stop("The following files do not exist: ", paste(missing_files, collapse = ", "))
  }

  # Internal function to parse a single file
  parse_one_file <- function(path) {
    data <- readr::read_tsv(path,
                            skip = 2,
                            show_col_types = FALSE,
                            col_types = 'ccccdccccddddddddddddddcccdddddddcccc')
    id <- sub(basename(path), pattern = '\\.txt$', replacement = '')
    return(dplyr::bind_cols(sample_id = rep(id, nrow(data)), data))
  }

  # Use purrr to parse all files and combine them into one data frame
  return(purrr::map_dfr(paths, parse_one_file))
}


