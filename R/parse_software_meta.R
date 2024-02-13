#' Parse Software Version Metadata from YAML File
#'
#' Reads a YAML file containing software version information and transforms it into a tibble.
#' Each software module and program, along with its version, is extracted and organized.
#'
#' @param version_path A character string specifying the path to the YAML file with version metadata.
#' @return A tibble with columns `module`, `program`, and `version` detailing the software versions.
#' @importFrom readr read_yaml
#' @importFrom tibble tibble
#' @importFrom purrr map_chr
#' @importFrom stringr str_split
#' @export
#' @examples
#' version_path <- file.path(params$inputs, "inputs", "versions.yml")
#' version_data <- parse_software_meta(version_path)
parse_software_meta <- function(version_path) {

  # Input validation
  if (!is.character(version_path) || length(version_path) != 1) {
    stop("version_path must be a single character string.")
  }

  # Check if file exists
  if (!file.exists(version_path)) {
    stop("The specified version_path does not exist: ", version_path)
  }

  # Reading the YAML file
  raw_version_data <- tryCatch(
    {
      unlist(readr::read_yaml(version_path))
    },
    error = function(e) {
      stop("Failed to read YAML file: ", e$message)
    }
  )

  # Transforming the data into a tibble
  version_data <- tibble::tibble(
    module = purrr::map_chr(stringr::str_split(names(raw_version_data), pattern = '\\.', n = 2), `[`, 1),
    program = purrr::map_chr(stringr::str_split(names(raw_version_data), pattern = '\\.', n = 2), `[`, 2),
    version = unname(raw_version_data),
  )

  return(version_data)
}
