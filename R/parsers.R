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
  path_data <- status_message_path_data(paths)
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
#' @return A [tibble::tibble()] with the sample metadata
#'
#' @export
sample_meta_parsed <- function(paths) {
  dplyr::bind_rows(lapply(sample_meta_path(paths), readr::read_csv, show_col_types = FALSE))
}

#' Get parsed reference metadata
#'
#' Return a [tibble::tibble()] (table) with reference metadata. The contents of
#' all reference metadata files found in the given paths will be combined. This
#' only contains data about references used by at least one step of the pipeline.
#'
#' @param paths The path to one or more folders that contain
#'   pathogensurveillance output.
#'
#' @return A [tibble::tibble()] with reference metadata
#'
#' @export
ref_meta_parsed <- function(paths) {
  dplyr::bind_rows(lapply(ref_meta_path(paths), readr::read_csv, show_col_types = FALSE))
}


#' Get estimated ANI distance matrix
#'
#' Return a list of [base::data.frame()]s with the estimated pairwise ANI values
#' calculated by sourmash.
#'
#' @param paths The path to one or more folders that contain
#'   pathogensurveillance output.
#' @return a `list` of ANI matrices
#'
#' @export
estimated_ani_matrix_parsed <- function(paths) {
  ani_matrix_paths <- estimated_ani_matrix_path(paths)
  output <- lapply(ani_matrix_paths, function(path) {
    ani_matrix <- read.csv(path, check.names = FALSE)
    rownames(ani_matrix) <- colnames(ani_matrix)
    return(ani_matrix)
  })
  names(output) <- ani_matrix_paths
  return(output)
}

#' Get POCP matrix
#'
#' Return a list of [base::data.frame()]s with the POCP values based on a core
#' gene analysis using pirate.
#'
#' @param paths The path to one or more folders that contain
#'   pathogensurveillance output.
#' @return a `list` of POCP matrices
#'
#' @export
pocp_matrix_parsed <- function(paths) {
  matrix_paths <- pocp_matrix_path(paths)
  output <- lapply(matrix_paths, function(path) {
    pocp_matrix <- read.csv(path, check.names = FALSE, sep = '\t')
    rownames(pocp_matrix) <- colnames(pocp_matrix)
    return(pocp_matrix)
  })
  names(output) <- matrix_paths
  return(output)
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

#' Parse Software Version Metadata from YAML File
#'
#' Reads a YAML file containing software version information and transforms it into a tibble.
#' Each software module and program, along with its version, is extracted and organized.
#'
#' @param paths The path to one or more folders that contain
#'   pathogensurveillance output.
#' @return A [tibble::tibble()] with columns `module`, `program`, and `version` detailing the software versions.
#' @importFrom yaml read_yaml
#' @importFrom tibble tibble
#' @importFrom purrr map_chr
#' @importFrom stringr str_split
#' @export
#' @examples
#' version_path <- file.path(params$inputs, "inputs", "versions.yml")
#' version_data <- parse_software_meta(version_path)
software_version_parsed <- function(paths) {

  parse_one <- function(version_path) {
    # Reading the YAML file
    raw_version_data <- tryCatch({
      unlist(yaml::read_yaml(version_path))
    },
    error = function(e) {
      stop("Failed to read YAML file: ", e$message)
    })
    # Transforming the data into a tibble
    version_data <- tibble::tibble(
      module = purrr::map_chr(stringr::str_split(
        names(raw_version_data), pattern = '\\.', n = 2
      ), `[`, 1),
      program = purrr::map_chr(stringr::str_split(
        names(raw_version_data), pattern = '\\.', n = 2
      ), `[`, 2),
      version = unname(raw_version_data),
    )
    return(version_data)
  }

  unique(dplyr::bind_rows(lapply(software_version_path(paths), parse_one)))
}
