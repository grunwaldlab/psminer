#' Find sample metadata path data
#'
#' Return the file path data to the CSV containing the sample metadata for a given
#' pathogensurveillance output folder.
#'
#' @param paths The path to one or more folders that contain pathogensurveillance output.
#' @return `tibble` with `report_group_id` and `path` columns
#' @family path tables
#' @export
sample_meta_path_table <- function(paths) {
  out_paths <- sample_meta_path(paths)
  tibble::tibble(
    report_group_id = find_path_report_group(out_paths),
    path = out_paths
  )
}

#' Find reference metadata path data
#'
#' Return the file path data to the CSV containing the reference metadata for a given
#' pathogensurveillance output folder.
#'
#' @param paths The path to one or more folders that contain pathogensurveillance output.
#' @return `tibble` with `report_group_id` and `path` columns
#' @family path tables
#' @export
ref_meta_path_table <- function(paths) {
  out_paths <- ref_meta_path(paths)
  tibble::tibble(
    report_group_id = find_path_report_group(out_paths),
    path = out_paths
  )
}

#' Find the BUSCO tree path data
#'
#' Return a table with file paths to the Newick formatted tree produced by
#' comparing BUSCO genes of samples and references for a given
#' pathogensurveillance output folder.
#'
#' @param paths The path to one or more folders that contain pathogensurveillance
#'   output.
#' @return `tibble` with `report_group_id` and `path` columns
#' @family path tables
#' @export
busco_tree_path_table <- function(paths) {
  out_paths <- busco_tree_path(paths)
  tibble::tibble(
    report_group_id = find_path_report_group(out_paths),
    path = out_paths
  )
}

#' Find the BUSCO analysis reference path data
#'
#' Return a table with the file path to the CSV with a list of references used
#' in the BUSCO analysis for a given pathogensurveillance output folder.
#'
#' @param paths The path to one or more folders that contain pathogensurveillance
#'   output.
#' @return `tibble` with `report_group_id` and `path` columns
#' @family path tables
#' @export
busco_ref_path_table <- function(paths) {
  out_paths <- busco_ref_path(paths)
  tibble::tibble(
    report_group_id = find_path_report_group(out_paths),
    path = out_paths
  )
}

#' Find the core gene analysis reference path data
#'
#' Return a table with the file path to the CSV with the list of references used in the core
#' gene analysis for a given pathogensurveillance output folder.
#'
#' @param paths The path to one or more folders that contain pathogensurveillance output.
#' @return `tibble` with `report_group_id` and `path` columns
#' @family path tables
#' @export
core_ref_path_table <- function(paths) {
  out_paths <- core_ref_path(paths)
  tibble::tibble(
    report_group_id = find_path_report_group(out_paths),
    path = out_paths
  )
}

#' Find the report group file path data
#'
#' Return a table with the file containing the name of the report group for a
#' given pathogensurveillance output folder.
#'
#' @param paths The path to one or more folders that contain
#'   pathogensurveillance output.
#' @return `tibble` with `report_group_id` and `path` columns
#' @family path tables
#' @export
report_group_path_table <- function(paths) {
  out_paths <- report_group_path(paths)
  tibble::tibble(
    report_group_id = find_path_report_group(out_paths),
    path = out_paths
  )
}


#' Determine the report group for files
#'
#' @keywords internal
find_path_report_group <- function(paths) {
  find_one <- function(path) {
    while (path != dirname(path)) {
      if ("group_id.txt" %in% list.files(path)) {
        return(readLines(file.path(path, "group_id.txt")))
      }
      path <- dirname(path)
    }
    return(NA_character_)
  }
  unlist(lapply(paths, find_one))
}
