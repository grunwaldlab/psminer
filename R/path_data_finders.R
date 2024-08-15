#' Find sample metadata path data
#'
#' Return the file path data to the CSV containing the sample metadata for a given
#' pathogensurveillance output folder.
#'
#' @param paths The path to one or more folders that contain pathogensurveillance output.
#' @return `tibble` with `report_group_id` and `path` columns
#' @family path tables
#' @export
sample_meta_path_data <- function(paths) {
  make_path_data_with_group(paths, sample_meta_path)
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
ref_meta_path_data <- function(paths) {
  make_path_data_with_group(paths, ref_meta_path)
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
busco_tree_path_data <- function(paths) {
  output <- make_path_data_with_group(paths, busco_tree_path)
  output$ref_id_key <- file.path(dirname(output$path), 'read2tree_ref_meta.csv')
  return(output)
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
busco_ref_path_data <- function(paths) {
  make_path_data_with_group(paths, busco_ref_path)
}

#' Find the core gene analysis reference path data
#'
#' Return a table with the file path to the CSV with the list of references used
#' in the core gene analysis for a given pathogensurveillance output folder.
#'
#' @param paths The path to one or more folders that contain
#'   pathogensurveillance output.
#' @return `tibble` with `report_group_id` and `path` columns
#' @family path tables
#' @export
core_ref_path_data <- function(paths) {
  make_path_data_with_group(paths, core_ref_path)
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
report_group_path_data <- function(paths) {
  make_path_data_with_group(paths, report_group_path)
}

#' Find the variant analysis reference path data
#'
#' Return a table with the file path to the CSV with the list of references used
#' in the variant analysis for a given pathogensurveillance output folder.
#'
#' @param paths The path to one or more folders that contain
#'   pathogensurveillance output.
#' @return `tibble` with `report_group_id` and `path` columns
#' @family path tables
#' @export
variant_ref_path_data <- function(paths) {
  make_path_data_with_group(paths, variant_ref_path)
}

#' Find the status message CSV path data
#'
#' Return a table with the file path to the CSV with the status reports,
#' warnings, and errors for a given pathogensurveillance output folder.
#'
#' @param paths The path to one or more folders that contain
#'   pathogensurveillance output.
#' @return `tibble` with `report_group_id` and `path` columns
#' @family path tables
#' @export
status_message_path_data <- function(paths) {
  make_path_data_with_group(paths, status_message_path)
}

#' Find the POCP matrix path data
#'
#' Return a table with the file path to the CSV with the POCP (percent of
#' conserved protein matrix for a given pathogensurveillance output folder.
#'
#' @param paths The path to one or more folders that contain
#'   pathogensurveillance output.
#' @return `tibble` with `report_group_id` and `path` columns
#' @family path tables
#' @export
pocp_matrix_path_data <- function(paths) {
  make_path_data_with_group(paths, pocp_matrix_path)
}

#' Find the estimated ANI matrix path data
#'
#' Return a table with the file path to the CSV with the approximate ANI
#' (average nucleotide identity) matrix estimated by sourmash for a given
#' pathogensurveillance output folder.
#'
#' @param paths The path to one or more folders that contain
#'   pathogensurveillance output.
#' @return `tibble` with `report_group_id` and `path` columns
#' @family path tables
#' @export
estimated_ani_matrix_path_data <- function(paths) {
  make_path_data_with_group(paths, estimated_ani_matrix_path)
}

#' Find the software version path data
#'
#' Return a table with the file path to the YAML file with the versions of
#' software used for a given pathogensurveillance output folder.
#'
#' @param paths The path to one or more folders that contain
#'   pathogensurveillance output.
#' @return `tibble` with `report_group_id` and `path` columns
#' @family path tables
#' @export
software_version_path_data <- function(paths) {
  make_path_data_with_group(paths, software_version_path)
}

#' Find the core gene analysis path data
#'
#' Return a table with file paths to the Newick formatted trees produced by the core
#' gene analysis for a given pathogensurveillance output folder.
#'
#' @param paths The path to one or more folders that contain pathogensurveillance output.
#' @return `tibble` with `report_group_id` and `path` columns
#' @family path tables
#' @export
core_tree_path_data <- function(paths) {
  output <- make_path_data_with_group(paths, core_tree_path)
  output$cluster_id <- unlist(lapply(1:nrow(output), function(index) {
    id <- sub(basename(output$path[index]), pattern = paste0('^', output$report_group_id[index], '_cluster_'), replacement = '')
    sub(id, pattern = '\\.treefile$', replacement = '')
  }))
  return(output)
}

#' Find the considered NCBI reference metadata path data
#'
#' Return a table with the file paths of the CSVs of metadata for references
#' considered by the pipeline for download for a given pathogensurveillance
#' output folder.
#'
#' @param paths The path to one or more folders that contain
#'   pathogensurveillance output.
#' @return `tibble` with `report_group_id` and `path` columns
#' @family path tables
#' @export
considered_ref_meta_path_data <- function(paths) {
  output <- make_path_data_with_group(paths, considered_ref_meta_path)
  output$family <- unlist(lapply(basename(output$path), function(file_name) {
    sub(file_name, pattern = '\\.tsv$', replacement = '')
  }))
  return(output)
}

#' Find the downloaded reference metadata path data
#'
#' Return a table with the file paths of the CSVs of metadata for references
#' selected and downloaded for each sample for a given pathogensurveillance
#' output folder.
#'
#' @param paths The path to one or more folders that contain
#'   pathogensurveillance output.
#' @return `tibble` with `report_group_id` and `path` columns
#' @family path tables
#' @export
selected_ref_meta_path_data <- function(paths) {
  output <- make_path_data_with_group(paths, selected_ref_meta_path)
  output$sample_id <- unlist(lapply(basename(output$path), function(file_name) {
    sub(file_name, pattern = '\\.tsv$', replacement = '')
  }))
  return(output)
}

#' Find the sendsketch result path data
#'
#' Return a table with the file paths of the sendsketch results for each sample
#' for a given pathogensurveillance output folder.
#'
#' @param paths The path to one or more folders that contain
#'   pathogensurveillance output.
#' @return `tibble` with `report_group_id` and `path` columns
#' @family path tables
#' @export
sendsketch_path_data <- function(paths) {
  output <- make_path_data_with_group(paths, sendsketch_path)
  output$sample_id <- unlist(lapply(basename(output$path), function(file_name) {
    sub(file_name, pattern = '\\.txt$', replacement = '')
  }))
  return(output)
}

#' Find the SNP alignment path data
#'
#' Return a table with the file paths of the SNP alignments for each reference
#' used in the variant analysis for a given pathogensurveillance output folder.
#'
#' @param paths The path to one or more folders that contain
#'   pathogensurveillance output.
#' @return `tibble` with `report_group_id` and `path` columns
#' @family path tables
#' @export
variant_align_path_data <- function(paths) {
  output <- make_path_data_with_group(paths, variant_align_path)
  file_names <- basename(output$path)
  output$ref_id <- unlist(lapply(1:nrow(output), function(index) {
    id <- sub(file_names[index], pattern = paste0('^', output$report_group_id[index], '_'), replacement = '')
    sub(id, pattern = '\\.fasta$', replacement = '')
  }))
  return(output)
}

#' Find the SNP tree path data
#'
#' Return a table with the file paths of the SNP tree for each reference used in
#' the variant analysis for a given pathogensurveillance output folder.
#'
#' @param paths The path to one or more folders that contain
#'   pathogensurveillance output.
#' @return `tibble` with `report_group_id` and `path` columns
#' @family path tables
#' @export
variant_tree_path_data <- function(paths) {
  output <- make_path_data_with_group(paths, variant_tree_path)
  file_names <- basename(output$path)
  output$ref_id <- unlist(lapply(1:nrow(output), function(index) {
    id <- sub(file_names[index], pattern = paste0('^', output$report_group_id[index], '_'), replacement = '')
    sub(id, pattern = '\\.treefile$', replacement = '')
  }))
  output$cluster_id <- unlist(lapply(1:nrow(output), function(index) {
    id <- sub(file_names[index], pattern = paste0('^', output$report_group_id[index], '_', output$ref_id[index], '_'), replacement = '')
    sub(id, pattern = '\\.treefile$', replacement = '')
  }))
  return(output)
}

#' Make a table with paths and report group
#'
#' @keywords internal
make_path_data_with_group <- function(paths, path_func) {
  out_paths <- path_func(paths)

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
  if (length(paths) == 0) {
    return(character(0))
  } else {
    return(unlist(lapply(paths, find_one)))
  }
}
