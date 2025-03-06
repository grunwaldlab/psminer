#' Find sample metadata file path
#'
#' Return the file path to the TSV containing the sample metadata for a given
#' pathogensurveillance output folder.
#'
#' @param paths The path to one or more folders that contain pathogensurveillance output.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'psminer')
#' sample_meta_path(path)
#'
#' @export
sample_meta_path <- function(paths) {
  find_static_file_path(paths, 'sample_data.tsv', 'sample metadata file')
}

#' Find reference metadata file path
#'
#' Return the file path to the TSV containing the reference metadata for a given
#' pathogensurveillance output folder.
#'
#' @param paths The path to one or more folders that contain pathogensurveillance output.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'psminer')
#' ref_meta_path(path)
#'
#' @export
ref_meta_path <- function(paths) {
  find_static_file_path(paths, 'reference_data.tsv', 'reference metadata file')
}

#' Find the BUSCO tree file path
#'
#' Return the file path to the Newick formatted tree produced by comparing BUSCO
#' genes of samples and references for a given pathogensurveillance output
#' folder.
#'
#' @param paths The path to one or more folders that contain pathogensurveillance output.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'psminer')
#' busco_tree_path(path)
#'
#' @export
busco_tree_path <- function(paths) {
  find_static_dir_path(paths, 'busco_trees', 'busco gene tree', file_required = FALSE, dir_required = FALSE)
}

#' Find the BUSCO analysis reference data file path
#'
#' Return the file path to the TSV with a list of references used in the BUSCO
#' analysis for a given pathogensurveillance output folder.
#'
#' @param paths The path to one or more folders that contain pathogensurveillance output.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'psminer')
#' busco_ref_path(path)
#'
#' @export
busco_ref_path <- function(paths) {
  find_static_file_path(paths, 'busco_tree_references.tsv', 'busco anaylsis reference file', file_required = FALSE)
}

#' Find the core gene analysis reference data file path
#'
#' Return the file path to the TSV with the list of references used in the core
#' gene analysis for a given pathogensurveillance output folder.
#'
#' @param paths The path to one or more folders that contain pathogensurveillance output.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'psminer')
#' core_ref_path(path)
#'
#' @export
core_ref_path <- function(paths) {
  find_static_file_path(paths, 'core_gene_tree_references.tsv', 'core gene analysis reference file', file_required = FALSE)
}

#' Find the report group file path
#'
#' Return the file containing the name of the report group for a given
#' pathogensurveillance output folder.
#'
#' @param paths The path to one or more folders that contain pathogensurveillance output.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'psminer')
#' run_info_path(path)
#'
#' @export
run_info_path <- function(paths) {
  find_static_file_path(paths, 'pathogensurveillance_run_info.yml', 'report group name file')
}

#' Find the variant analysis reference data file path
#'
#' Return the file path to the TSV with the list of references used in the
#' variant analysis for a given pathogensurveillance output folder.
#'
#' @param paths The path to one or more folders that contain pathogensurveillance output.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'psminer')
#' variant_ref_path(path)
#'
#' @export
variant_ref_path <- function(paths) {
  find_static_file_path(paths, 'mapping_references.tsv', 'variant analysis reference file')
}

#' Find the status message TSV file path
#'
#' Return the file path to the TSV with the status reports, warnings, and errors
#' for a given pathogensurveillance output folder.
#'
#' @param paths The path to one or more folders that contain pathogensurveillance output.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'psminer')
#' status_message_path(path)
#'
#' @export
status_message_path <- function(paths) {
  find_static_file_path(paths, 'messages.tsv', 'status message data file', file_required = FALSE)
}

#' Find the POCP matrix file path
#'
#' Return the file path to the TSV with the POCP (percent of conserved protein
#' matrix for a given pathogensurveillance output folder.
#'
#' @param paths The path to one or more folders that contain pathogensurveillance output.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'psminer')
#' pocp_matrix_path(path)
#'
#' @export
pocp_matrix_path <- function(paths) {
  find_static_file_path(paths, 'pocp.tsv', 'pocp matrix file', file_required = FALSE)
}

#' Find the estimated ANI matrix file path
#'
#' Return the file path to the CSV with the approximate ANI (average nucleotide identity)
#' matrix estimated by sourmash for a given pathogensurveillance output folder.
#'
#' @param paths The path to one or more folders that contain pathogensurveillance output.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'psminer')
#' estimated_ani_matrix_path(path)
#'
#' @export
estimated_ani_matrix_path <- function(paths) {
  find_static_file_path(paths, 'sourmash_ani_matrix.csv', 'estimated ANI matrix file')
}

#' Find the software version file path
#'
#' Return the file path to the YAML file with the versions of software used for
#' a given pathogensurveillance output folder.
#'
#' @param paths The path to one or more folders that contain pathogensurveillance output.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'psminer')
#' software_version_path(path)
#'
#' @export
software_version_path <- function(paths) {
  find_static_file_path(paths, 'versions.yml', 'software version file')
}

#' Find the core gene analysis tree paths
#'
#' Return the file paths to the Newick formatted trees produced by the core
#' gene analysis for a given pathogensurveillance output folder.
#'
#' @param paths The path to one or more folders that contain pathogensurveillance output.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'psminer')
#' core_tree_path(path)
#'
#' @export
core_tree_path <- function(paths) {
  find_static_dir_path(paths, 'core_gene_trees', 'core gene tree', file_required = FALSE, dir_required = FALSE)
}

#' Find the considered NCBI reference metadata paths
#'
#' Return the file paths of the TSVs of metadata for references considered by
#' the pipeline for download for a given pathogensurveillance output folder.
#'
#' @param paths The path to one or more folders that contain pathogensurveillance output.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'psminer')
#' considered_ref_meta_path(path)
#'
#' @export
considered_ref_meta_path <- function(paths) {
  find_static_dir_path(paths, 'ncbi_reference_data', 'considered NCBI reference metadata')
}

#' Find the downloaded reference metadata paths
#'
#' Return the file paths of the TSVs of metadata for references selected and
#' downloaded for each sample for a given pathogensurveillance output folder.
#'
#' @param paths The path to one or more folders that contain pathogensurveillance output.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'psminer')
#' selected_ref_meta_path(path)
#'
#' @export
selected_ref_meta_path <- function(paths) {
  find_static_dir_path(paths, 'selected_references', 'selected NCBI reference metadata')
}

#' Find the sendsketch result paths
#'
#' Return the file paths of the sendsketch results for each sample for a given
#' pathogensurveillance output folder.
#'
#' @param paths The path to one or more folders that contain pathogensurveillance output.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'psminer')
#' sendsketch_path(path)
#'
#' @export
sendsketch_path <- function(paths) {
  find_static_dir_path(paths, 'sendsketch', 'sendsketch result')
}

#' Find the SNP alignment paths
#'
#' Return the file paths of the SNP alignments for each reference used in the
#' variant analysis for a given pathogensurveillance output folder.
#'
#' @param paths The path to one or more folders that contain pathogensurveillance output.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'psminer')
#' variant_align_path(path)
#'
#' @export
variant_align_path <- function(paths) {
  find_static_dir_path(paths, 'snp_alignments', 'SNP alignment', dir_required = FALSE, file_required = FALSE)
}

#' Find the SNP tree paths
#'
#' Return the file paths of the SNP tree for each reference used in the
#' variant analysis for a given pathogensurveillance output folder.
#'
#' @param paths The path to one or more folders that contain pathogensurveillance output.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'psminer')
#' variant_tree_path(path)
#'
#' @export
variant_tree_path <- function(paths) {
  find_static_dir_path(paths, 'snp_trees', 'SNP tree', dir_required = FALSE, file_required = FALSE)
}

#' @keywords internal
find_static_file_path <- function(paths, file_name, file_description, file_required = TRUE) {

  group_result_paths <- find_group_result_path(paths)

  find_one <- function(path) {
    output_path <- file.path(path, file_name)
    if (! file.exists(output_path)) {
      if (file_required) {
        stop(
          call. = FALSE,
          'Cannot locate ', file_description, '. It should be located at "',
          output_path, '". Verify that "', path,
          '" is an pathogensurveillance output folder.'
        )
      } else {
        return(character(0))
      }
    }
    return(output_path)
  }
  output <- unlist(lapply(group_result_paths, find_one))
  if (is.null(output)) {
    output <- character(0)
  }
  return(output)
}

#' @keywords internal
find_static_dir_path <- function(paths, dir_name, file_description, dir_required = TRUE, file_required = TRUE, ...) {

  group_result_paths <- find_group_result_path(paths)

  # List files in target directories
  find_one <- function(path) {
    # Check that directory exists
    dir_path <- file.path(path, dir_name)
    if (! file.exists(dir_path)) {
      if (dir_required) {
        stop(
          call. = FALSE,
          'Cannot locate the ', file_description, ' directory. It should be located at "',
          dir_path, '". Verify that "', path,
          '" is an pathogensurveillance output folder.'
        )
      } else {
        return(character(0))
      }
    }
    # List files inside directory
    out_paths <- list.files(dir_path, full.names = TRUE, ...)
    if (file_required && length(out_paths) == 0) {
      stop(
        call. = FALSE,
        'No files found in the ', file_description, ' directory. It should be located at "',
        dir_path, '". Verify that "', path,
        '" is an pathogensurveillance output folder.'
      )
    }
    return(out_paths)
  }
  output <- unlist(lapply(group_result_paths, find_one))
  if (is.null(output)) {
    output <- character(0)
  }
  return(output)
}


#' Find subfolders with report group results
#'
#' For one or more paths, return the paths of subfolders (or the given folders)
#' that are the report group output of a pathogensurviellance run.
#'
#' @keywords internal
find_group_result_path <- function(paths) {
  subdir_paths <- list.dirs(paths)
  is_valid_dir <- function(path) {
    group_id_path <- file.path(path, 'pathogensurveillance_run_info.yml')
    return(file.exists(group_id_path))
  }
  subdir_paths <- subdir_paths[unlist(lapply(subdir_paths, is_valid_dir))]
  return(subdir_paths)
}
