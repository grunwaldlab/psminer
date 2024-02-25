#' Parse sample metadata files
#'
#'
#' @param sample_data_path The file path to the 'sample_data.csv' file
#' @param assigned_refs_path The file path to the 'assigned_refs.csv' file
#' @param group The name of the report group for specific set of samples.
#'
#' @return `data.frame` of metadata for samples, with one row for each sample in
#'   the group
#'
#' @export
parse_sample_meta <- function(sample_data_path, assigned_refs_path, group) {
  sample_data <- read.csv(sample_data_path, check.names = FALSE)
  assigned_refs <- read.csv(assigned_refs_path, col.names = c('sample_id', 'reference_id'), header = FALSE)
  # Subset to samples in the report group
  group_data <- strsplit(sample_data$report_group, split = ";")
  sample_data <- sample_data[map_lgl(group_data, function(x) group %in% x), ]
  # Add in which reference was used in variant calling
  sample_data[match(assigned_refs$sample_id, sample_data$sample_id), 'reference_id'] <- assigned_refs$reference_id
  return(sample_data)
}


#' Parse reference metadata files
#'
#'
#' @param ref_data_path The folder with files with the data on references
#'   assigned to each sample. The individual files are named by sample ID.
#' @param assigned_refs_path The file path to the 'assigned_refs.csv' file
#' @param sample_data_path The file path to the 'sample_data.csv' file
#' @param group The name of the report group for specific set of samples.
#'
#' @return Reformatted reference metadata table containing both refseq-defined
#'   reference assemblies and user-defined assemblies
#'
#' @export

parse_ref_meta <- function(ref_data_path, assigned_refs_path, sample_data_path, group) {
  sample_data <- read.csv(sample_data_path, check.names = FALSE)
  group_data <- strsplit(sample_data$report_group, split = ";")
  sample_data <- sample_data[map_lgl(group_data, function(x) group %in% x), ] # Subset to samples in the report group
  assigned_refs <- read.csv(assigned_refs_path, col.names = c('sample_id', 'reference_id'), header = FALSE)

  # Get data for the references selected by the pipeline
  ref_file_names <- list.files(ref_data_path)
  ref_data_per_sample <- lapply(ref_file_names, function(path) {
    sample_id <- gsub(path, pattern = '\\.tsv$', replacement = '')
    output <- read.csv(file.path(ref_data_path, path), sep = '\t')
    return(cbind(
      used_in_phylo_by = sample_id,
      source = 'pipeline',
      reference_name = output$Organism,
      output
    ))
  })
  ref_data <- do.call(rbind, ref_data_per_sample)

  # Combine sample IDs for each unique reference ID into a single row
  ref_data <- do.call(rbind, lapply(unique(ref_data$reference_id), function(id) {
    subset <- ref_data[ref_data$reference_id == id, ]
    subset$used_in_phylo_by <- paste0(unique(subset$used_in_phylo_by), collapse = ';')
    return(subset[1, ])
  }))

  # Add user-defined references
  user_ref_ids <- unique(sample_data$reference_id)
  user_ref_ids <- user_ref_ids[user_ref_ids != '' & !is.na(user_ref_ids)]
  user_ref_data <- ref_data[NA, ][seq_along(user_ref_ids), ] # Hack to make a data.frame with all NAs
  user_ref_data$reference_id <- user_ref_ids
  user_ref_data$reference_name <- sample_data$reference_name[match(user_ref_ids, sample_data$reference_name)]
  user_ref_data$source <- 'user'
  ref_data <- rbind(ref_data, user_ref_data)
  row.names(ref_data) <- NULL

  # Add in which reference was used in variant calling
  samps_per_ref <- split(assigned_refs$sample_id, assigned_refs$reference_id)
  samp_ref_key <- unlist(lapply(samps_per_ref, paste0, collapse = ';'))
  ref_data$used_in_var_by <- samp_ref_key[ref_data$reference_id]

  return(ref_data)
}

#' Parse and rename ANI matrix names for downstream steps
#'
#'
#' @param ani_matrix_path the file path to 'ani_matrix.csv' file
#'
#' @return ANI matrix with column and row names that are compatible with other data inputs
#'
#' @export
parse_ani_matrix <- function(ani_matrix_path) {
  ani_matrix <- read.csv(ani_matrix_path, check.names = FALSE)
  rownames(ani_matrix) <- colnames(ani_matrix)
  return(ani_matrix)
}


