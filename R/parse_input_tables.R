#' Parse sample metadata files
#'
#'
#' @param path the file path to 'sample_data.csv' file
#' @param group parsed information on report group for specific set of samples
#'
#' @return sample_data dataframe
#'
#' @export

parse_sample_meta <- function(path, group) {
  all_sample_data <- read_csv(path, show_col_types = FALSE)
  all_sample_data$modified_id <- gsub(all_sample_data$sample, pattern = "-", replacement = "_", fixed = TRUE) #TODO: check that this is needed with new data format
  group_data <- strsplit(all_sample_data$report_group, split = ";")
  sample_data <- all_sample_data[map_lgl(group_data, function(x) group %in% x), ]
  return(sample_data)
}


#' Parse reference metadata files
#'
#'
#' @param reference_data_path
#' @param ref_ids_path
#' @param sample_data
#'
#' @return Reformatted reference metadata table containing both refseq-defined reference assemblies and user-defined assemblies
#'
#' @export

parse_ref_meta <- function(reference_data_path, ref_ids_path, sample_data) {
  convert_id <- function(ids) gsub(ids, pattern = "[.-]", replacement = "_")

  reference_data <- read_tsv(reference_data_path , col_types = 'dcccccccccccccccddc')
  refs <- strsplit(read_lines(ref_ids_path), split = ';', fixed = TRUE)[[1]]

  reference_data$LastMajorReleaseAccession<-convert_id(reference_data$LastMajorReleaseAccession)

  reference_data$origin = "refseq"
  reference_data$reference_id = reference_data$LastMajorReleaseAccession
  reference_data$display_name = reference_data$Organism

  new_reference_ids <- unique(sample_data$reference_id[!is.na(sample_data$reference_id) & !sample_data$reference_id %in% reference_data$reference_id])

  if (length(new_reference_ids) > 0) {
    new_rows <- data.frame(
      display_name = new_reference_ids,
      reference_id = new_reference_ids,
      origin = "user"
    )
    reference_data <- bind_rows(new_rows,reference_data)
  }

  return(reference_data)
}
