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
#' @param ref_data_path
#' @param ref_ids_path
#' @param sample_data
#'
#' @return Reformatted reference metadata table containing both refseq-defined reference assemblies and user-defined assemblies
#'
#' @export

parse_ref_meta <- function(ref_data_path, ref_ids_path, sample_data) {
  convert_id <- function(ids) gsub(ids, pattern = "[.-]", replacement = "_")

  ref_data <- read_tsv(ref_data_path , col_types = 'dcccccccccccccccddc')
  refs <- strsplit(read_lines(ref_ids_path), split = ';', fixed = TRUE)[[1]]

  reference_data$LastMajorReleaseAccession<-convert_id(reference_data$LastMajorReleaseAccession)

  ref_data$origin = "refseq"
  ref_data$reference_id = ref_data$LastMajorReleaseAccession
  ref_data$display_name = ref_data$Organism

  new_reference_ids <- unique(sample_data$reference_id[!is.na(sample_data$reference_id) & !sample_data$reference_id %in% ref_data$reference_id])

  if (length(new_reference_ids) > 0) {
    new_rows <- data.frame(
      display_name = new_reference_ids,
      reference_id = new_reference_ids,
      origin = "user"
    )
    ref_data <- bind_rows(new_rows,ref_data)
  }

  return(ref_data)
}
