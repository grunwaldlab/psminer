#' Parse sample metadata files
#'
#'
#' @param path the file path to 'samp_data.csv' file
#' @param group parsed information on report group for specific set of samples
#'
#' @return samp_data dataframe
#'
#' @export

parse_sample_meta <- function(path, group) {
  all_samp_data <- read_csv(path, show_col_types = FALSE)
  all_samp_data$modified_id <- gsub(all_samp_data$sample, pattern = "-", replacement = "_", fixed = TRUE) #TODO: check that this is needed with new data format
  group_data <- strsplit(all_samp_data$report_group, split = ";")
  samp_data <- all_samp_data[map_lgl(group_data, function(x) group %in% x), ]
  return(samp_data)
}


#' Parse reference metadata files
#'
#'
#' @param refseq_path
#' @param user_ref_path
#'
#' @return Reformatted reference metadata table containing both refseq-defined reference assemblies and user-defined assemblies
#'
#' @export

parse_ref_meta <- function(refseq_path, user_ref_path) {
  ref_data <- read_tsv(ref_data_path , col_types = 'dcccccccccccccccddc')
  refs <- strsplit(read_lines(ref_ids_path), split = ';', fixed = TRUE)[[1]]
  ref_data$Origin = "refseq"
  ref_data$Reference_id = ref_data$LastMajorReleaseAccession
  ref_data$Display_name = ref_data$Organism

  new_reference_ids <- unique(samp_data$reference_id[!is.na(samp_data$reference_id) & !samp_data$reference_id %in% ref_data$Reference_id])

  if (length(new_reference_ids) > 0) {
    new_rows <- data.frame(
      display_name = new_reference_ids,
      reference_id = new_reference_ids,
      origin = "user"
    )
  }
    ref_data <- bind_rows(new_rows, ref_data)
    return(ref_data)
}
