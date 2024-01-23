#' Make an ANI heatmap and dendrogram.
#'
#' The ANI heatmap will either be static when PDF is rendered or interactive a NJ tree from Sourmash approximate ANI values
#'
#' @param sourmash_ani_matrix Approximate ANI matrix output from Sourmash analysis
#' @param sample_data A tibble/data.frame with the sample metadata
#' @param reference_data A tibble/data/frame with information on references used in analysis
#' @param interactive Whether to use an HTML-based interactive format or not (default: TRUE)
#' @return A heatmap and dendrogram
#' @export
#' @examples
#' make_ani_heatmap(ani_matrix, ref_data, samp_data, interactive=FALSE)

make_ani_heatmap <- function(sourmash_ani_matrix, reference_data, sample_data, interactive = TRUE) {
  ani_matrix_format <-as.matrix(sourmash_ani_matrix)
  colnames(ani_matrix_format) <- convert_id(colnames(ani_matrix_format))
  rownames(ani_matrix_format) <- convert_id(rownames(ani_matrix_format))
  assembly_entries <- rownames(ani_matrix_format)[grepl("_assembly$", rownames(ani_matrix_format))]

  name_key <- c(
    setNames(reference_data$Organism, convert_id(reference_data$LastMajorReleaseAccession)),
    setNames(sample_data$sample, convert_id(sample_data$sample)),
    setNames(assembly_entries, assembly_entries)
  )

  # Relabel the tree tips
  colnames(ani_matrix_format) <- name_key[colnames(ani_matrix_format)]
  rownames(ani_matrix_format) <- name_key[rownames(ani_matrix_format)]

  if (interactive) {
    #Trying to find good solution for interactive plot, but based on complicated renaming of columns and rows, I am just putting this as place holder
    heatmap_ani <- pheatmap(ani_matrix_format, show_rownames = T, labels_row =colnames(ani_matrix_format))
  } else {
    heatmap_ani <- pheatmap(ani_matrix_format, show_rownames = T, labels_row =colnames(ani_matrix_format))
  }
}
