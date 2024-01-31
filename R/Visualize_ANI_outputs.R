#' Make an ANI heatmap and dendrogram.
#'
#' The ANI heatmap will either be static when PDF is rendered or interactive a NJ tree from Sourmash approximate ANI values
#'
#' @param sourmash_ani_matrix Approximate ANI matrix output from Sourmash analysis
#' @param sample_data A tibble/data.frame with the sample metadata
#' @param reference_data A tibble/data/frame with information on references used in analysis
#' @param interactive Whether to use an HTML-based interactive format or not (default: TRUE)
#'
#' @return A heatmap and dendrogram
#'
#' @export
#'
#' @examples
#' make_ani_heatmap(ani_matrix, ref_data, samp_data, interactive=FALSE)

make_ani_heatmap <- function(sourmash_ani_matrix, reference_data, sample_data, interactive = TRUE) {

  assembly_entries <- rownames(ani_matrix_format)[grepl("_assembly$", rownames(ani_matrix_format))]

  name_key <- c(
    setNames(reference_data$Organism, reference_data$LastMajorReleaseAccession),
    setNames(sample_data$modified_id, sample_data$modified_id),
    setNames(assembly_entries, assembly_entries)
  )

  colnames(ani_matrix_format) <- name_key[colnames(ani_matrix_format)]
  rownames(ani_matrix_format) <- name_key[rownames(ani_matrix_format)]

  if (interactive) {
    # Trying to find a good solution for an interactive plot, but based on complicated renaming of columns and rows,
    # Still figuring this out. Defaulting to pheatmap for both for now.
    heatmap_ani <- pheatmap(ani_matrix_format, show_rownames = TRUE, labels_row = colnames(ani_matrix_format))
  } else {
    heatmap_ani <- pheatmap(ani_matrix_format, show_rownames = TRUE, labels_row = colnames(ani_matrix_format))
  }
}
