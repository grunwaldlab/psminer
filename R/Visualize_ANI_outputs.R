#' Make an ANI heatmap and dendrogram.
#'
#' The ANI heatmap is based on approximate ANI similarity matrix output by Sourmash
#'
#' @param formatted_ani_matri Approximate ANI matrix output from Sourmash analysis
#' @param sample_data A tibble/data.frame with the sample metadata
#' @param reference_data A tibble/data.frame with information on references used in analysis
#' @param interactive Whether to use an HTML-based interactive format or not (default: TRUE)
#'
#' @return A heatmap and dendrogram
#'
#' @export
#'
#' @examples
#' make_ani_heatmap(ani_matrix, ref_data, samp_data, interactive=FALSE)

make_ani_heatmap <- function(formatted_ani_matrix, reference_data, sample_data, interactive = TRUE) {

  assembly_entries <- rownames(formatted_ani_matrix)[grepl("_assembly$", rownames(formatted_ani_matrix))]

  name_key <- c(
    setNames(reference_data$display_name_shorter, reference_data$LastMajorReleaseAccession),
    setNames(sample_data$modified_id, sample_data$modified_id),
    setNames(reference_data$display_name_shorter, assembly_entries)
  )

  colnames(formatted_ani_matrix) <- name_key[colnames(formatted_ani_matrix)]
  rownames(formatted_ani_matrix) <- name_key[rownames(formatted_ani_matrix)]

  if (interactive) {
    heatmap_ani <- heatmaply(formatted_ani_matrix,
                             fontsize_row = 8, fontsize_col = 8, width = 1000, height = 800)
  } else {
    heatmap_ani <- pheatmap(formatted_ani_matrix, show_rownames = TRUE, labels_row = colnames(formatted_ani_matrix))
  }
  return(heatmap_ani)
}
