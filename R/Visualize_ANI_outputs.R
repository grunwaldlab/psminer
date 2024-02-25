#' Make an ANI heatmap and dendrogram.
#'
#' The ANI heatmap is based on approximate ANI similarity matrix output by Sourmash
#'
#' @param ani_matrix Approximate ANI matrix output from Sourmash analysis
#' @param sample_data A tibble/data.frame with the sample metadata
#' @param ref_data A tibble/data.frame with information on references used in analysis
#' @param interactive Whether to use an HTML-based interactive format or not (default: TRUE)
#'
#' @return A heatmap and dendrogram
#'
#' @export
#'
#' @examples
#' make_ani_heatmap(ani_matrix, ref_data, samp_data, interactive=FALSE)

make_ani_heatmap <- function(ani_matrix, ref_data, sample_data, interactive = knitr::is_html_output()) {
  # Rename rows/columns for plotting
  name_key <- c(
    setNames(ref_data$reference_name, ref_data$reference_id),
    setNames(sample_data$sample_name, sample_data$sample_id)
  )
  colnames(ani_matrix) <- name_key[colnames(ani_matrix)]
  rownames(ani_matrix) <- name_key[rownames(ani_matrix)]

  if (interactive) {
    heatmap_ani <- heatmaply(ani_matrix,
                             fontsize_row = 8, fontsize_col = 8, width = 1000, height = 800)
  } else {
    heatmap_ani <- pheatmap(ani_matrix, show_rownames = TRUE, labels_row = colnames(ani_matrix))
  }
  return(heatmap_ani)
}
