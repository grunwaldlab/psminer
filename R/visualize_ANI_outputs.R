#' Make an ANI heatmap and dendrogram.
#'
#' The ANI heatmap is based on approximate ANI similarity matrix output by Sourmash
#'
#' @param ani_matrix Approximate ANI matrix output from Sourmash analysis
#' @param sample_data A tibble/data.frame with the sample metadata
#' @param ref_data A tibble/data.frame with information on references used in analysis
#' @param interactive Whether to use an HTML-based interactive format or not (default: TRUE)
#' @param height The height in pixels. If not `interactive`, this is divided by `dpi` to convert it to inches.
#' @param width The width in pixels. If not `interactive`, this is divided by `dpi` to convert it to inches.
#' @param dpi How pixels are converted to inches
#'
#' @return A heatmap and dendrogram
#'
#' @export
#'
#' @examples
#' make_ani_heatmap(ani_matrix, ref_data, samp_data, interactive=FALSE)

make_ani_heatmap <- function(ani_matrix, ref_data, sample_data, interactive = FALSE, height = 1000, width = 1000, dpi = 100) {
  # Rename rows/columns for plotting
  name_key <- c(
    stats::setNames(ref_data$ref_name, ref_data$ref_id),
    stats::setNames(sample_data$name, sample_data$sample_id)
  )
  name_key <- stats::setNames(make.unique(name_key, sep = ' '), names(name_key))
  colnames(ani_matrix) <- name_key[colnames(ani_matrix)]
  rownames(ani_matrix) <- name_key[rownames(ani_matrix)]
  if (interactive) {
    heatmap_ani <- heatmaply::heatmaply(ani_matrix, fontsize_row = 8, fontsize_col = 8, width = width, height = height)
  } else {
    heatmap_ani <- pheatmap::pheatmap(ani_matrix, show_rownames = TRUE, labels_row = colnames(ani_matrix), width = width / dpi, height = height / dpi)
  }
  return(heatmap_ani)
}
