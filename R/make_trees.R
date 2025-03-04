#' Plot a phylogeny
#'
#' @param phylo_path Path to the .treefile
#' @param ids The ids used to identify tips in `phylo_path`. The tree will be subset to these tips.
#' @param labels The labels to print on the plot. Must be the same length as `ids`.
#' @param colors The color or the tips. Must be the same length as `ids`.
#' @param interactive Whether to use an HTML-based interactive format or not (default: TRUE)
#'
#' @export
plot_phylogeny <- function(phylo_path, ids, labels = ids, colors = "#000000", interactive = knitr::is_html_output(), ...) {
  # Read tree file
  tree <- ape::read.tree(phylo_path)

  # Root tree based on sample IDs
  tree <- phytools::midpoint_root(tree)

  # Subset to ids included
  tree <- treeio::tree_subset(tree, )

  if (interactive) {
    plot_phylogeny_phylocanvas(phylo_path = phylo_path, ids = ids, labels = labels, colors = colors)
  } else {
    plot_phylogeny_ggtree(phylo_path = phylo_path, ids = ids, labels = labels, colors = colors)
  }
}
