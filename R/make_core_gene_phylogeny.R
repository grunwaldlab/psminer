#' Make core gene phylogeny
#'
#' The tree will either be made using ggtree, if output is a PDF, or phylo
#'
#' @param core_phylo_path
#' @param sample_data
#' @param ani_matrix
#' @param ref_data
#' @param interactive
#'
#' @return Core gene phylogeny
#'
#' @export
make_phylogeny <- function(core_phylo_path, sample_data, ref_data, interactive = TRUE) {
  core_tree <- ape::read.tree(core_phylo_path)
  tip_ids <- core_tree$tip.label
  core_tree <- phangorn::midpoint(core_tree) # Root tree

  # Set tip labels to taxon names
  name_key <- set_names(c(ref_data$reference_name, sample_data$sample_name),
                        c(ref_data$reference_id, sample_data$sample_id))
  core_tree$tip.label <- name_key[core_tree$tip.label]

  if (interactive) {
    # Create phylocanvas for interactive visualization
    phycanv <- phylocanvas(core_tree, treetype = "rectangular", alignlabels = TRUE, showscalebar = TRUE, width = "100%")

    # Style nodes for sample IDs
    sample_ids <- tip_ids[tip_ids %in% sample_data$sample_id]
    for (x in name_key[sample_ids]) {
      phycanv <- style_node(phycanv, x, labelcolor = "green", labeltextsize = 30)
    }
    phycanv
  } else {
    print("in progress")
  }
}



#' Plot a phylogeny
#'
#' @param phylo_path Path to the .treefile
#' @param ids The ids used to identify tips in `phylo_path`. The tree will be subset to these tips.
#' @param labels The labels to print on the plot. Must be the same length as `ids`.
#' @param colors The color or the tips. Must be the same length as `ids`.
#' @param interactive Whether to use an HTML-based interactive format or not (default: TRUE)
#'
#' @export
plot_phylogeny <- function(phylo_path, ids, labels = ids, colors = "#000000", interactive = TRUE, ...) {
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
