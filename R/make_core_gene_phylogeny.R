#' Make core gene phylogeny
#'
#' The tree will either be made using ggtree, if output is a PDF, or phylo
#'
#' @param input Either a file path to a .treefile or a "phylo" object
#' @param sample_data
#' @param ani_matrix
#' @param ref_data
#' @param interactive
#'
#' @return Core gene phylogeny
#'
#' @export
make_phylogeny <- function(input, sample_data, ref_data, interactive = knitr::is_html_output()) {
  # Check if the input is a path or a parsed phylo object
  if (class(input) != "phylo") {
    tree <- ape::read.tree(input)
  } else {
    tree <- input
  }

  tip_ids <- tree$tip.label
  tree <- phangorn::midpoint(tree) # Root tree
  name_key <- set_names(c(ref_data$reference_name, sample_data$sample_name),
                        c(ref_data$reference_id, sample_data$sample_id))

  if (interactive) {
    # Convert characters that are not allowed in the format used by phylocanvas
    name_key <- gsub(name_key, pattern = ')', replacement = ']', fixed = TRUE)
    name_key <- gsub(name_key, pattern = '(', replacement = '[', fixed = TRUE)
    name_key <- gsub(name_key, pattern = ':', replacement = '-', fixed = TRUE)
    name_key <- gsub(name_key, pattern = ',', replacement = '.', fixed = TRUE)
    tree$tip.label <- name_key[tree$tip.label]

    # Create phylocanvas for interactive visualization
    phycanv <- phylocanvas(tree, treetype = "rectangular", alignlabels = TRUE, showscalebar = TRUE, width = "100%")

    # Remove underscores from tip labels (phylocanvas converts spaces to underscores)
    phycanv$x$tree <- gsub(phycanv$x$tree, pattern = '_', replacement = ' ', fixed = TRUE)

    # Style nodes for sample IDs
    sample_ids <- tip_ids[tip_ids %in% sample_data$sample_id]
    for (x in name_key[sample_ids]) { # NOTE: using the names instead of IDs like this could cause bugs if the name is not unique
      phycanv <- style_node(phycanv, x, labelcolor = "green", labeltextsize = 30)
    }
    return(phycanv)
  } else {
    tree <- groupOTU(tree, .node = tip_ids[tip_ids %in% sample_data$sample_id])
    plotted_tree <- ggtree(tree) +
      geom_tiplab(aes(color = group), show.legend = FALSE) +
      scale_color_manual(values = c("black", "green"), name = '')
    return(plotted_tree)
  }
}

