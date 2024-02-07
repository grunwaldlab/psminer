#' Make core gene phylogeny
#'
#' The tree will either be made using ggtree, if output is a PDF, or phylo
#'
#' @param core_phylogeny_path
#' @param sample_data
#' @param formatted_ani_matrix
#' @param reference_data
#' @param interactive
#'
#' @return Core gene phylogeny
#'
#' @export
make_phylogeny <- function(core_phylogeny_path, sample_data, formatted_ani_matrix, reference_data, interactive = TRUE) {
  core_tree <- ape::read.tree(core_phylogeny_path)
  sample_ids <- core_tree$tip.label[core_tree$tip.label %in% sample_data$modified_id]

  # Root tree based on sample IDs
  group_ani <- formatted_ani_matrix[rownames(formatted_ani_matrix) %in% core_tree$tip.label, colnames(formatted_ani_matrix) %in% core_tree$tip.label]
  core_tree <- root(core_tree, names(which.min(colMeans(group_ani[sample_ids, ]))))

  # Set tip labels to taxon names
  name_key <- set_names(c(reference_data$display_name_shorter, sample_data$modified_id), c(reference_data$LastMajorReleaseAccession, sample_data$modified_id))
  core_tree$tip.label <- name_key[core_tree$tip.label]

  if (interactive) {
    # Create phylocanvas for interactive visualization
    phycanv <- phylocanvas(core_tree, treetype = "rectangular", alignlabels = TRUE, showscalebar = TRUE, width = "100%")

    # Style nodes for sample IDs
    for (x in name_key[sample_ids]) {
      phycanv <- style_node(phycanv, x, labelcolor = "green", labeltextsize = 30)
    }
    phycanv
  } else {
    print("in progress")
    #ggtree_core_tree <- ggtree(core_tree)

    #ggtree_core_tree <- ggtree_core_tree +
      #geom_tiplab(aes(label = name_key[core_tree$tip.label], color = ifelse(core_tree$tip.label %in% name_key[sample_ids], "green", "black")), size = 3, hjust = 0, vjust = 0.5, fontface = "bold") +
      #scale_color_identity() +
      #theme(legend.position = "none")

    #ggtree_core_tree
  }
}
