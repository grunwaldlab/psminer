#' Make core gene phylogeny
#'
#' The tree will either be made using ggtree, if output is a PDF, or phylo
#'
#' @param core_phylogeny_path
#' @param sample_data
#' @param sourmash_sourmash_ani_matrix
#' @param reference_data
#' @param interactive
#'
#' @return Core gene phylogeny
#'
#' @export
make_phylogeny <- function(core_phylogeny_path, sample_data, sourmash_ani_matrix, reference_data, interactive = TRUE) {
  # Helper function to convert IDs
  convert_id <- function(ids) gsub(ids, pattern = "[.-]", replacement = "_")

  # Read core tree
  core_tree <- ape::read.tree(core_phylogeny_path)

  # Identify sample IDs in the tree
  sample_ids <- core_tree$tip.label[core_tree$tip.label %in% convert_id(sample_data$sample)]

  # Root tree based on sample IDs
  colnames(sourmash_ani_matrix) <- convert_id(colnames(sourmash_ani_matrix))
  rownames(sourmash_ani_matrix) <- colnames(sourmash_ani_matrix)
  group_ani <- sourmash_ani_matrix[rownames(sourmash_ani_matrix) %in% core_tree$tip.label, colnames(sourmash_ani_matrix) %in% core_tree$tip.label]
  core_tree <- root(core_tree, names(which.min(colMeans(group_ani[sample_ids, ]))))

  # Set tip labels to taxon names
  name_key <- set_names(c(reference_data$Organism, sample_data$sample), c(convert_id(reference_data$LastMajorReleaseAccession), convert_id(sample_data$sample)))
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
    print("In progress")
  }
}
