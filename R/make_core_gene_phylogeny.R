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
#'
#' @examples
#'
make_phylogeny <- function(core_phylogeny_path, sample_data, sourmash_sourmash_ani_matrix, reference_data, interactive = TRUE) {
  #This helpfer function is used quite a bit, so perhaps it is worth pulling out of individual fxns?
  convert_id <- function(ids) {
    gsub(ids, pattern = "[.-]", replacement = "_")
  }

  core_tree <- ape::read.tree(core_phylogeny_path)

  # Identify which tips are samples and references
  sample_ids <- core_tree$tip.label[core_tree$tip.label %in% convert_id(sample_data$sample)]

  # Root tree
  colnames(sourmash_ani_matrix) <- convert_id(colnames(sourmash_ani_matrix))
  rownames(sourmash_ani_matrix) <- colnames(sourmash_ani_matrix)
  group_ani <- sourmash_ani_matrix[rownames(sourmash_ani_matrix) %in% core_tree$tip.label, colnames(sourmash_ani_matrix) %in% core_tree$tip.label]
  core_tree <- root(core_tree, names(which.min(colMeans(group_ani[sample_ids, ]))))

  # Set tip labels to taxon names for reference sequences
  # TODO-we need to generalize
  name_key <- c(
    reference_data$Organism,
    sample_data$sample
  )
  names(name_key) <- c(
    convert_id(reference_data$LastMajorReleaseAccession),
    convert_id(sample_data$sample)
  )
  core_tree$tip.label <- name_key[core_tree$tip.label]


  if (interactive) {
    phycanv <- phylocanvas(core_tree, treetype = "rectangular", alignlabels = T, showscalebar = T, width = "100%")
    for (x in name_key[sample_ids]) {
      phycanv <- style_node(phycanv, x, labelcolor = "green", labeltextsize = 30)
  } else {
    print("In progress")
    }
  }
}
