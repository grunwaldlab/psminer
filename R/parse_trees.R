#' Read SNP trees
#'
#' Reads .treefile files produced by variant calling into R.
#'
#' @param snp_tree_paths One or more paths to .treefiles
#' @param group The name of the report group this tree is for
#'
#' @export
parse_snp_trees <- function(snp_tree_paths, group) {
  file_name <- sub(basename(snp_tree_paths), pattern = '\\.treefile$', replacement = '')
  ref_id <- sub(file_name, pattern = paste0('^', group, '_'), replacement = '')
  output <- lapply(seq_along(snp_tree_paths), function(i) {
    tree <- ape::read.tree(snp_tree_paths[i])
    tree$tip.label <- sub(tree$tip.label, pattern = paste0('^', ref_id[i], '_'), replacement = '')
    tree$tip.label[tree$tip.label == "REF"] <- ref_id[i]
    return(tree)
  })
  names(output) <- ref_id
  return(output)
}
