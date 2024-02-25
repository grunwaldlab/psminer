#' Read SNP trees
#' 
#' Reads .treefile files produced by variant calling into R.
#' 
#' @param snp_tree_paths One or more paths to .treefiles
#' @param group The name of the report group this tree is for
#' 
#' @export
parse_snp_trees <- function(snp_tree_paths, group) {
    lapply(seq_along(snp_tree_paths), function(i) {
        file_name <- sub(basename(snp_tree_paths[i]), pattern = '\\.treefile$', replacement = '')
        ref_id <- sub(file_name, pattern = paste0('^', group, '_'), replacement = '')
        tree <- ape::read.tree(snp_tree_paths[i])
        tree$tip.label <- sub(tree$tip.label, pattern = paste0('^', ref_id, '_'), replacement = '')
        tree$tip.label[tree$tip.label == "REF"] <- ref_id
        return(tree)
    })
}
