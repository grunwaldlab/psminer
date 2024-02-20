#' Make Minimum spanning network
#'
#' @param tree_path
#' @param snp_align_path
#' @param sample_data
#' @param population
#' @param interactive
#' @param snp_threshold
#' @param show_MLG_table
#' @return minimum spanning network
#'
#' @export


make_MSN <- function(tree_path, snp_alignment_path, sample_data, population=NULL, interactive = TRUE, snp_threshold=NULL, show_MLG_table=FALSE) {
  snp_trees <- ape::read.tree(tree_path)
  snp_alignment <- ape::read.dna(snp_alignment_path, format =  "fasta")
  snp_aln.gi <- DNAbin2genind(snp_alignment)
  snp_aln.gi <- snp_aln.gi[indNames(snp_aln.gi) != "REF"]

  genind_names <- indNames(snp_aln.gi)
  cleaned_names <- sub(".*assembly_", "", genind_names)
  indNames(snp_aln.gi) <- cleaned_names

  mat <- match(indNames(snp_aln.gi), samp_data$sample)
  samp_data <- samp_data[mat, ]
  snp_genclone <- as.genclone(snp_aln.gi)

  if (is.null(snp_threshold)) {
    snpdist_stats <- filter_stats(snp_genclone)
    average_thresh <- cutoff_predictor(snpdist_stats$average$THRESHOLDS)
    mlg.filter(snp_genclone, distance = bitwise.dist, percent = FALSE, threshold=average_thresh)

    } else {
    mlg.filter(snp_genclone, distance = bitwise.dist, percent = FALSE, threshold=snp_threshold)
    }

    if (!is.null(population) && population %in% names(samp_data)) {
      # Extract population from samp_data based on the specified column
      user_factor <- samp_data[[population]]
      node_color <- as.factor(ifelse(is.na(user_factor), "Unknown", user_factor))
      myColors <- rainbow(length(unique(node_color)))
      names(myColors) <- levels(node_color)
      num_columns <- ncol(samp_data)
      strata(snp_genclone) <- cbind(samp_data[, c(1:num_columns)], node_color)
      names(strata(snp_genclone))[num_columns + 1] <- "color_node_by"
      setPop(snp_genclone) <- ~color_node_by

      ms.loc <- poppr.msn(snp_genclone,
                          distmat = bitwise.dist(snp_genclone, percent = FALSE),
                          include.ties = TRUE,
                          showplot = FALSE)
      the_edges <- igraph::E(ms.loc$graph)$weight
      edges <- as.list(the_edges)

      set.seed(8)
      plot_poppr_msn(
        snp_genclone,
        poppr_msn = ms.loc,
        palette = myColors,
        mlg = FALSE,
        quantiles = FALSE,
        wscale = FALSE,
        inds = "None",
        layfun = igraph::layout_with_lgl,
        edge.label = the_edges,
        edge.label.font = 2,
        edge.label.cex = 1,
        edge.label.family = "Helvetica",
        edge.label.color = "darkslateblue")

    } else {
      print("in progress")
    }

  if (interactive) {
    print("in progress")
  }
}




