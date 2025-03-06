#' Plot core gene phylogeny
#'
#' Plot the core gene phylogenies present in the output of a pathogensurveillance
#' run.
#'
#' @param input The path to one or more folders that contain
#'   pathogensurveillance output or paths to tree files.
#' @param collapse_by_tax A [base::character()] vector of taxonomic
#'   classifications, each delimited by `;`, and named by sample or reference
#'   ids present in `sample_meta` or `ref_meta`. These are used to provide a
#'   taxonomic tree that functions as a backbone to combine multiple trees into
#'   a single one. (Default: return a list of plots instead of a single plot)
#' @param interactive Whether to use an HTML-based interactive format or not
#'   (default: TRUE)
#'
#' @return  A list of plots, unless `collapse_by_tax` is used, in which case a single plot is returned.
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'psminer')
#' core_tree_plot(path)
#'
#' @export
core_tree_plot <- function(input, collapse_by_tax = NULL, interactive = FALSE) {
  generalized_tree_plot(input, core_tree_parsed, collapse_by_tax = collapse_by_tax, interactive = interactive)
}

#' Plot busco gene phylogeny
#'
#' Plot the busco gene phylogenies present in the output of a pathogensurveillance
#' run.
#'
#' @param input The path to one or more folders that contain
#'   pathogensurveillance output or paths to tree files.
#' @param collapse_by_tax A [base::character()] vector of taxonomic
#'   classifications, each delimited by `;`, and named by sample or reference
#'   ids present in `sample_meta` or `ref_meta`. These are used to provide a
#'   taxonomic tree that functions as a backbone to combine multiple trees into
#'   a single one. (Default: return a list of plots instead of a single plot)
#' @param interactive Whether to use an HTML-based interactive format or not
#'   (default: TRUE)
#'
#' @return  A list of plots, unless `collapse_by_tax` is useded, in which case a single plot is returned.
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'psminer')
#' busco_tree_plot(path)
#'
#' @export
busco_tree_plot <- function(input, collapse_by_tax = NULL, interactive = FALSE) {
  generalized_tree_plot(input, busco_tree_parsed, collapse_by_tax = collapse_by_tax, interactive = interactive)
}

#' Plot multigene phylogeny
#'
#' Plot the any multigene phylogenies present in the output of a pathogensurveillance
#' run. This includes core gene phylogenies and busco gene phylogenes.
#'
#' @param input The path to one or more folders that contain
#'   pathogensurveillance output or paths to tree files.
#' @param collapse_by_tax A [base::character()] vector of taxonomic
#'   classifications, each delimited by `;`, and named by sample or reference
#'   ids present in `sample_meta` or `ref_meta`. These are used to provide a
#'   taxonomic tree that functions as a backbone to combine multiple trees into
#'   a single one. (Default: return a list of plots instead of a single plot)
#' @param interactive Whether to use an HTML-based interactive format or not
#'   (default: TRUE)
#'
#' @return  A list of plots, unless `collapse_by_tax` is used, in which case a single plot is returned.
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'psminer')
#' multigene_tree_plot(path)
#'
#' @export
multigene_tree_plot <- function(input, collapse_by_tax = NULL, interactive = FALSE) {
  generalized_tree_plot(input, multigene_tree_parsed, collapse_by_tax = collapse_by_tax, interactive = interactive)
}

#' Plot generic phylogeny
#'
#' Plot phylogenies present in the output of a pathogensurveillance run.
#'
#' @param input The path to one or more folders that contain
#'   pathogensurveillance output or paths to tree files.
#' @param parser A function that takes directory paths as input and returns
#'   parsed trees found in those directories, such as
#'   [psminer::core_tree_parsed()].
#' @param collapse_by_tax A [base::character()] vector of taxonomic
#'   classifications, each delimited by `;`, and named by sample or reference
#'   ids present in `sample_meta` or `ref_meta`. These are used to provide a
#'   taxonomic tree that functions as a backbone to combine multiple trees into
#'   a single one. (Default: return a list of plots instead of a single plot)
#' @param interactive Whether to use an HTML-based interactive format or not
#'   (default: TRUE)
#'
#' @return  A list of plots, unless `collapse_by_tax` is used, in which case a
#'   single plot is returned.
#'
#' @keywords internal
generalized_tree_plot <- function(input, parser, collapse_by_tax = NULL, interactive = FALSE) {
  # If no trees are found, return an empty list
  trees <-  parser(input)
  if (length(trees) == 0) {
    return(list())
  }

  # Find and parse needed data
  sample_meta <- sample_meta_parsed(input)
  ref_meta <- ref_meta_parsed(input)
  sendsketch <- sendsketch_taxonomy_data_parsed(input, only_best = TRUE, only_shared = TRUE)

  # Find which columns are used to provide colors to the trees, if any
  ids_in_trees <- unique(unlist(lapply(trees, function(t) t$tip.label)))
  color_by_cols <- unique(unlist(strsplit(sample_meta$color_by[sample_meta$sample_id %in% ids_in_trees], split = ';')))
  color_by_cols <- color_by_cols[!is.na(color_by_cols)]
  color_by_col_names <- c(color_by_cols, 'Default')
  color_by_cols <- c(as.list(color_by_cols), list(NULL))  # NULL ensures that the default color scheme is also used

  # Plot one tree for each color_by column
  tree_plots <- lapply(color_by_cols, function(color_by) {
    plot_phylogeny(
      trees,
      sample_meta,
      ref_meta,
      color_by,
      sendsketch, # NOTE: A better source for the taxonomy should be found, perhaps based on reference NCBI taxon IDs.
      interactive = interactive
    )
  })
  names(tree_plots) <- color_by_col_names
  return(tree_plots)
}

#' Plot SNP trees
#'
#' Plot the SNP trees from the variant analysis present in the output of a pathogensurveillance
#' run.
#'
#' @param input The path to one or more folders that contain
#'   pathogensurveillance output or paths to tree files.
#' @param collapse_by_tax A [base::character()] vector of taxonomic
#'   classifications, each delimited by `;`, and named by sample or reference
#'   ids present in `sample_meta` or `ref_meta`. These are used to provide a
#'   taxonomic tree that functions as a backbone to combine multiple trees into
#'   a single one. (Default: return a list of plots instead of a single plot)
#' @param interactive Whether to use an HTML-based interactive format or not
#'   (default: TRUE)
#'
#' @return  A list of plots, unless `collapse_by_tax` is used, in which case a single plot is returned.
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'psminer')
#' variant_tree_plot(path)
#'
#' @export
variant_tree_plot <- function(input, collapse_by_tax = NULL, interactive = FALSE) {
  # If no trees are found, return an empty list
  trees <-  variant_tree_parsed(input)
  if (length(trees) == 0) {
    return(list())
  }

  # Find and parse needed data
  sample_meta <- sample_meta_parsed(input)
  ref_meta <- ref_meta_parsed(input)
  sendsketch <- sendsketch_taxonomy_data_parsed(input, only_best = TRUE, only_shared = TRUE)

  # Find which columns are used to provide colors to the trees, if any
  ids_in_trees <- unique(unlist(lapply(trees, function(t) t$tip.label)))
  color_by_cols <- unique(unlist(strsplit(sample_meta$color_by[sample_meta$sample_id %in% ids_in_trees], split = ';')))
  color_by_cols <- color_by_cols[!is.na(color_by_cols)]
  color_by_col_names <- c(color_by_cols, 'Default')
  color_by_cols <- c(as.list(color_by_cols), list(NULL))  # NULL ensures that the default color scheme is also used

  # Plot one tree for each color_by column
  tree_plots <- lapply(color_by_cols, function(color_by) {
    plot_phylogeny(
      trees,
      sample_meta,
      ref_meta,
      color_by,
      sendsketch, # NOTE: A better source for the taxonomy should be found, perhaps based on reference NCBI taxon IDs.
      interactive = interactive
    )
  })
  names(tree_plots) <- color_by_col_names
  return(tree_plots)
}


#' Plot a phylogeny
#'
#' @param trees One or more [ape::phylo()] objects. Inputting multiple trees
#'   only makes sense when `collapse_by_tax` is used.
#' @param sample_meta A table containing metadata for samples in `tree_paths`.
#' @param ref_meta A table containing metadata for references in `tree_paths`.
#' @param collapse_by_tax A [tibble::tibble()] with the columns `sample_id` and
#'   those named by taxonomic ranks. These are ids present in `sample_meta`.
#'   These are used to provide a taxo`nomic tree that functions as a backbone to
#'   combine multiple trees into a single one. (Default: return a list of plots
#'   instead of a single plot)
#' @param color_by The name of the column in the metadata used to color samples.
#' @param interactive Whether to use an HTML-based interactive format or not
#'   (default: TRUE)
#'
#' @import ggtree
#'
#' @keywords internal
plot_phylogeny <- function(trees, sample_meta, ref_meta, color_by = NULL, collapse_by_tax = NULL, interactive = FALSE, ...) {
  if (interactive) {
    stop('Interactive tree plotting is not yet supported.')
  }

  # If a single tree is provided, convert it to a list of a single tree to simplify code below
  if (inherits(trees, "phylo")) {
    trees <- list(trees)
  }

  # Combine trees by connecting to the taxonomy-derived tree
  if (is.null(collapse_by_tax) || length(trees) == 1) {
    if (length(trees) == 1) {
      combined_tree <- trees[[1]]
    } else {
      stop(call. = FALSE, 'If "collapse_by_tax" is not provided, then only a single tree can be plotted at a time. Multiple trees were given.')
    }
  } else {
    # If taxonomy does not have at least one conserved rank, then add one
    ranks <- colnames(collapse_by_tax)[colnames(collapse_by_tax) != 'sample_id']
    n_unique_taxa <- unlist(lapply(collapse_by_tax[ranks], function(x) length(unique(x))))
    if (! any(n_unique_taxa == 1)) {
      collapse_by_tax <- cbind(sample_id = collapse_by_tax$sample_id, root = 'Life', collapse_by_tax[ranks])
      ranks <- c('root', ranks)
    }
    # Subset taxonomy data to just the part used by each tree
    collapse_by_tax[ranks] <- lapply(collapse_by_tax[ranks], as.factor)
    tree_tax_data <- lapply(trees, function(tree) {
      sample_ids <- tree$tip.label[tree$tip.label %in% sample_meta$sample_id]
      tax_subset <- collapse_by_tax[collapse_by_tax$sample_id %in% sample_ids, colnames(collapse_by_tax) != 'sample_id']
      tax_subset <- tax_subset[, apply(tax_subset, MARGIN = 2, function(col) length(unique(col)) == 1)]
      unique(tax_subset)
    })
    shared_cols <- table(unlist(lapply(tree_tax_data, colnames)))
    shared_cols <- names(shared_cols[as.numeric(shared_cols) == max(as.numeric(shared_cols))])
    tree_tax <- do.call(rbind, lapply(tree_tax_data, function(x) {
      x[colnames(x) %in% shared_cols]
    }))
    rownames(tree_tax) <- NULL
    # tree_tax <- tree_tax[, 1:ncol(tree_tax) < which(colnames(tree_tax) == 'g')]
    keep_rank <- apply(tree_tax, MARGIN = 2, function(col) length(unique(col)) != 1)
    keep_rank[1] <- TRUE # Always include the first rank
    tree_tax <- tree_tax[, keep_rank]
    tree_tax[] <- lapply(tree_tax, function(col) {
      col <- as.character(col)
      col[is.na(col)] <- 'test'
      as.factor(col)
    })
    # Combine trees
    base_tree <- ape::as.phylo(stats::as.formula(paste0('~', paste0(colnames(tree_tax), collapse = '/'))), data = tree_tax)
    mean_edge_len <- mean(unlist(lapply(trees, function(x) x$edge.length)))
    base_tree$edge.length <- rep(mean_edge_len, nrow(base_tree$edge))
    combined_tree <- base_tree
    index_key <- match(base_tree$tip.label, tree_tax[[ncol(tree_tax)]])
    for (index in rev(seq_along(trees))) {
      combined_tree <- ape::bind.tree(combined_tree, trees[[index_key[index]]], where = index)
    }
  }

  # Prepare labels
  label_key <- c(
    stats::setNames(sample_meta$description, sample_meta$sample_id),
    stats::setNames(ref_meta$ref_description, ref_meta$ref_id)
  )
  tip_labels <- label_key[combined_tree$tip.label]
  tip_labels <- stats::setNames(make.unique(tip_labels, sep = ' '), names(tip_labels))

  # Prepare colors
  if (is.null(color_by)) {
    color_by <- '_sequence_type_'
    sample_meta[['_sequence_type_']] <- 'Sample'
    ref_meta[['_sequence_type_']] <- 'Reference'
  }
  ids_in_trees <- combined_tree$tip.label
  color_key <- character(0)
  if (color_by %in% colnames(sample_meta)) {
    color_key <- c(color_key, stats::setNames(sample_meta[[color_by]], sample_meta$sample_id))
  }
  if (color_by %in% colnames(ref_meta)) {
    color_key <- c(color_key, stats::setNames(ref_meta[[color_by]], ref_meta$ref_id))
  }
  tip_colors <- color_key[combined_tree$tip.label]

  # Make data associated with tree tips
  tip_data <- tibble::tibble(
    newick_label = combined_tree$tip.label,
    tip_color = unname(tip_colors),
    tip_label = tip_labels[combined_tree$tip.label]
  )

  if (is.null(color_by) || length(trees) <= 1) {
    base_tree_node_labels <- character(0)
  } else {
    base_tree_node_labels <- base_tree$node.label
  }

  if (is.null(combined_tree$node.label)) {
    all_labels = c(rep(NA, combined_tree$Nnode), combined_tree$tip.label)
  } else {
    all_labels = c(combined_tree$node.label, combined_tree$tip.label)
  }

  node_data <- tibble::tibble(
    newick_label = all_labels,
    node_label = ifelse(all_labels %in% c(base_tree_node_labels, "Root"), "", newick_label),
    branch_color = ifelse(all_labels %in% c(base_tree_node_labels, "Root"), "grey", "black"),
    branch_type = ifelse(all_labels %in% c(base_tree_node_labels, "Root"), "dashed", "solid")
  )

  legend_title <- tools::toTitleCase(trimws(gsub(color_by, pattern = '_', replacement = ' ')))

  plotted_tree <- ggtree(combined_tree, ggplot2::aes_string(color = 'branch_color', linetype = 'branch_type')) %<+% node_data +
    geom_nodelab(ggplot2::aes_string(label = 'node_label'), hjust = 1.3, nudge_y = 0.3, size = 3) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_linetype_identity()

  plotted_tree <- plotted_tree %<+% tip_data +
    scale_x_continuous(expand = ggplot2::expansion(add = c(0.1, 0.1 + max(nchar(tip_data$tip_label)) * 0.1 / sqrt(nrow(tip_data))))) +
    ggnewscale::new_scale_color() +
    geom_tiplab(ggplot2::aes_string(label = 'tip_label', color = 'tip_color')) +
    geom_tippoint(ggplot2::aes_string(color = 'tip_color'), alpha = 0) + # Invisible tips just there to make override.aes below change the legend color shapes
    ggplot2::scale_color_viridis_d(end = 0.8, na.value = "black") +
    ggplot2::guides(color = guide_legend(title = legend_title, override.aes = list(label = "", size = 3, alpha = 1, shape = 15))) +
    theme(legend.position = "bottom")

  return(plotted_tree)
}


#' Plot MSN of variant data
#'
#' Plot a minimum spanning network for each group of samples aligned to a
#' reference found in `pathogensurveillance` output.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @param combine If `TRUE` combine multiple MSNs into a single figure and
#'   return a single figure. If `FALSE`, return a list of figures.
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'psminer')
#' variant_msn_plot(path)
#'
#' @export
variant_msn_plot <- function(path, combine = TRUE) {
  align_data <- variant_align_path_data(path)
  alignments <- variant_align_parsed(path)
  sample_data <- sample_meta_parsed(path)
  ref_data <- ref_meta_parsed(path)

  # Find which columns are used to provide colors to the trees, if any
  ids_used <- unique(unlist(lapply(alignments, function(a) {
    if (is.null(a)) {
      return(character(0))
    } else {
      return(rownames(a))
    }
  })))
  color_by_cols <- unique(unlist(strsplit(sample_data$color_by[sample_data$sample_id %in% ids_used], split = ';')))
  color_by_cols <- color_by_cols[! is.na(color_by_cols)]
  color_by_col_names <- c(color_by_cols, 'Default')
  color_by_cols <- c(as.list(color_by_cols), list(NULL))  # NULL ensures that the default color scheme is also used

  # Plot MSNs
  graphics::plot.new()
  output <- lapply(seq_len(nrow(align_data)), function(i) {
    align_without_ref <- alignments[[i]][rownames(alignments[[i]]) != align_data$ref_id[i], ]
    if (is.null(align_without_ref)) {
      return(NULL)
    }
    plot_data <- make_MSN(align_without_ref, sample_data, user_seed = 1, snp_diff_prop = 0.1, population = NULL)
    if (is.null(plot_data)) {
      return(NULL)
    }
    out <- grDevices::recordPlot()
    graphics::plot.new()
    return(out)
  })


}


#' Make Minimum spanning network
#'
#' @param snp_fasta_alignment A DNA alignment in fasta format.
#' @param sample_data A data frame containing information about samples.
#' @param population A character string specifying the column name in sample_data to be used for stratification.
#' @param interactive A logical value indicating whether Whether or not to produce an interactive HTML/javascript-based figures or tables or a static ones.
#' @param snp_threshold An integer specifying the number of SNPs to be used as the threshold for filtering. User can specify whole numbers or alternatively, a relative proportion using 'snp_diff_prop', but not both. Default is NULL.
#' @param snp_diff_prop A numeric value specifying the proportion of SNPs to be used as the threshold for filtering. If user prefers to specify whole numbers, use 'snp_threshold' instead. Default is NULL.
#' @param use_cutoff_predictor A logical value indicating whether to use cutoff predictor for determining SNP threshold.
#' @param show_MLG_table A logical value indicating whether to display a table of multi-locus genotypes.
#' @param user_seed An optional integer specifying the seed for reproducibility.
#' @param ... Additional arguments to be passed to internal functions.
#' @return Minimum spanning network
#'
#' @keywords internal
make_MSN <- function(snp_fasta_alignment, sample_data, population = NULL, interactive = FALSE, snp_threshold = NULL, snp_diff_prop = NULL, use_cutoff_predictor = FALSE, show_MLG_table = FALSE, user_seed = NULL, ...) {

  set.seed(user_seed)
  snp_aln.gi <- adegenet::DNAbin2genind(snp_fasta_alignment)
  if (is.null(snp_aln.gi)) {
    return(NULL)
  }
  snp_genclone <- poppr::as.genclone(snp_aln.gi)

  if (use_cutoff_predictor) {
    snpdist_stats <- poppr::filter_stats(snp_genclone)
    average_thresh <- poppr::cutoff_predictor(snpdist_stats$average$THRESHOLDS)
    poppr::mlg.filter(snp_genclone, distance = poppr::bitwise.dist, percent = TRUE) <- average_thresh
  } else if (!is.null(snp_threshold)) {
    poppr::mlg.filter(snp_genclone, distance = poppr::bitwise.dist, percent = FALSE) <- snp_threshold
  } else if (!is.null(snp_diff_prop)) {
    poppr::mlg.filter(snp_genclone, distance = poppr::bitwise.dist, percent = TRUE) <- snp_diff_prop
  }

  if (is.null(population)) {
    sample_data$no_factor_provided <- 'Sample (No factor provided)'
    population <- 'no_factor_provided'
  }

  sample_data <- sample_data[sample_data$sample_id %in% rownames(snp_fasta_alignment), , drop = FALSE]
  user_factor <- sample_data[[population]]
  node_color <- as.factor(ifelse(is.na(user_factor) | user_factor == "", "Unknown", user_factor))
  myColors <- grDevices::rainbow(length(unique(node_color)))
  names(myColors) <- levels(node_color)
  num_columns <- ncol(sample_data)
  adegenet::strata(snp_genclone) <- cbind(sample_data[, c(1:num_columns)], color_node_by = node_color)
  adegenet::setPop(snp_genclone) <- ~color_node_by

  if (!is.null(snp_threshold)) {
    ms.loc <- poppr::poppr.msn(snp_genclone,
                               distmat = poppr::bitwise.dist(snp_genclone, percent = FALSE),
                               include.ties = TRUE,
                               showplot = FALSE)
  }
  else {
    ms.loc <- poppr::poppr.msn(snp_genclone,
                               distmat = poppr::bitwise.dist(snp_genclone, percent = TRUE),
                               include.ties = TRUE,
                               showplot = FALSE)
  }

  the_edges <- igraph::E(ms.loc$graph)$weight
  edges <- as.list(the_edges)

  if (length(igraph::V(ms.loc$graph)) > 1) {
    output <- poppr::plot_poppr_msn(
      snp_genclone,
      poppr_msn = ms.loc,
      palette = myColors,
      mlg = FALSE,
      quantiles = FALSE,
      wscale = FALSE,
      inds = "None",
      ...
    )
  } else {
    output <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 4, y = 25, size=8, label = "All samples are in the same multilocus genotype.") +
      ggplot2::theme_void()
  }


  # if (show_MLG_table) {
  #   idlist <- mlg.id(snp_genclone)
  #   mlglist <- data.frame("MLG", "strain")
  #   colnames(mlglist) <- c("V1", "V2")
  #
  #   for (name in names(idlist)) {
  #     newframe <- as.data.frame(cbind(paste0("MLG", "_", name), idlist[[name]]))
  #     mlglist <- rbind(mlglist, newframe)
  #   }
  #
  #   colnames(mlglist) <- c("Multi-locus genotype", "Strain")
  #   mlglist <- mlglist[mlglist$Strain != "strain",]
  #
  #   if (interactive) {
  #     DT::datatable(mlglist, class = "display nowrap", ...) %>%
  #       formatStyle(colnames(mlglist), "white-space" = "nowrap")
  #
  #   } else {
  #     print(mlglist)
  #   }
  # }

  return(output)
}



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
#' @examples
#' make_ani_heatmap(ani_matrix, ref_data, samp_data, interactive=FALSE)
#'
#' @export
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


#' Make sunburst plot of sendsketch taxonomy
#'
#' Converts classifications of top hits in sendsketch output into an interactive
#' sunburst plot.
#'
#' @param input The path to one or more folders that contain
#'   pathogensurveillance output or a table in the format of the
#'   [sendsketch_parsed()] output.
#' @param interactive Whether or not to produce an interactive
#'   HTML/javascript-based plot or a static one.
#' @param ... Passed to `sendsketch_best_hits`
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'psminer')
#' sendsketch_taxonomy_plot(path)
#' sendsketch_taxonomy_plot(path, interactive = TRUE)
#'
#' @export
sendsketch_taxonomy_plot <- function(input, interactive = TRUE, ...) {

  # Parse the input if it is a file/folder path
  if (is.data.frame(input)) {
    sketch_data <- input
  } else {
    sketch_data <- sendsketch_parsed(input, only_best = TRUE)
  }

  # Sort and filter data
  top_hits <- sendsketch_best_hits(sketch_data, ...)

  # Make table with values duplicated for each taxon in all classifications
  split_taxa <- strsplit(top_hits$taxonomy, split = ';')
  tax_data <- do.call(rbind, lapply(seq_len(nrow(top_hits)), function(i) {
    n_taxa <- length(split_taxa[[i]])
    out <- top_hits[rep(i, n_taxa), , drop = FALSE]
    out$taxon <- vapply(seq_len(n_taxa), FUN.VALUE = character(1), function(n) {
      paste0(split_taxa[[i]][1:n], collapse = ';')
    })
    out$supertaxon <- c(NA_character_, out$taxon[1:(n_taxa - 1)])
    out$tip <- vapply(seq_len(n_taxa), FUN.VALUE = character(1), function(n) {
      split_taxa[[i]][n]
    })
    rownames(out) <- NULL
    return(out)
  }))

  # Convert to edge list
  plot_data <- unique(tax_data[, c('supertaxon', 'taxon', 'tip'), drop = FALSE])
  plot_data$count <- vapply(plot_data$taxon, FUN.VALUE = numeric(1), function(t) {
    sum(t == tax_data$taxon)
  })

  # Make unique IDs for each taxon
  id_key <- as.character(seq_len(nrow(plot_data)))
  names(id_key) <- plot_data$taxon
  plot_data$from <- id_key[plot_data$supertaxon]
  plot_data$to <- id_key[plot_data$taxon]

  # Rename taxon to just the tip taxon name and include rank if needed to make unique
  plot_data$rank <- sub(plot_data$tip, pattern = '^([a-z]+?):(.+)$', replacement = '\\1')
  plot_data$name <- sub(plot_data$tip, pattern = '^([a-z]+?):(.+)$', replacement = '\\2')
  plot_data$is_unique <- vapply(plot_data$name, FUN.VALUE = logical(1), function(x) {
    sum(x == plot_data$name) <= 1
  })
  plot_data$name <- ifelse(plot_data$is_unique, plot_data$name, plot_data$tip)

  output <- plotly::plot_ly(
    type = 'sunburst',
    ids = plot_data$to,
    labels = plot_data$name,
    parents = plot_data$from,
    values = plot_data$count,
    domain = list(column = 0),
    branchvalues = 'total'
  )

  if (! interactive) {
    temp_file_html <- tempfile(fileext = ".html")
    temp_file_png <- tempfile(fileext = ".png")
    htmlwidgets::saveWidget(widget = plotly::config(output, displayModeBar = FALSE), file = temp_file_html)
    output <- webshot2::webshot(url = temp_file_html, file = temp_file_png,
                                delay = 1, vheight = 750, vwidth = 750, zoom = 2)
  }

  return(output)
}
