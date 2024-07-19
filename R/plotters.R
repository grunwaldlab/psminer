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
#' @return  A list of plots, unless `collapse_by_tax` is useded, in which case a single plot is returned.
#'
#' @export
core_tree_plot <- function(input, collapse_by_tax = NULL, interactive = knitr::is_html_output()) {
  # Find and parse needed data
  sample_meta <- sample_meta_parsed(input)
  ref_meta <- ref_meta_parsed(input)
  trees <-  core_tree_parsed(input)
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
#' @return  A list of plots, unless `collapse_by_tax` is useded, in which case a single plot is returned.
#'
#' @export
variant_tree_plot <- function(input, collapse_by_tax = NULL, interactive = knitr::is_html_output()) {
  # Find and parse needed data
  sample_meta <- sample_meta_parsed(input)
  ref_meta <- ref_meta_parsed(input)
  trees <-  variant_tree_parsed(input)
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
#' @keywords internal
plot_phylogeny <- function(trees, sample_meta, ref_meta, color_by = NULL, collapse_by_tax = NULL, interactive = knitr::is_html_output(), ...) {
  # If a single tree is provided, convert it to a list of a single tree to simplify code below
  if (class(trees) == "phylo") {
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
    # Subset taxonomy data to just the part used by each tree
    ranks <- colnames(collapse_by_tax)[colnames(collapse_by_tax) != 'sample_id']
    collapse_by_tax[ranks] <- lapply(collapse_by_tax[ranks], as.factor)
    tree_tax <- dplyr::bind_rows(lapply(trees, function(tree) {
      sample_ids <- tree$tip.label[tree$tip.label %in% sample_meta$sample_id]
      tax_subset <- collapse_by_tax[collapse_by_tax$sample_id %in% sample_ids, colnames(collapse_by_tax) != 'sample_id']
      tax_subset <- tax_subset[, apply(tax_subset, MARGIN = 2, function(col) length(unique(col)) == 1)]
      unique(tax_subset)
    }))
    tree_tax <- tree_tax[, 1:ncol(tree_tax) < which(colnames(tree_tax) == 'g')]
    keep_rank <- apply(tree_tax, MARGIN = 2, function(col) length(unique(col)) != 1)
    keep_rank[1] <- TRUE # Always include the first rank
    tree_tax <- tree_tax[, keep_rank]
    tree_tax[] <- lapply(tree_tax, function(col) {
      col <- as.character(col)
      col[is.na(col)] <- 'test'
      as.factor(col)
    })
    # Combine trees
    base_tree <- ape::as.phylo(as.formula(paste0('~', paste0(colnames(tree_tax), collapse = '/'))), data = tree_tax)
    mean_edge_len <- mean(unlist(lapply(trees, function(x) x$edge.length)))
    base_tree$edge.length <- rep(mean_edge_len, nrow(base_tree$edge))
    combined_tree <- base_tree
    for (index in rev(seq_along(trees))) {
      combined_tree <- ape::bind.tree(combined_tree, trees[[index]], where = index)
    }
  }

  # Prepare labels
  label_key <- c(
    stats::setNames(sample_meta$description, sample_meta$sample_id),
    stats::setNames(ref_meta$ref_description, ref_meta$ref_id)
  )
  tip_labels <- label_key[combined_tree$tip.label]
  tip_labels <- stats::setNames(make.unique(tip_labels, sep = ' '), names(tip_labels)) # phylocanvas does not plot anything if tip labels are not unique.

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

  # Plot tree
  if (interactive) {
    is_color <- function(x) {
      unlist(lapply(x, function(y) {
        tryCatch(is.matrix(col2rgb(y)), error = function(e) FALSE)
      }))
    }
    if (! all(is_color(tip_colors))) {
      factors <- unique(tip_colors)
      factors <- factors[! is.na(factors)]
      factor_key <- stats::setNames(viridis::viridis(length(factors), end = 0.8), factors)
      tip_colors <- stats::setNames(factor_key[tip_colors], names(tip_colors))
    }
    # Convert characters that are not allowed in the format used by phylocanvas
    tip_labels <- gsub(tip_labels, pattern = ',', replacement = '.', fixed = TRUE)
    tip_labels <- gsub(tip_labels, pattern = '[():; ]+', replacement = '_')
    tip_labels <- trimws(tip_labels, whitespace = '_')
    combined_tree$tip.label <- tip_labels
    names(tip_colors) <- tip_labels
    phycanv <- phylocanvas(combined_tree, treetype = "rectangular", alignlabels = TRUE, showscalebar = TRUE, width = "100%", height = "10in")
    # Add label colors
    phycanv$x$nodestyles <- lapply(tip_labels, function(x) {
      list(
        highlighted = FALSE,
        colour = "black",
        shape = "circle",
        size = 0.8,
        leafStyle = list(strokeStyle = unname(tip_colors[x]), fillStyle = unname(tip_colors[x]), lineWidth = 1),
        labelStyle = list(colour = unname(tip_colors[x]), textSize = 25, font = "Arial", format = "bold")
      )
    })
    names(phycanv$x$nodestyles) <- unname(tip_labels)
    # Remove underscores from tip labels (phylocanvas converts spaces to underscores)
    phycanv$x$tree <- gsub(phycanv$x$tree, pattern = '_', replacement = ' ', fixed = TRUE)
    names(phycanv$x$nodestyles) <- gsub(names(phycanv$x$nodestyles), pattern = '_', replacement = ' ', fixed = TRUE)
    return(phycanv)
  } else {

    tip_data <- tibble::tibble(
      newick_label = combined_tree$tip.label,
      tip_color = unname(tip_colors),
      tip_label = tip_labels[combined_tree$tip.label]
    )
    # tip_data$tip_color[is.na(tip_data$tip_color)] <- 'Undefined'

    if (is.null(color_by) || length(trees) <= 1) {
      base_tree_node_labels <- character(0)
    } else {
      base_tree_node_labels <- base_tree$node.label
    }

    node_data <- tibble::tibble(
      newick_label = c(combined_tree$node.label, combined_tree$tip.label),
      node_label = ifelse(newick_label %in% c(base_tree_node_labels, "Root"), "", newick_label),
      branch_color = ifelse(newick_label %in% c(base_tree_node_labels, "Root"), "grey", "black"),
      branch_type = ifelse(newick_label %in% c(base_tree_node_labels, "Root"), "dashed", "solid")
    )

    legend_title <- tools::toTitleCase(trimws(gsub(color_by, pattern = '_', replacement = ' ')))

    plotted_tree <- ggtree(combined_tree, aes(color = branch_color, linetype = branch_type)) %<+% node_data +
      geom_nodelab(aes(label = node_label), nudge_x = -0.04, nudge_y = 0.45) +
      scale_color_identity() +
      scale_linetype_identity()

    plotted_tree <- plotted_tree %<+% tip_data +
      scale_x_continuous(expand = expansion(mult = c(0.1, .2 + max(nchar(tip_data$tip_label)) * 0.03))) +
      ggnewscale::new_scale_color() +
      geom_tiplab(aes(label = tip_label, color = tip_color)) +
      geom_tippoint(aes(color = tip_color), alpha = 0) + # Invisible tips just there to make override.aes below change the legend color shapes
      scale_color_viridis_d(end = 0.8, na.value = "black") +
      guides(color = guide_legend(title = legend_title, override.aes = list(label = "", size = 3, alpha = 1, shape = 15))) +
      theme(legend.position = "bottom")

    return(plotted_tree)
  }
}

