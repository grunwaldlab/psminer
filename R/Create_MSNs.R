#' Make Minimum spanning network
#'
#' @param snp_fasta_alignment A DNA alignment in fasta format.
#' @param sample_data A data frame containing information about samples.
#' @param population A character string specifying the column name in sample_data to be used for stratification.
#' @param interactive A logical value indicating whether Whether or not to produce an interactive HTML/javascript-based figures or tables or a static ones (default is determined by knitr::is_html_output()).
#' @param snp_threshold An integer specifying the number of SNPs to be used as the threshold for filtering. User can specify whole numbers or alternatively, a relative proportion using 'snp_diff_prop', but not both. Default is NULL.
#' @param snp_diff_prop A numeric value specifying the proportion of SNPs to be used as the threshold for filtering. If user prefers to specify whole numbers, use 'snp_threshold' instead. Default is NULL.
#' @param use_cutoff_predictor A logical value indicating whether to use cutoff predictor for determining SNP threshold.
#' @param show_MLG_table A logical value indicating whether to display a table of multi-locus genotypes.
#' @param user_seed An optional integer specifying the seed for reproducibility.
#' @param ... Additional arguments to be passed to internal functions.
#' @return Minimum spanning network
#' @export
make_MSN <- function(snp_fasta_alignment, sample_data, population = NULL, interactive = knitr::is_html_output(), snp_threshold = NULL, snp_diff_prop = NULL, use_cutoff_predictor = FALSE, show_MLG_table = FALSE, user_seed = NULL, ...) {

  set.seed(user_seed)
  snp_aln.gi <- DNAbin2genind(snp_fasta_alignment)
  snp_aln.gi <- snp_aln.gi[indNames(snp_aln.gi) != "REF"]

  genind_names <- indNames(snp_aln.gi)
  cleaned_names <- sub(".*assembly_", "", genind_names)
  indNames(snp_aln.gi) <- cleaned_names

  mat <- match(indNames(snp_aln.gi), sample_data$sample_id)
  sample_data <- sample_data[mat, ]
  snp_genclone <- as.genclone(snp_aln.gi)


  if (use_cutoff_predictor) {
    snpdist_stats <- filter_stats(snp_genclone)
    average_thresh <- cutoff_predictor(snpdist_stats$average$THRESHOLDS)
    mlg.filter(snp_genclone, distance = bitwise.dist, percent = TRUE) <- average_thresh
  } else if (!is.null(snp_threshold)) {
    mlg.filter(snp_genclone, distance = bitwise.dist, percent = FALSE) <- snp_threshold
  } else if (!is.null(snp_diff_prop)) {
    mlg.filter(snp_genclone, distance = bitwise.dist, percent = TRUE) <- snp_diff_prop
  }

  if (!is.null(population) && population %in% names(sample_data)) {
    user_factor <- sample_data[[population]]
    node_color <- as.factor(ifelse(is.na(user_factor) | user_factor == "", "Unknown", user_factor))
    myColors <- rainbow(length(unique(node_color)))
    names(myColors) <- levels(node_color)
    num_columns <- ncol(sample_data)
    strata(snp_genclone) <- cbind(sample_data[, c(1:num_columns)], color_node_by = node_color)
    setPop(snp_genclone) <- ~color_node_by

    if (!is.null(snp_threshold)) {
      ms.loc <- poppr.msn(snp_genclone,
                          distmat = bitwise.dist(snp_genclone, percent = FALSE),
                          include.ties = TRUE,
                          showplot = FALSE)
    }
    else {
      ms.loc <- poppr.msn(snp_genclone,
                          distmat = bitwise.dist(snp_genclone, percent = TRUE),
                          include.ties = TRUE,
                          showplot = FALSE)
    }

    the_edges <- igraph::E(ms.loc$graph)$weight
    edges <- as.list(the_edges)

    plot_poppr_msn(
      snp_genclone,
      poppr_msn = ms.loc,
      palette = myColors,
      mlg = FALSE,
      quantiles = FALSE,
      wscale = FALSE,
      inds = "None",
      ...
    )
  } else if (!is.null(sample_data$color_by)) {
    unique_factors <- unique(sample_data$color_by)
    unique_factors <- unlist(strsplit(unique_factors, ";"))

    for (factor in unique_factors) {
      factor_column <- sample_data[[factor]]
      node_color <- as.factor(ifelse(is.na(factor_column) | factor_column == "", "Unknown", factor_column))
      myColors <- rainbow(length(unique(node_color)))
      names(myColors) <- levels(node_color)
      num_columns <- ncol(sample_data)
      strata(snp_genclone) <- cbind(sample_data[, c(1:num_columns)], color_node_by = node_color)
      setPop(snp_genclone) <- ~color_node_by

      if (!is.null(snp_threshold)) {
        ms.loc <- poppr.msn(snp_genclone,
                            distmat = bitwise.dist(snp_genclone, percent = FALSE),
                            include.ties = TRUE,
                            showplot = FALSE)
      }
      else {
        ms.loc <- poppr.msn(snp_genclone,
                            distmat = bitwise.dist(snp_genclone, percent = TRUE),
                            include.ties = TRUE,
                            showplot = FALSE)
      }

      plot_poppr_msn(
        snp_genclone,
        poppr_msn = ms.loc,
        palette = myColors,
        mlg = FALSE,
        quantiles = FALSE,
        wscale = FALSE,
        inds = "None",
        ...
      )
    }
  } else {
    node_color <- as.factor(rep("No_Factor_Provided", length(indNames(snp_genclone))))
    myColors <- rainbow(length(unique(node_color)))
    names(myColors) <- levels(node_color)
    strata(snp_genclone) <- list(color_node_by = node_color)

    if (!is.null(snp_threshold)) {
      ms.loc <- poppr.msn(snp_genclone,
                          distmat = bitwise.dist(snp_genclone, percent = FALSE),
                          include.ties = TRUE,
                          showplot = FALSE)
    }
    else {
      ms.loc <- poppr.msn(snp_genclone,
                          distmat = bitwise.dist(snp_genclone, percent = TRUE),
                          include.ties = TRUE,
                          showplot = FALSE)
    }

    plot_poppr_msn(
      snp_genclone,
      poppr_msn = ms.loc,
      palette = myColors,
      mlg = FALSE,
      quantiles = FALSE,
      wscale = FALSE,
      inds = "None",
      ...
    )
  }

  if (show_MLG_table) {
    idlist <- mlg.id(snp_genclone)
    mlglist <- data.frame("MLG", "strain")
    colnames(mlglist) <- c("V1", "V2")

    for (name in names(idlist)) {
      newframe <- as.data.frame(cbind(paste0("MLG", "_", name), idlist[[name]]))
      mlglist <- rbind(mlglist, newframe)
    }

    colnames(mlglist) <- c("Multi-locus genotype", "Strain")
    mlglist <- mlglist[mlglist$Strain != "strain",]

    if (interactive) {
      DT::datatable(mlglist, class = "display nowrap", ...) %>%
        formatStyle(colnames(mlglist), "white-space" = "nowrap")

    } else {
      print(mlglist)
    }
  }
}
