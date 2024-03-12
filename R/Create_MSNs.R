make_MSN <- function(snp_fasta_alignment, sample_data, population=NULL, interactive = knitr::is_html_output(), snp_threshold=NULL, use_cutoff_predictor = TRUE, show_MLG_table=FALSE, user_seed=NULL, ...) {

  snp_aln.gi <- DNAbin2genind(snp_fasta_alignment)
  snp_aln.gi <- snp_aln.gi[indNames(snp_aln.gi) != "REF"]

  genind_names <- indNames(snp_aln.gi)
  cleaned_names <- sub(".*assembly_", "", genind_names)
  indNames(snp_aln.gi) <- cleaned_names

  mat <- match(indNames(snp_aln.gi), sample_data$sample_id)
  sample_data <- sample_data[mat, ]
  snp_genclone <- as.genclone(snp_aln.gi)

  # Define threshold options
  threshold_options <- c(0.0001, 0.001, 0.01, 0.1)

  if (use_cutoff_predictor) {
    snpdist_stats <- filter_stats(snp_genclone)
    average_thresh <- cutoff_predictor(snpdist_stats$average$THRESHOLDS)
    cat("Predicted SNP threshold, using cutoff_predictor function from poppr is:", average_thresh, "\n")
    mlg.filter(snp_genclone, distance = bitwise.dist, percent = TRUE) <- average_thresh
  } else if (!is.null(snp_threshold)) {
    mlg.filter(snp_genclone, distance = bitwise.dist, percent = TRUE) <- snp_threshold
    cat("User-defined SNP threshold is:", snp_threshold, "\n")
  } else {
    # Loop through threshold options and choose the first one
    for (thresh in threshold_options) {
      mlg.filter(snp_genclone, distance = bitwise.dist, percent = TRUE) <- thresh
      cat("Using threshold:", thresh, "\n")
      break
    }
  }

  # Create MSN based on population information
  if (!is.null(population) && population %in% names(sample_data)) {
    user_factor <- sample_data[[population]]
    node_color <- as.factor(ifelse(is.na(user_factor) | user_factor == "", "Unknown", user_factor))
    myColors <- rainbow(length(unique(node_color)))
    names(myColors) <- levels(node_color)
    num_columns <- ncol(sample_data)
    strata(snp_genclone) <- cbind(sample_data[, c(1:num_columns)], color_node_by = node_color)
    setPop(snp_genclone) <- ~color_node_by

    set.seed(user_seed)
    ms.loc <- poppr.msn(snp_genclone,
                        distmat = bitwise.dist(snp_genclone, percent = FALSE),
                        include.ties = TRUE,
                        showplot = FALSE)

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
      #layfun = igraph::layout_with_lgl,
      #edge.label = igraph::E(ms.loc$graph)$weight,
      #edge.label.font = 2,
      #edge.label.cex = 1,
      #edge.label.family = "Helvetica",
      #edge.label.color = "darkslateblue"
      ...
    )
  }
  # Create MSN based on color_by information
  else if (!is.null(sample_data$color_by)) {
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

      set.seed(user_seed)
      ms.loc <- poppr.msn(snp_genclone,
                          distmat = bitwise.dist(snp_genclone, percent = FALSE),
                          include.ties = TRUE,
                          showplot = FALSE)

      plot_poppr_msn(
        snp_genclone,
        poppr_msn = ms.loc,
        palette = myColors,
        mlg = FALSE,
        quantiles = FALSE,
        wscale = FALSE,
        inds = "None",
        #layfun = igraph::layout_with_lgl,
        #edge.label = igraph::E(ms.loc$graph)$weight,
        #edge.label.font = 2,
        #edge.label.cex = 1,
        #edge.label.family = "Helvetica",
        #edge.label.color = "darkslateblue"
      )
    }
  } else {
    # No color_by provided, create "No_Factor_Provided" label
    node_color <- as.factor(rep("No_Factor_Provided", length(indNames(snp_genclone))))
    myColors <- rainbow(length(unique(node_color)))
    names(myColors) <- levels(node_color)
    strata(snp_genclone) <- list(color_node_by = node_color)

    set.seed(user_seed)
    ms.loc <- poppr.msn(snp_genclone,
                        distmat = bitwise.dist(snp_genclone, percent = FALSE),
                        include.ties = TRUE,
                        showplot = FALSE)

    plot_poppr_msn(
      snp_genclone,
      poppr_msn = ms.loc,
      palette = myColors,
      mlg = FALSE,
      quantiles = FALSE,
      wscale = FALSE,
      inds = "None",
      #layfun = igraph::layout_with_lgl,
      #edge.label = igraph::E(ms.loc$graph)$weight,
      #edge.label.font = 2,
      #edge.label.cex = 1,
      #edge.label.family = "Helvetica",
      #edge.label.color = "darkslateblue"
    )
  }

  if (show_MLG_table) {
    idlist <- mlg.id(snp_genclone)
    mlglist <- data.frame("MLG","strain")
    colnames(mlglist) <- c("V1","V2")

    for (name in names(idlist)) {
      newframe <- as.data.frame(cbind(paste0("MLG","_",name),idlist[[name]]))
      mlglist <- rbind(mlglist,newframe)
    }


    colnames(mlglist) <- c("Multi-locus genotype","Strain")
    mlglist <- mlglist[mlglist$Strain != "strain",]

    # Print table
    if (interactive) {
      DT::datatable(mlglist, class = "display nowrap", ...) %>%
        formatStyle(colnames(mlglist), "white-space" = "nowrap")

    } else {
      print(mlglist)
    }
  }
}
