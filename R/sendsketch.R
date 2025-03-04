#' Create an summary table from sendsketch output table
#'
#' This function processes sendsketch data from pathogen surveillance nextflow
#' pipeline and generates an interactive HTML table displaying the top hits from
#' the sketch-based analysis. The resulting table includes metrics such as
#' Weighted Kmer IDentity (WKID), Average Nucleotide Identity (ANI), and
#' Completeness, providing insights into genomic similarity and representational
#' completeness between the query and reference genomes. The table visually
#' encodes the percentage values using horizontal bars for an intuitive and
#' accessible presentation.
#'
#' @param input The path to one or more folders that contain
#'   pathogensurveillance output or a table in the format of the
#'   [sendsketch_parsed()] output. Only the best hits are returned, based on the
#'   default behavior of [sendsketch_best_hits()]. To change this behavior, pass
#'   in the results of running [sendsketch_best_hits()].
#' @param interactive Whether or not to produce an interactive
#'   HTML/javascript-based table or a static one.
#'
#' @return An interactive DT::datatable object that renders as an HTML table
#'   when used in an R Markdown document or Shiny application. The table will
#'   have interactive features such as sorting and search enabled.
#'
#' @export
#'
#' @examples
#' # Assuming `sketch_data` is your dataframe with the appropriate structure:
#' sketch_idtb(sketch_data, sort_columns = c("sample_id", "WKID", "ANI", "Complt"), top_n = 1)
#'
#' # If you want to sort by Completeness and then WKID, keeping the top 2 entries for each sample_id:
#' sketch_idtb(sketch_data, sort_columns = c("sample_id", "Complt", "WKID"), top_n = 2)
sendsketch_table <- function(input, interactive = knitr::is_html_output()) {
  # Parse the input if it is a file/folder path
  if (is.data.frame(input)) {
    sketch_data <- input
  } else {
    sketch_data <- sendsketch_parsed(input, only_best = TRUE)
  }

  # Sort and filter data
  final_table <- sketch_data[, c('sample_id', 'WKID', 'ANI', 'Complt', 'taxName'), drop = FALSE]
  names(final_table) <- c('Sample', 'WKID (%)', 'WKID (%)', 'c', 'Completeness (%)', 'Top Hit')

  if (interactive) {
    # Define a function called 'bordered_bar'
    bordered_bar <- function(value, color) {
      sprintf('<div style="border: 1px solid gray; width: 100%%; border-radius: 12px;">
              <div style="width: %s%%; background-color: %s; border-radius: 12px; text-align: center;">%s</div>
             </div>',
              value, color, value)
    }

    # Apply the bordered_bar function
    final_table$`WKID (%)` <- sapply(final_table$`WKID (%)`, function(x) bordered_bar(x, 'lightblue'))
    final_table$`ANI (%)` <- sapply(final_table$`ANI (%)`, function(x) bordered_bar(x, 'lightgreen'))
    final_table$`Completeness (%)` <- sapply(final_table$`Completeness (%)`, function(x) bordered_bar(x, 'lightpink'))

    # Render the table using DT for HTML output
    return(
      DT::datatable(final_table,
                    options = list(
                      pageLength = 10,
                      autoWidth = TRUE,
                      columnDefs = list(
                        list(width = '150px', targets = c(1, 2, 3))
                      ),
                      searchHighlight = TRUE
                    ),
                    rownames = FALSE,
                    escape = FALSE
      )
    )
  } else {
    return(print_static_table(final_table))
  }
}



#' Pick best sendsketch hits for each sample
#'
#' This function processes sendsketch data from pathogen surveillance nextflow pipeline and
#' filters the output so that only the best hits are present
#'
#' @param sketch_data A dataframe containing sketch analysis results.
#' @param sort_columns Character vector; specifies the columns to sort the table by in descending order.
#'   The default is `c("sample_id", "WKID", "ANI", "Complt")`. The `sample_id` is expected to be
#'   unique for each row, and the table will be sorted by this column first, followed by the other
#'   metrics in the order provided.
#' @param top_n Integer; the number of top entries to retain for each unique sample ID after sorting.
#'   The default is `1`, meaning only the top entry is kept.
#'
#' @return A `data.frame`
#'
#' @export
#'
#' @examples
#' # Assuming `sketch_data` is your dataframe with the appropriate structure:
#' sendsketch_best_hits(sketch_data, sort_columns = c("sample_id", "WKID", "ANI", "Complt"), top_n = 1)
#'
#' # If you want to sort by Completeness and then WKID, keeping the top 2 entries for each sample_id:
#' sendsketch_best_hits(sketch_data, sort_columns = c("sample_id", "Complt", "WKID"), top_n = 2)
sendsketch_best_hits <- function(sketch_data, sort_columns = c("WKID", "ANI", "Complt"), top_n = 1) {
  sketch_data <- sketch_data[do.call(order, sketch_data[sort_columns]), , drop = FALSE]
  split_data <- split(sketch_data, sketch_data[c('sample_id', 'report_group_id')])
  final_table <- do.call(rbind, lapply(split_data, function(x) {
    x[seq_len(top_n)]
  }))
  return(final_table)
}



#' Make sunburst plot of sendsketch taxonomy
#'
#' Converts classifications of top hits in sendsketch output into an interactive sunburst plot.
#'
#' @param input The path to one or more folders that contain
#'   pathogensurveillance output or a table in the format of the
#'   [sendsketch_parsed()] output.
#' @param interactive Whether or not to produce an interactive HTML/javascript-based plot or a static one.
#'
#' @export
sendsketch_taxonomy_plot <- function(input, interactive = knitr::is_html_output(), ...) {

  # Parse the input if it is a file/folder path
  if (is.data.frame(input)) {
    sketch_data <- input
  } else {
    sketch_data <- sendsketch_parsed(input, only_best = TRUE)
  }

  # Sort and filter data
  top_hits <- sendsketch_best_hits(sketch_data, ...)

  # Parse taxonomic data
  x <- metacoder::parse_tax_data(tax_data = top_hits,
                                 class_cols = 'taxonomy',
                                 class_key = c("taxon_rank", "taxon_name"),
                                 class_regex = "([a-z]+)?:?([a-zA-Z0-9.-_, ]+)",
                                 class_sep = ";")
  x <- metacoder::filter_taxa(x, taxon_ranks == "s", supertaxa = TRUE)

  # Replace duplicated names with their name + rank
  duplicated_names <- metacoder::taxon_names(x)[duplicated(metacoder::taxon_names(x))]
  unique_taxon_names <- paste0(metacoder::taxon_names(x), " (", metacoder::taxon_ranks(x), ")")
  names(unique_taxon_names) <- metacoder::taxon_ids(x)
  unique_tax_names <- ifelse(metacoder::taxon_names(x) %in% duplicated_names, unique_taxon_names, metacoder::taxon_names(x))
  names(unique_tax_names) <- metacoder::taxon_ids(x)

  # Convert to an edge list for plotting
  plot_data <- x$edge_list
  plot_data$count <- metacoder::n_obs(x)[plot_data$to]
  plot_data$from <- unique_tax_names[plot_data$from]
  plot_data$to <- unique_tax_names[plot_data$to]

  output <- plotly::plot_ly(
    type = 'sunburst',
    ids = plot_data$to,
    labels = plot_data$to,
    parents = plot_data$from,
    values = plot_data$count,
    domain = list(column = 0),
    branchvalues = 'total'
  )

  if (! interactive) {
    temp_file_html <- tempfile(fileext = ".html")
    temp_file_png <- tempfile(fileext = ".png")
    htmlwidgets::saveWidget(widget = config(output, displayModeBar = FALSE), file = temp_file_html)
    output <- webshot2::webshot(url = temp_file_html, file = temp_file_png,
                                delay = 1, vheight = 750, vwidth = 750, zoom = 2)
  }

  return(output)
}
