#' Create an summary table from sendsketch output table
#'
#' This function processes sendsketch data from pathogen surveillance nextflow pipeline and
#' generates an interactive HTML table displaying the top hits from the sketch-based analysis.
#' The resulting table includes metrics such as Weighted Kmer IDentity (WKID),
#' Average Nucleotide Identity (ANI), and Completeness, providing insights into
#' genomic similarity and representational completeness between the query and reference genomes.
#' The table can be sorted by any of these metrics, and only the top results as specified
#' by `top_n` are displayed for each unique sample ID. Additionally, the table
#' visually encodes the percentage values using horizontal bars for an intuitive
#' and accessible presentation.
#'
#' @param sketch_data A dataframe containing sketch analysis results with the following columns:
#' @param sample_id Character; unique identifiers for each sample.
#' @param WKID Numeric; Weighted Kmer IDentity percentages indicating the kmer-based similarity.
#' @param ANI Numeric; Average Nucleotide Identity percentages representing genomic similarity.
#' @param Complt Numeric; Completeness percentages reflecting the coverage of the genome.
#' @param taxName Character; taxonomic names providing the preliminary identification of the sample.
#' @param sort_columns Character vector; specifies the columns to sort the table by in descending order.
#'   The default is `c("sample_id", "WKID", "ANI", "Complt")`. The `sample_id` is expected to be
#'   unique for each row, and the table will be sorted by this column first, followed by the other
#'   metrics in the order provided.
#' @param top_n Integer; the number of top entries to retain for each unique sample ID after sorting.
#'   The default is `1`, meaning only the top entry is kept.
#'
#' @return An interactive DT::datatable object that renders as an HTML table when used in an R Markdown document or Shiny application.
#' The table will have interactive features such as sorting and search enabled.
#'
#' @importFrom dplyr arrange group_by slice_head ungroup select rename
#' @importFrom magrittr %>%
#' @import DT
#' @export
#'
#' @examples
#' # Assuming `sketch_data` is your dataframe with the appropriate structure:
#' sketch_idtb(sketch_data, sort_columns = c("sample_id", "WKID", "ANI", "Complt"), top_n = 1)
#'
#' # If you want to sort by Completeness and then WKID, keeping the top 2 entries for each sample_id:
#' sketch_idtb(sketch_data, sort_columns = c("sample_id", "Complt", "WKID"), top_n = 2)
sketch_idtb <- function(sketch_data, sort_columns = c("sample_id", "WKID", "ANI", "Complt"), top_n = 1) {
  # Convert percentage fields from character to numeric
  sketch_data$WKID <- as.numeric(gsub("%", "", sketch_data$WKID))
  sketch_data$ANI <- as.numeric(gsub("%", "", sketch_data$ANI))
  sketch_data$Complt <- as.numeric(gsub("%", "", sketch_data$Complt))

  # Prepare sorting expressions
  sorting_exprs <- lapply(sort_columns, function(col) {
    expr <- paste0("desc(", col, ")")
    rlang::parse_expr(expr)
  })

  # Sort and filter data
  final_table <- sketch_data %>%
    arrange(!!!sorting_exprs) %>%
    group_by(sample_id) %>%
    slice_head(n = top_n) %>%
    ungroup() %>%
    select(sample_id, WKID, ANI, Complt, taxName) %>%
    rename(Sample = sample_id,
           `WKID (%)` = WKID,
           `ANI (%)` = ANI,
           `Completeness (%)` = Complt,
           `Top Hit` = taxName)

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
}
