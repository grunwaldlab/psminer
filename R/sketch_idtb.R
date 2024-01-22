#' Create an Identification Table from Sketch Data
#'
#' This function processes sketch data from pathogen surveillance reports and
#' generates an interactive table displaying the top hits from the sketch-based analysis.
#' It includes metrics like Weighted Kmer IDentity (WKID), Average Nucleotide Identity (ANI),
#' and Completeness, providing insights into genomic similarity and representational completeness.
#'
#' @param sketch_data A dataframe containing sketch analysis results, with columns
#' for sample IDs, WKID, ANI, Completeness, and taxonomic names.
#' @param sort_columns A vector of column names to sort by, in descending order.
#' @param top_n The number of top entries to retain for each sample after sorting.
#' @return An interactive table rendered using the DT package, suitable for HTML output.
#' @import dplyr
#' @import DT
#' @export
#' @examples
#' # Assuming `sketch_data` is your dataframe
#' create_identification_table(sketch_data, sort_columns = c("sample_id", "WKID", "ANI", "Complt"), top_n = 1)
sketch_idtb <- function(sketch_data, sort_columns = c("sample_id", "WKID", "ANI", "Complt"), top_n = 1) {
  # Convert percentage fields from character to numeric
  sketch_data$WKID <- as.numeric(gsub("%", "", sketch_data$WKID))
  sketch_data$ANI <- as.numeric(gsub("%", "", sketch_data$ANI))
  sketch_data$Complt <- as.numeric(gsub("%", "", sketch_data$Complt))

  # Sort and filter data
  final_table <- sketch_data %>%
    arrange(desc(sample_id), desc(WKID), desc(ANI), desc(Complt)) %>%  # Sort by sample_id, WKID, ANI, and Complt in descending order
    group_by(sample_id) %>%  # Group by sample_id
    slice_head(n = 1) %>%  # Take the first entry per group
    ungroup() %>%  # Remove grouping
    select(sample_id, WKID, ANI, Complt, taxName) %>%  # Select required columns
    rename(Sample = sample_id,
           `WKID (%)` = WKID,
           `ANI (%)` = ANI,
           `Completeness (%)` = Complt,
           `Top Hit` = taxName)  # Rename columns

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
