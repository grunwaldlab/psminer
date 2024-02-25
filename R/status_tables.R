#' Produce Summary and Detailed Tables for Pipeline Messages
#'
#' This function takes a dataframe containing pipeline messages and produces three tables:
#' a simple summary table, a detailed table for samples, and a detailed table for report groups.
#' It allows toggling between interactive and static outputs.
#'
#' @param messages A dataframe containing the columns `sample_id`, `group_id`, `workflow`, `message_level`, and `message_text`.
#' @param interactive Logical; whether to produce interactive tables (TRUE) or static tables (FALSE).
#'                    Default is TRUE if running in an environment that supports HTML output, otherwise FALSE.
#' @return A list containing three objects: simpleTable, detailedSamplesTable, and detailedGroupsTable.
#'         Each object is either a DT::datatable (for interactive output) or a knitr::kable (for static output).
#' @export
#' @examples
#' messages <- readr::read_tsv("path/to/messages.tsv") %>%
#'   tidyr::separate(message, into = c("message_level", "message_text"), sep = "\t")
#' tables <- status_tables(messages, interactive = TRUE)
#' tables$simpleTable
#' tables$detailedSamplesTable
#' tables$detailedGroupsTable
status_tables <- function(messages, interactive = knitr::is_html_output()) {
  library(dplyr)
  library(tidyr)
  library(DT)
  library(purrr)
  library(knitr)

  # Mapping workflow to major steps
  workflow_to_step <- list(
    genetic_diversity = c("ASSIGN_REFERENCES", "VARIANT_ANALYSIS"),
    initial_identification = "COARSE_SAMPLE_TAXONOMY",
    rigorous_identification = c("CORE_GENOME_PHYLOGENY", "DOWNLOAD_REFERENCES", "GENOME_ASSEMBLY"),
    report = c("CUSTOM_DUMPSOFTWAREVERSIONS", "MAIN_REPORT"),
    qc = c("FASTQC", "INPUT_CHECK", "MULTIQC", "RECORD_MESSAGES")
  )

  # Flatten the list to make it easier to use for mapping
  step_flat <- purrr::imap(workflow_to_step, function(val, name) setNames(rep(name, length(val)), val)) %>%
    purrr::reduce(c)

  # Map workflows to major steps
  messages$major_step <- step_flat[messages$workflow]

  # Adjusted summarization logic to use 'message_level' for calculating Failures and identifying group issues
  summary_data <- messages %>%
    group_by(major_step) %>%
    summarise(Failures = sum(message_level == "WARNING" | message_level == "ERROR"), .groups = 'drop')

  # Calculate report group issues
  report_group_issues <- messages %>%
    filter(is.na(sample_id)) %>%
    group_by(major_step) %>%
    summarise(Group_Failures = n(), .groups = 'drop')

  # Merge report group issues into summary data
  summary_data_full <- summary_data %>%
    left_join(report_group_issues, by = "major_step") %>%
    replace_na(list(Group_Failures = 0)) %>%
    mutate(Sample_Issues = Failures - Group_Failures) %>%
    select(major_step, Group_Failures, Sample_Issues)

  # Customize column names for the simple summary table
  colnames(summary_data_full) <- c("Pipeline step", "Group issues", "Sample issues")

  # Mapping message level to emojis for detailed tables
  messages$emoji <- case_when(
    messages$message_level == "WARNING" ~ "⚠️WARNING",
    messages$message_level == "ERROR" ~ "❌ERROR",
    TRUE ~ "✅" # Default to a success symbol, adjust as necessary
  )

  # Detailed table for samples
  detailed_data_samples <- messages %>%
    filter(!is.na(sample_id)) %>%
    select(sample_id, major_step, emoji, summarized_message = message_text)

  # Customize column names for the detailed table for samples
  colnames(detailed_data_samples) <- c("Sample", "Pipeline step", "Issue type", "Message")

  # Detailed table for report groups
  detailed_data_groups <- messages %>%
    filter(is.na(sample_id)) %>%
    select(group_id, major_step, emoji, summarized_message = message_text)

  # Customize column names for the detailed table for report groups
  colnames(detailed_data_groups) <- c("Group", "Pipeline step", "Issue type", "Message")

  # Conditional output function and options
  if (interactive) {
    simple_table_output <- DT::datatable(summary_data_full, options = list(pageLength = 5, autoWidth = TRUE), escape = FALSE)
    detailed_samples_output <- DT::datatable(detailed_data_samples, options = list(pageLength = 10, autoWidth = TRUE), escape = FALSE)
    detailed_groups_output <- DT::datatable(detailed_data_groups, options = list(pageLength = 10, autoWidth = TRUE), escape = FALSE)
  } else {
    # Using kable and kableExtra for static output
    simple_table_output <- summary_data_full %>%
      kable("html", escape = FALSE) %>%
      kable_styling(bootstrap_options = c("striped", "hover"))
    detailed_samples_output <- detailed_data_samples %>%
      kable("html", escape = FALSE) %>%
      kable_styling(bootstrap_options = c("striped", "hover"))
    detailed_groups_output <- detailed_data_groups %>%
      kable("html", escape = FALSE) %>%
      kable_styling(bootstrap_options = c("striped", "hover"))
  }

  # Return the tables as a list
  return(
    list(
      simpleTable = simple_table_output,
      detailedSamplesTable = detailed_samples_output,
      detailedGroupsTable = detailed_groups_output
    )
  )
}
