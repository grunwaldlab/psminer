#' Produce Summary and Detailed Tables for Pipeline Messages
#'
#' This function analyzes a dataframe of pipeline messages, producing three
#' types of tables: a simple summary table, a detailed table for samples, and a
#' detailed table for report groups. It supports both interactive and static
#' outputs, allowing users to choose the best format for their needs.
#' Interactive tables are generated using DT::datatable, providing a dynamic way
#' to explore the data. Static tables are produced with knitr::kable, suitable
#' for static reports and documentation.
#'
#'
#' @param input The path to one or more folders that contain
#'   pathogensurveillance output or a table in the format of the
#'   [status_message_parsed()] output.
#' @param interactive Logical; indicates whether to produce interactive tables
#'   (TRUE) or static tables (FALSE). Defaults to TRUE if the environment
#'   supports HTML output, otherwise FALSE. Interactive tables offer enhanced
#'   browsing capabilities, while static tables are best for printed pdf
#'   reports.
#'
#'
#' @return A list containing three tables: `simpleTable`,
#'   `detailedSamplesTable`, and `detailedGroupsTable`. `simpleTable` is a
#'   high-level summary of pipeline steps and associated issues.
#'   `detailedSamplesTable` and `detailedGroupsTable` provide in-depth insights
#'   into issues at the sample and group levels, respectively. Tables are
#'   presented as DT::datatable objects for interactive output or knitr::kable
#'   objects for static output.
#'
#' @importFrom dplyr filter group_by summarize left_join select mutate case_when
#' @importFrom DT datatable formatStyle
#' @importFrom magrittr %>%
#' @importFrom knitr kable is_html_output
#'
#' @export
status_message_tables <- function(input, interactive = knitr::is_html_output()) {
  # Parse the input if it is a file/folder path
  if (is.data.frame(input)) {
    messages <- input
  } else {
    messages <- status_message_parsed(input)
  }

  # Mapping workflow to major steps
  workflow_to_step <- list(
    genetic_diversity = c("ASSIGN_REFERENCES", "VARIANT_ANALYSIS"),
    initial_identification = "COARSE_SAMPLE_TAXONOMY",
    rigorous_identification = c("CORE_GENOME_PHYLOGENY", "DOWNLOAD_REFERENCES", "GENOME_ASSEMBLY"),
    report = c("CUSTOM_DUMPSOFTWAREVERSIONS", "MAIN_REPORT"),
    qc = c("FASTQC", "INPUT_CHECK", "MULTIQC", "RECORD_MESSAGES")
  )

  # Flatten the list to make it easier to use for mapping
  step_flat <- rep(names(workflow_to_step), vapply(workflow_to_step, FUN.VALUE = numeric(1), length))
  names(step_flat) <- unlist(workflow_to_step)

  # Map workflows to major steps
  messages$major_step <- step_flat[messages$workflow]

  # Adjusted summarization logic to use 'level' for calculating Failures and identifying group issues
  summary_data <- dplyr::group_by(messages, major_step) %>%
    dplyr::summarise(Failures = sum(level == "WARNING" | level == "ERROR"), .groups = 'drop')

  # Calculate report group issues
  report_group_issues <- dplyr::filter(messages, is.na(sample_id)) %>%
    dplyr::group_by(major_step) %>%
    dplyr::summarise(Group_Failures = n(), .groups = 'drop')
  report_group_issues$Group_Failures[is.na(report_group_issues$Group_Failures)] <- 0

  # Merge report group issues into summary data
  summary_data_full <- dplyr::left_join(summary_data, report_group_issues, by = "major_step") %>%
    dplyr::mutate(Sample_Issues = Failures - Group_Failures) %>%
    dplyr::select(major_step, Group_Failures, Sample_Issues)

  # Customize column names for the simple summary table
  colnames(summary_data_full) <- c("Pipeline step", "Group issues", "Sample issues")

  # Mapping message level to emojis for detailed tables
  messages$emoji <- dplyr::case_when(
    messages$level == "WARNING" ~ "⚠️ WARNING",
    messages$level == "ERROR" ~ "❌ ERROR",
    TRUE ~ "✅" # Default to a success symbol, adjust as necessary
  )

  # Detailed table for samples
  detailed_data_samples <- dplyr::filter(messages, !is.na(sample_id)) %>%
    dplyr::select(emoji, sample_id, major_step, summarized_message = message)

  # Customize column names for the detailed table for samples
  colnames(detailed_data_samples) <- c("Issue type", "Sample", "Pipeline step", "Message")

  # Detailed table for report groups
  detailed_data_groups <- messages %>%
    dplyr::filter(is.na(sample_id)) %>%
    dplyr::select(emoji, major_step, summarized_message = message)

  # Customize column names for the detailed table for report groups
  colnames(detailed_data_groups) <- c("Issue type", "Pipeline step", "Message")

  # Conditional output function and options
  if (interactive) {
    simple_table_output <- DT::datatable(summary_data_full, options = list(pageLength = 5, autoWidth = TRUE), escape = FALSE) %>%
      DT::formatStyle(colnames(summary_data_full), "white-space" = "nowrap")
    detailed_samples_output <- DT::datatable(detailed_data_samples, options = list(pageLength = 10, autoWidth = TRUE), escape = FALSE) %>%
      DT::formatStyle(colnames(detailed_data_samples), "white-space" = "nowrap")
    detailed_groups_output <- DT::datatable(detailed_data_groups, options = list(pageLength = 10, autoWidth = TRUE), escape = FALSE) %>%
      DT::formatStyle(colnames(detailed_data_groups), "white-space" = "nowrap")
  } else {
    # Using kable and kableExtra for static output
    simple_table_output <- print_static_table(summary_data_full)
    detailed_samples_output <- print_static_table(detailed_data_samples)
    detailed_groups_output <- print_static_table(detailed_data_groups)
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
