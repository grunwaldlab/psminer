#' Parse one or more inputs to `pathogensurveillance`
#'
#' Parse and validate one or more input table of sample/reference metadata to
#' the `pathogensurveillance` pipeline and puts the parsed data in a single
#' table.
#'
#' @param input_paths The paths to one or more tables in TSV, CSV, or ODS
#'   (LibreOffice Calc), xls, or xlsx files (Excel) format.
#'
#' @examples
#'
#' input_paths <- system.file(package = 'psminer',
#'                            file.path('extdata', c(
#'                              'feature_test.tsv',
#'                              'feature_test_refs.csv'
#'                            )))
#'
#' @export
parse_path_surveil_input <- function(input_paths) {

}

#' Fix formatting errors
#'
#' Removes empty columns and rows, cleans up known column names, remove
#' duplicate columns, trims excess whitespace, etc
#'
#' @param table The input table
#' @param known_columns Column names to convert similar column names to
#'
#' @keywords internal
clean_input_table <- function(table, known_columns = known_input_columns()) {

  # Remove all whitespace
  table[] <- lapply(table, trimws)

  # Replace empty stings with NA
  table[] <- lapply(table, function(x) {
    x[x == ''] <- NA
    return(x)
  })

  # Remove empty rows
  is_empty <- apply(table, MARGIN = 1, function(row) all(is.na(row)))
  table <- table[! is_empty, , drop = FALSE]

  # Check that there is data
  if (nrow(table) == 0) {
    return(table)
  }

  # Fix names that appear to be known column names
  modified_names <- trimws(colnames(table))
  modified_names <- tolower(modified_names)
  modified_names <- gsub('[ -]+', '_', modified_names)
  underscore_replace_key <- known_columns
  names(underscore_replace_key) <- gsub('_+', '', known_columns)
  colnames_to_replace <- underscore_replace_key[modified_names]
  modified_names[!is.na(colnames_to_replace)] <- colnames_to_replace[! is.na(colnames_to_replace)]
  colnames(table)[modified_names %in% known_columns] <- modified_names[modified_names %in% known_columns]

  # Remove empty columns and columns with no header
  is_empty <- apply(table, MARGIN = 2, function(col) all(is.na(col)))
  is_headerless <- colnames(table) == '' | is.na(colnames(table))
  if (any(is_headerless & ! is_empty)) {
    stop(call. = FALSE, paste0(
      'The following columns in the ', csv_name, ' input CSV have values but no header: ',
      paste0(which(is_headerless & ! is_empty), collapse = ', ')
    ))
  }
  table <- table[, ! is_headerless, drop = FALSE]

  # Check for duplicated columns
  present_known_cols <- colnames(table)[colnames(table) %in% known_columns]
  duplicated_cols <- unique(present_known_cols[duplicated(present_known_cols)])
  if (length(duplicated_cols) > 0) {
    stop(call. = FALSE,
         'The following columns occur more than once in the ', csv_name, ' CSV: ',
         paste0('"', duplicated_cols, '"', collapse = ', ')
    )
  }

  # Reorder columns and add any missing columns
  empty_columns <- lapply(known_columns, function(column) {
    rep(NA, nrow(table))
  })
  names(empty_columns) <- known_columns
  reordered <- as.data.frame(empty_columns)
  reordered[colnames(table)] <- table
  table <- reordered

  # Add default values for some columns
  defaults <- default_column_values()
  table[names(defaults)] <- lapply(names(defaults), function(col) {
    ifelse(! is.na(table[[col]]), table[[col]], defaults[col])
  })

  return(table)
}

#' Read input tables from supported formats
#'
#' When a spreadsheet format such as .ods or .xlsx is supplied, each sheet is
#' returned as its own table.
#'
#' @param input_paths The paths to one or more files containing tables
#' @param read_all_sheets If `TRUE`, each sheet in a file is read as its own
#'   table. Otherwise, only the first sheet is read.
#'
#' @return A `list` of `tibble`.
#'
#' @param read_all_sheets If `TRUE`, read all sheets in spreadsheet files like
#'   .ods and .xlsx
#'
#' @examples
#' input_paths <- system.file(c('extdata/xls_example.xls', 'extdata/feature_test.tsv'), package = 'psminer')
#' psminer:::read_input_tables(input_paths)
#'
#' @keywords internal
read_input_tables <- function(input_paths, read_all_sheets = TRUE) {
  if (length(input_paths) == 0) {
    return(list())
  }
  input_paths <- unique(input_paths)
  read_one <- function(path) {
    if (endsWith(path, '.csv')) {
      output <- list(tibble::as_tibble(read.csv(path, check.names = FALSE)))
      names(output) <- path
    } else if (endsWith(path, '.tsv')) {
      output <- list(tibble::as_tibble(read.csv(path, check.names = FALSE, sep = '\t')))
      names(output) <- path
    } else if (endsWith(path, '.ods')) {
      sheets <- readODS::list_ods_sheets(path)
      output <- lapply(sheets, function(sheet) {
        readODS::read_ods(path, sheet = sheet)
      })
      names(output) <- paste0(path, '::', sheets)
    } else if (endsWith(path, '.xls') || endsWith(path, '.xlsx')) {
      sheets <- readxl::excel_sheets(path)
      output <- lapply(sheets, function(sheet) {
        readxl::read_excel(path, sheet = sheet)
      })
      names(output) <- paste0(path, '::', sheets)
    } else {
      stop('Input file extension not supported. Must be .tsv, .csv, .ods, .xls, or .xlsx')
    }
    return(output)
  }
  unlist(lapply(input_paths, read_one), recursive = FALSE)
}


#' Columns used by `pathogensurveillance`
#'
#' Column in input tables that are used by the `pathogensurveillance` pipeline.
#' Other columns can be given, but only these are used.
#'
#' @return A `vector` of column names
#'
#' @examples
#' known_input_columns()
#'
#' @export
known_input_columns <- function() {
  c(
    'id',
    'ref_group_id',
    'report_id',
    'name',
    'description',
    'input_type',
    'data_type',
    'data_source',
    'enabled',
    'ref_primary_usage',
    'ref_contextual_usage',
    'color_by',
    'ploidy',
    'latitude',
    'longitude',
    'location',
    'country',
    'region',
    'subregion',
    'place',
    'district',
    'date',
    'year',
    'month',
    'day',
    'hour',
    'minute',
    'second',
    'link'
  )
}


#' Defaults for known columns
#'
#' These are used if a cell is left blank.
#'
#' @return A `vector` of default values named by column names.
#'
#' @examples
#' default_column_values()
#'
#' @keywords internal
default_column_values <- function() {
  c(
    input_type = 'sample',
    report_group_id = 'miscellaneous',
    enabled = TRUE,
    query_max = 10,
    ref_primary_usage = 'optional',
    ref_contextual_usage = 'optional',
    color_by = 'input_type'
  )
}


#' Types of data supported by the pipeline
#'
#' These are the valid values that can be supplied to the `data_type` column.
#'
#' @return A character `vector`
#'
#' @examples
#' known_data_types()
#'
#' @export
known_data_types <- function() {
  c(
    'illumina',
    'nanopore',
    'pacbio',
    'assembly',
    'ncbi accession',
    'ncbi sra query',
    'ncbi assembly query',
    'image',
    'observation'
  )
}

#' Characters that cant be in IDs
#'
#' Regular expression for characters that cannot appear in IDs
#'
#' @return A character `vector` of length 1
#'
#' @examples
#' invalid_id_char_pattern()
#'
#' @export
invalid_id_char_pattern <- function() {
  '[\\/:;*?"<>| .()-]+'
}

