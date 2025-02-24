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
#'                              'feature_test_refs.csv',
#'                              'ods_example.ods'
#'                            )))
#' parse_path_surveil_input(input_paths)
#'
#' @export
parse_path_surveil_input <- function(input_paths) {

  # Read and standardize all inputs into a single table without any data interpretation
  metadata <- read_input_tables(input_paths)



}


#' Read input tables from supported formats
#'
#' When a spreadsheet format such as .ods or .xlsx is supplied, each sheet is
#' returned as its own table.
#'
#' @param input_paths The paths to one or more files containing tables
#' @param read_all_sheets If `TRUE`, each sheet in a file is read as its own
#'   table. Otherwise, only the first sheet is read.
#' @param add_input_metadata Whether to add `_input_path_`, `_sheet_name_` and
#'   `_row_index_` columns to each table to record where each row came from.
#'
#' @return A `list` of `tibble`.
#'
#' @param read_all_sheets If `TRUE`, read all sheets in spreadsheet files like
#'   .ods and .xlsx
#'
#' @examples
#' input_paths <- system.file(package = 'psminer',
#'                            file.path('extdata', c(
#'                              'feature_test.tsv',
#'                              'feature_test_refs.csv',
#'                              'ods_example.ods'
#'                            )))
#' psminer:::read_input_tables(input_paths)
#'
#' @keywords internal
read_input_tables <- function(input_paths, read_all_sheets = TRUE, add_input_metadata = TRUE) {

  # Read data from whatever format it is stored in
  input_paths <- unique(input_paths)
  read_one <- function(path) {
    if (endsWith(path, '.csv')) {
      output <- list(tibble::as_tibble(read.csv(path, check.names = FALSE)))
      names(output) <- path
      sheets <- NA
    } else if (endsWith(path, '.tsv')) {
      output <- list(tibble::as_tibble(read.csv(path, check.names = FALSE, sep = '\t')))
      names(output) <- path
      sheets <- NA
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
    if (add_input_metadata) {
      for (i in seq_along(output)) {
        output[[i]]['_input_path_'] <- path
        output[[i]]['_sheet_name_'] <- sheets[i]
        output[[i]]['_row_index_'] <- seq_len(nrow(output[[i]]))
      }
    }
    return(lapply(output, clean_input_table))
  }
  tables <- unlist(lapply(input_paths, read_one), recursive = FALSE)

  # Reorder columns, add any missing columns, and combine
  user_columns <- unique(unlist(lapply(tables, colnames)))
  all_columns <- unique(c(known_input_columns(), user_columns))
  tables <- lapply(tables, function(table) {
    empty_columns <- lapply(all_columns, function(column) {
      rep(NA, nrow(table))
    })
    names(empty_columns) <- all_columns
    reordered <- tibble::as_tibble(empty_columns)
    reordered[colnames(table)] <- table
    return(reordered)
  })
  combined_tables <- do.call(rbind, tables)

  # Add default values for some columns
  defaults <- default_column_values()
  combined_tables[names(defaults)] <- lapply(names(defaults), function(col) {
    ifelse(is.na(combined_tables[[col]]), rep(defaults[col], nrow(combined_tables)), combined_tables[[col]])
  })

  # Clean up columns with multiple inputs
  multi_cols <- multi_input_columns()
  combined_tables[multi_cols] <- lapply(combined_tables[multi_cols], function(values) {
    gsub(values, pattern = ' *; *', replacement = ';')
  })

  # Replace any characters in ID columns that cannot appear in file names
  id_cols <- c('id', 'ref_group_id', 'report_id')
  replace_id_chars <- function(values) {
    gsub(values, pattern = invalid_id_char_pattern(), replacement = '_')
  }
  combined_tables[id_cols] <- lapply(id_cols, function(col) {
    if (col %in% multi_cols) {
      parts <- strsplit(combined_tables[[col]], split = ';')
      parts <- lapply(parts, replace_id_chars)
      return(vapply(parts, FUN.VALUE = character(1), function(x) {
        if (length(x) == 1 && is.na(x)) {
          return(x)
        } else {
          return(paste0(x, collapse = ';'))
        }
      }))
    } else {
      return(replace_id_chars(combined_tables[[col]]))
    }
  })

  return(combined_tables)
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
  modified_names <- gsub('[ -_]+', '_', modified_names)
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

  return(table)
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
    'part',
    'time',
    'year',
    'month',
    'day',
    'hour',
    'minute',
    'second',
    'link',
    'query_max',
    '_input_path_',
    '_sheet_name_',
    '_row_index_'
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
    report_id = 'miscellaneous',
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


#' Columns that can have multiple values
#'
#' Columns that can have multiple values separated by ';'.
#'
#' @return A character `vector`
#'
#' @examples
#' multi_input_columns()
#'
#' @export
multi_input_columns <- function() {
  c(
    'ref_group_id',
    'report_id',
    'input_type',
    'data_source',
    'color_by',
    'link'
  )
}
