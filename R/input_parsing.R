#' Parse one or more inputs to `pathogensurveillance`
#'
#' Parse and validate one or more input table of sample/reference metadata to
#' the `pathogensurveillance` pipeline and puts the parsed data in a single
#' table.
#'
#' @param input_paths The paths to one or more tables in TSV, CSV, or ODS
#'   (LibreOffice Calc), xls, or xlsx files (Excel) format.
#' @param fail_on_warning If `TRUE`, things that will normally cause a warning
#'   will instead cause an error.
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
parse_path_surveil_input <- function(input_paths, fail_on_warning = FALSE) {

  # Read and standardize all inputs into a single table without any data interpretation
  metadata <- read_input_tables(input_paths)

  # Make standard way to apply changes made by validation functions
  messages <- make_warning_message(metadata, numeric(0), 'Not used.')
  apply_validation <- function(validation_func) { # NOTE: causes side effects, unlike most R code
    result <- validation_func(metadata)
    metadata <<- result$metadata
    messages <<- rbind(messages, result$messages)
  }

  # Apply validation functions
  # NOTE: These modify `metadata` and `messages` and the order matters
  apply_validation(validate_enabled)
  apply_validation(validate_usage)
  apply_validation(validate_type)
  apply_validation(validate_query_max)
  apply_validation(validate_source_illumina)
  apply_validation(validate_source_nanopore)
  apply_validation(validate_source_pacbio)
  apply_validation(validate_source_assembly)
  apply_validation(validate_source_ncbi_accession)
  apply_validation(validate_source_sra_query)
  apply_validation(validate_source_assembly_query)
  apply_validation(validate_source_image)
  apply_validation(validate_source_observation)

  return(metadata)
}

#' @keywords internal
validate_query_max <- function(metadata) {
  can_be_numeric <- function(values) {
    ! suppressWarnings(is.na(as.numeric(as.character(values))))
  }
  is_percentage <- endsWith(metadata$query_max, '%')
  metadata$query_max <- gusb(metadata$query_max, pattern = ' *%+$', replacement = '')
  is_invalid <- ! can_be_numeric(metadata$query_max)
  messages <- make_warning_message(
    metadata,
    which(is_invalid),
    paste0('Value in `query_max` column does not appear to be an integer, proportion, or percentage (i.e. ending with %).')
  )
  metadata$enabled[is_invalid] <- FALSE
  metadata$query_max[is_invalid] <- NA
  metadata$query_max <- as.numeric(metadata$query_max)
  metadata$query_max[! is_invalid & is_percentage] <- metadata$query_max[! is_invalid & is_percentage] / 100
  list(metadata = metadata, messages = messages)
}

#' @keywords internal
validate_source_observation <- function(metadata) {
  validate_source_generic(
    metadata,
    type = 'observation',
    must_exist = FALSE,
    min_count = NULL
  )
}

#' @keywords internal
validate_source_assembly_query <- function(metadata) {
  validate_source_generic(
    metadata,
    type = 'NCBI assembly query',
    must_exist = TRUE,
    min_count = NULL
  )
}


#' @keywords internal
validate_source_sra_query <- function(metadata) {
  validate_source_generic(
    metadata,
    type = 'NCBI SRA query',
    must_exist = TRUE,
    min_count = NULL
  )
}


#' @keywords internal
validate_source_ncbi_accession <- function(metadata) {
  validate_source_generic(
    metadata,
    type = 'NCBI accession',
    must_exist = TRUE,
    max_count = NULL,
    starts_with = c(
      'ERR', 'SRR', 'DRR',    # SRA runs
      'SAMD', 'SAME', 'SAMN', # Biosamples
      'GCA_', 'GCF_'          # Assemblies
    )
  )
}


#' @keywords internal
validate_source_image <- function(metadata) {
  validate_source_generic(
    metadata,
    type = 'image',
    must_exist = TRUE,
    max_count = NULL,
    ends_with = c('.png', '.jpg', '.jpeg', '.tiff', '.webp')
  )
}


#' @keywords internal
validate_source_assembly <- function(metadata) {
  validate_source_generic(
    metadata,
    type = 'Assembly',
    must_exist = TRUE,
    max_count = 1,
    ends_with = c('.fastq', '.fastq.gz', 'fq', 'fq.gz')
  )
}

#' @keywords internal
validate_source_pacbio <- function(metadata) {
  validate_source_generic(
    metadata,
    type = 'Pacbio',
    must_exist = TRUE,
    max_count = 1,
    ends_with = c('.fastq', '.fastq.gz', 'fq', 'fq.gz')
  )
}

#' @keywords internal
validate_source_nanopore <- function(metadata) {
  validate_source_generic(
    metadata,
    type = 'Nanopore',
    must_exist = TRUE,
    max_count = 1,
    ends_with = c('.fastq', '.fastq.gz', 'fq', 'fq.gz')
  )
}

#' @keywords internal
validate_source_illumina <- function(metadata) {
  validate_source_generic(
    metadata,
    type = 'Illumina',
    must_exist = TRUE,
    max_count = 2,
    ends_with = c('.fastq', '.fastq.gz', 'fq', 'fq.gz')
  )
}


#' Check that `source` counts and values
#'
#' Check that `source` has the right number of values and starts/ends with
#' the right values.
#'
#' @param metadata the input table
#' @param type the value of the `type` column to validate
#' @param must_exist If `TRUE`, it cannot be `NA` or blank. This check does not
#'   consider the presence of `;` so it is safe for values where `;` is not used
#'   a separater.
#' @param min_count The minimum number of values that can be present. Set to
#'   `NULL` to disable this check.
#' @param max_count The maximum number of values that can be present. Set to
#'   `NULL` to disable this check.
#' @param starts_with A `vector` of values that a valid value can have as a
#'   prefix. Set to `NULL` to disable this check.
#' @param ends_with A `vector` of values that a valid value can have as a
#'   suffix. Set to `NULL` to disable this check.
#' @param default What is put if a value is missing.
#'
#' @keywords internal
validate_source_generic <- function(metadata, type, must_exist = TRUE, min_count = 1, max_count = NULL, starts_with = NULL, ends_with = NULL, default = NULL) {
  # Make empty message table
  messages <- make_warning_message(metadata, index = NULL, message = 'not used')

  # Check for missing values
  subset_meta <- metadata[! is.na(metadata$type) & metadata$type == type, , drop = FALSE]
  if (! must_exist) {
    is_missing <- rep(FALSE, nrow(subset_meta))
  } else {
    is_missing <- is.na(subset_meta$source) | subset_meta$source == ''
    messages <- rbind(messages, make_warning_message(
      subset_meta,
      which(is_missing),
      paste0('Missing or value in `source` column. For `', type, '` input a value must be supplied.')
    ))
  }

  # Check for too few values
  source_parts <- strsplit(subset_meta$source, split = ';')
  source_counts <- vapply(source_parts, FUN.VALUE = numeric(1), length)
  if (is.null(min_count) || min_count <= 0) {
    is_too_few <- rep(FALSE, nrow(subset_meta))
  } else {
    is_too_few <- is.na(subset_meta$source) | source_counts < min_count
    messages <- rbind(messages, make_warning_message(
      subset_meta,
      which(is_too_few),
      paste0('Too few values in `source` column. For `', type, '` input, at least ', min_count, ' values must be supplied.')
    ))
  }

  # Check for too many values
  if (is.null(max_count)) {
    is_too_many <- rep(FALSE, nrow(subset_meta))
  } else {
    is_too_many <- source_counts > max_count
    messages <- rbind(messages, make_warning_message(
      subset_meta,
      which(is_too_many),
      paste0('Too many values in the `source` column. For `', type, '` input, at most ', max_count, ' value(s) must be supplied.')
    ))
  }

  # Check for invalid starts
  if (is.null(starts_with)) {
    is_invalid_start <- rep(FALSE, nrow(subset_meta))
  } else {
    is_invalid_start <- vapply(tolower(subset_meta$source), FUN.VALUE = logical(1), function(path) {
      ! is.na(path) && all(! startsWith(path, starts_with))
    })
    messages <- rbind(messages, make_warning_message(
      subset_meta,
      which(is_invalid_start),
      pasteo('Invalid prefix for value in the `source` column. For `', type, '` input, one of the following prefixs are required: ',
             pasteo(starts_with, collapse = ', '))
    ))
  }

  # Check for invalid endings
  if (is.null(ends_with)) {
    is_invalid_end <- rep(FALSE, nrow(subset_meta))
  } else {
    is_invalid_end <- vapply(tolower(subset_meta$source), FUN.VALUE = logical(1), function(path) {
      ! is.na(path) && all(! endsWith(path, ends_with))
    })
    messages <- rbind(messages, make_warning_message(
      subset_meta,
      which(is_invalid_end),
      pasteo('Invalid suffix for value in the `source` column. For `', type, '` input, one of the following suffixes are required: ',
             pasteo(ends_with, collapse = ', '))
    ))
  }

  is_invalid <- is_missing | is_too_few | is_too_many | is_invalid_start | is_invalid_end
  metadata$enabled[metadata$`_row_index_` %in% subset_meta$`_row_index_`[is_invalid]] <- FALSE
  list(metadata = metadata, messages = messages)
}


#' @keywords internal
validate_enabled <- function(metadata) {
  is_enabled <- as.logical(metadata$enabled)
  is_invalid <- is.na(is_enabled)
  metadata$enabled <- as.logical(metadata$enabled)
  metadata$enabled[! is_invalid] <- FALSE
  messages <- make_warning_message(
    metadata,
    which(is_invalid),
    'Invalid value in `enabled` column. Must be `TRUE`, `FALSE`, or left empty.'
  )
  list(metadata = metadata, messages = messages)
}

#' @keywords internal
validate_type <- function(metadata) {
  valid_values <- valid_types()
  metadata$type <- gsub(metadata$type, pattern = ' +', replacement = ' ')
  cleaned_values <- valid_values[match(tolower(metadata$type), tolower(valid_values))]
  is_invalid <- is.na(cleaned_values)
  metadata$enabled[is_invalid] <- FALSE
  messages <- make_warning_message(
    metadata,
    which(is_invalid),
    paste0('Invalid value in `type` column. Must be one of: ', paste0(valid_values, collapse = ', '))
  )
  metadata$type <- cleaned_values
  list(metadata = metadata, messages = messages)
}

#' @keywords internal
validate_usage <- function(metadata) {
  metadata$usage <- gsub(metadata$usage, pattern = ' +', replacement = ' ')
  is_invalid <- ! metadata$usage %in% valid_usages()
  metadata$enabled[! is_invalid] <- FALSE
  messages <- make_warning_message(
    metadata,
    which(is_invalid),
    paste0('Invalid value in `type` column. Must be one of: ', paste0(valid_usages(), ', '))
  )
  list(metadata = metadata, messages = messages)
}

#' @keywords internal
make_warning_message <- function(metadata, index, message) {
  if (length(index) == 0) {
    return(
      data.frame(check.names = FALSE,
        '_row_index_' = numeric(0),
        message_type = character(0),
        description = character(0)
      )
    )
  } else {
    return(
      data.frame(check.names = FALSE,
        '_row_index_' = metadata$`_row_index_`[index],
        message_type = 'WARNING',
        description = message
      )
    )
  }
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
    'data_id',
    'bio_id',
    'report_id',
    'usage',
    'name',
    'description',
    'enabled',
    'type',
    'source',
    'query_max',
    'reference',
    'ref_group_id',
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
    'date',
    'year',
    'month',
    'day',
    'hour',
    'minute',
    'second',
    'link',
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
    usage = 'sample',
    report_id = 'miscellaneous',
    enabled = TRUE,
    query_max = 10,
    ref_primary_usage = 'optional',
    ref_contextual_usage = 'optional',
    color_by = 'usage'
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
    'usage',
    'source',
    'color_by',
    'link'
  )
}


#' Values for the `type` column
#'
#' Valid values for the `type` column of `pathogensurveillance` input
#' tables.
#'
#' @return A character `vector`
#'
#' @examples
#' valid_types()
#'
#' @export
valid_types <- function() {
  c(
    'Illumina',
    'Nanopore',
    'Pacbio',
    'Assembly',
    'NCBI accession',
    'NCBI SRA query',
    'NCBI assembly query',
    'image',
    'observation'
  )
}


#' Values for the `usage` column
#'
#' Valid values for the `usage` column of `pathogensurveillance` input
#' tables.
#'
#' @return A character `vector`
#'
#' @examples
#' valid_usages()
#'
#' @export
valid_usages <- function() {
  c(
    'sample',
    'reference',
    'sample metadata',
    'reference metadata'
  )
}
