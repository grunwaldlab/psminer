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
  apply_validation(validate_ref_primary_usage)
  apply_validation(validate_ref_contextual_usage)
  apply_validation(validate_ploidy)
  apply_validation(validate_color_by)
  apply_validation(validate_date_cols)
  apply_validation(validate_source_illumina)
  apply_validation(validate_source_nanopore)
  apply_validation(validate_source_pacbio)
  apply_validation(validate_source_assembly)
  apply_validation(validate_source_image)
  apply_validation(validate_source_observation)
  apply_validation(validate_source_sra_query)
  apply_validation(validate_source_assembly_query)
  apply_validation(validate_source_ncbi_accession)
  apply_validation(validate_location_cols)


  # Save/print messages
  messages$description
  data.frame(
    data_id = metadata$data_id[messages$`_row_index_`],
    usage = metadata$usage[messages$`_row_index_`],
  )


  list(metadata = metadata, messages = messages)
}


#' @keywords internal
validate_date_cols <- function(metadata) {

  # If both `date` and date part columns are used, give preference to date part columns
  date_part_func <- list(
    year = lubridate::year,
    month = lubridate::month,
    day = lubridate::day,
    hour = lubridate::hour,
    minute = lubridate::minute,
    second = lubridate::second
  )
  has_parts <- apply(metadata[, names(date_part_func), drop = FALSE], MARGIN = 1, function(x) any(! is.na(x)))
  metadata$date[has_parts] <- apply(metadata[has_parts, names(date_part_func), drop = FALSE], MARGIN = 1, function(row) {
    paste0(row[!is.na(row)], collapse = ' ')
  })

  # Parse dates
  date_time_formats <- c(
    'Y-Om-d H:M:S',
    'Y-Om-d H:M',
    'Y-Om-d H',
    'Y-Om-d I:M:SOp',
    'Y-Om-d I:MOp',
    'Y-Om-d IOp',
    'Y-Om-d',
    'Y-Om',
    'Y',
    'Om/d/Y H:M:S',
    'Om/d/Y H:M',
    'Om/d/Y H',
    'Om/d/Y I:M:SOp',
    'Om/d/Y I:MOp',
    'Om/d/Y IOp',
    'Om/d/Y',
    'Om/Y'
  )
  parsed_dates <- lubridate::parse_date_time(metadata$date, orders = date_time_formats, quiet = TRUE)
  is_invalid <- is.na(parsed_dates) & ! is.na(metadata$date)
  messages <- make_warning_message(
    metadata,
    which(is_invalid),
    paste0('Samples do not have valid values for the `date` column.')
  )
  metadata$date <- lubridate::format_ISO8601(parsed_dates)

  # Split dates into columns for parts of date/time
  metadata[, names(date_part_func)] <- lapply(date_part_func, function(func) {
    func(parsed_dates)
  })

  list(metadata = metadata, messages = messages)
}

#' @keywords internal
validate_location_cols <- function(metadata) {
  # Parse location, latitude, and longitude columns
  loc_is_coord <- is_coordinate(metadata$location) & ! is.na(metadata$location)
  metadata$location[loc_is_coord] <- convert_coord_to_decimal_degrees(metadata$location[loc_is_coord])
  loc_is_address <- ! loc_is_coord & ! is.na(metadata$location)
  location_data <- tidygeocoder::geo(address = metadata$location[loc_is_address],
                                     full_results = TRUE, progress_bar = FALSE,
                                     custom_query = list("accept-language" = "en-US"))
  location_data <- unique(location_data[! is.na(location_data$lat), , drop = FALSE])
  metadata$latitude <- convert_coord_part_to_decimal_degrees(metadata$latitude)
  metadata$longitude <- convert_coord_part_to_decimal_degrees(metadata$longitude)

  # Use data from the location if latitude and longitude columns are missing values
  missing_long_lat <- is.na(metadata$latitude) | is.na(metadata$longitude)
  replace_with_coord_loc <- loc_is_coord & missing_long_lat
  split_loc <- strsplit(metadata$location, split = ', ')
  metadata$latitude[replace_with_coord_loc] <- vapply(split_loc[replace_with_coord_loc], FUN.VALUE = numeric(1), function(x) {
    as.numeric(x[1])
  })
  metadata$longitude[replace_with_coord_loc] <- vapply(split_loc[replace_with_coord_loc], FUN.VALUE = numeric(1), function(x) {
    as.numeric(x[2])
  })
  replace_with_lookup_loc <- missing_long_lat & metadata$location %in% location_data$address
  metadata$latitude[replace_with_lookup_loc] <- vapply(metadata$location[replace_with_lookup_loc], FUN.VALUE = numeric(1), function(x) {
    location_data$lat[location_data$address == x]
  })
  metadata$longitude[replace_with_lookup_loc] <- vapply(metadata$location[replace_with_lookup_loc], FUN.VALUE = numeric(1), function(x) {
    location_data$long[location_data$address == x]
  })

  # Infer location from latitude and longitude, giving preference to parsed version of user location
  loc_to_lookup <- ! is.na(metadata$latitude) & ! is.na(metadata$longitude)
  coord_data <- tidygeocoder::reverse_geo(lat = metadata$latitude[loc_to_lookup], long = metadata$longitude[loc_to_lookup],
                                          full_results = TRUE, custom_query = list("accept-language" = "en-US"),
                                          progress_bar = FALSE)
  parsed_locations <- location_data$display_name[match(metadata$location, location_data$address)]
  metadata$location[loc_to_lookup] <- coord_data$address
  metadata$location[! is.na(parsed_locations)] <- parsed_locations[! is.na(parsed_locations)]
  metadata$location[is.na(metadata$location)] <- ''

  # Add individual columns for major place name hierarchy elements
  grouped_level_key <- list(
    country = c('country'),
    region = c('state', 'region'),
    subregion = c('county', 'municipality'),
    place = c('city', 'town', 'village', 'hamlet'),
    part = c('district', 'city_district', 'subdistrict', 'place')
  )
  add_place_names <- apply(metadata[loc_to_lookup, names(grouped_level_key)], MARGIN = 1, function(x) {
    all(is.na(x))
  })
  metadata[loc_to_lookup, names(grouped_level_key)] <- lapply(seq_along(grouped_level_key), function(i) {
    cols <- grouped_level_key[[i]][grouped_level_key[[i]] %in% colnames(coord_data)]
    downloaded_names <- unlist(apply(coord_data[, cols], MARGIN = 1, simplify = FALSE, function(row) {
      row[! is.na(row)][1]
    }))
    ifelse(add_place_names, downloaded_names, metadata[[ names(grouped_level_key)[i]]][loc_to_lookup])
  })

  # Report problems
  is_invalid <- ! is_valid_lat(metadata$latitude) | ! is_valid_long(metadata$longitude)
  messages <- make_warning_message(
    metadata,
    which(is_invalid),
    paste0('Samples do not have valid values for the `latitude` or `longitude` columns.')
  )
  metadata$enabled[is_invalid] <- FALSE

  list(metadata = metadata, messages = messages)
}


#' @keywords internal
is_coordinate <- function(text) {
  grepl(text, pattern = '^[NSWE0-9,.\u00B0"\u201D\u201C\' -]+$') # \u00B0 = degree symbol
}


#' @keywords internal
split_coord <- function(text) {
  split_text <- strsplit(trimws(text), split = '[, ]+')
  lengths <- vapply(split_text, length, FUN.VALUE = numeric(1))
  is_odd <- lengths %% 2 != 0
  if (any(is_odd)) {
    warning('Could not split coordinates into equal parts. Replaceing with NA.')
    split_text[is_odd] <- rep(list(c(NA_character_, NA_character_)), sum(is_odd))
  }
  lapply(split_text, function(x) {
    c(
      paste0(x[1:(length(x) / 2)], collapse = ''),
      paste0(x[((length(x) / 2) + 1):length(x)], collapse = '')
    )
  })
}


#' @keywords internal
convert_coord_part_to_decimal_degrees <- function(coord_parts) {
  vapply(coord_parts, FUN.VALUE = numeric(1), function(coord_part) {
    if (is.na(coord_part) | coord_part == '') {
      return(as.numeric(coord_part))
    }
    coord_part <- trimws(coord_part)
    is_negative <- grepl(coord_part, pattern = '[SW-]+')
    coord_part <- gsub(coord_part, pattern = '[NWSE-]+', replacement = '')
    coord_part <- gsub(coord_part, pattern = '\u201D', replacement = '"') # replace fancy quotes
    coord_part <- gsub(coord_part, pattern = '\u201C', replacement = '"') # replace fancy quotes
    coord_part <- gsub(coord_part, pattern = "''", replacement = '"')
    if (grepl(coord_part, pattern = '^[0-9.]+$')) {
      output <- as.numeric(coord_part)
    } else if (grepl(coord_part, pattern = "^[0-9]{1,3}\u00B0[0-9.]+'$")) {
      subparts <- strsplit(coord_part, split = '\u00B0')[[1]]
      degrees <- as.numeric(subparts[1])
      minutes <- as.numeric(gsub(subparts[2], pattern = "'", replacement = ''))
      output <- degrees + minutes / 60
    } else if (grepl(coord_part, pattern = "^[0-9]{1,3}\u00B0[0-9.]+'[0-9.]+\"$")) {
      subparts <- strsplit(coord_part, split = "[\u00B0']")[[1]]
      degrees <- as.numeric(subparts[1])
      minutes <- as.numeric(subparts[2])
      seconds <- as.numeric(gsub(subparts[3], pattern = '"', replacement = ''))
      output <- degrees + minutes / 60 + seconds / 3600
    }
    if (is_negative) {
      output <- output * -1
    }
    return(output)
  })
}


#' @keywords internal
convert_coord_to_decimal_degrees <- function(coord_text) {
  parts <- lapply(split_coord(coord_text), convert_coord_part_to_decimal_degrees)
  vapply(parts, FUN.VALUE = character(1), function(x) {
    if (any(is.na(x))) {
      return(NA_character_)
    }
    paste0(x, collapse = ', ')
  })
}


#' @keywords internal
can_be_numeric <- function(values) {
  ! suppressWarnings(is.na(as.numeric(as.character(values))))
}


#' @keywords internal
is_valid_lat <- function(values) {
  vapply(values, FUN.VALUE = logical(1), function(x) {
    is.na(x) || (can_be_numeric(x) && as.numeric(x) >= -90 && as.numeric(x) <= 90)
  })
}


#' @keywords internal
is_valid_long <- function(values) {
  vapply(values, FUN.VALUE = logical(1), function(x) {
    is.na(x) || (can_be_numeric(x) && as.numeric(x) >= -180 && as.numeric(x) <= 180)
  })
}


#' @keywords internal
split_address <- function(values) {
  strsplit(values, split = ' *, *')
}


#' @keywords internal
validate_ploidy <- function(metadata) {
  # Assume haploid if missing
  is_missing <- is.na(metadata$ploidy)
  messages <- make_warning_message(
    metadata,
    which(is_missing),
    paste0('Value in `ploidy` column is missing. Haploid (1) is assumed. Set the ploidy to avoid this warning.')
  )
  metadata$ploidy[is_missing] <- 1

  # Check for invalid user inputs
  is_invalid <- ! grepl(metadata$ploidy, pattern = '[0-9]+')
  messages <- make_warning_message(
    metadata,
    which(is_invalid),
    paste0('Value in `ploidy` column does not appear to be an integer.')
  )
  metadata$enabled[is_invalid] <- FALSE

  list(metadata = metadata, messages = messages)
}


#' @keywords internal
validate_query_max <- function(metadata) {
  is_percentage <- endsWith(metadata$query_max, '%')
  metadata$query_max <- gsub(metadata$query_max, pattern = ' *%+$', replacement = '')
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
is_latin_binomial <- function(x) {
  grepl(x, pattern = '^[a-zA-Z]+ [a-zA-Z]+($| ).*$')
}

#' @keywords internal
validate_source_sra_query <- function(metadata, prefer_unique = TRUE, prefer_binomial = TRUE) {
  is_sra_query <- ! is.na(metadata$type) &  ! is.na(metadata$source) & metadata$type == 'NCBI SRA query' & metadata$enabled
  unique_query_data <- unique(metadata[is_sra_query, c('source', 'query_max'), drop = FALSE])
  ncbi_result <- lapply(seq_len(nrow(unique_query_data)), function(i) {
    get_ncbi_sra_runs(unique_query_data$source[i], max_count = unique_query_data$query_max[i])
  })
  new_sample_data <- do.call(rbind, lapply(which(is_sra_query), function(i) {
    query_data <- ncbi_result[unique_query_data$source == metadata$source[i] & unique_query_data$query_max == metadata$query_max[i]][[1]]
    output <- metadata[rep(i, nrow(query_data)), , drop = FALSE]
    if (is.na(metadata$data_id[i])) {
      output$data_id <- query_data$ncbi_acc
    } else {
      output$data_id <- paste0(metadata$data_id[i], query_data$ncbi_acc)
    }
    if (is.na(metadata$bio_id[i])) {
      output$bio_id <- query_data$biosample
    } else {
      output$bio_id <- paste0(metadata$bio_id[i], query_data$biosample)
    }
    if (is.na(metadata$name[i])) {
      output$name <- query_data$scientific_name
    } else {
      output$name <- paste0(metadata$name[i], query_data$scientific_name)
    }
    if (is.na(metadata$description[i])) {
      output$description <- query_data$title
    } else {
      output$description <- paste0(metadata$description[i], query_data$title)
    }
    output$type <- "NCBI accession"
    output$source <- query_data$ncbi_acc
    rownames(output) <- NULL
    output
  }))
  metadata$enabled[is_sra_query] <- FALSE
  metadata <- rbind(
    metadata,
    new_sample_data
  )
  list(metadata = metadata, messages = NULL)
}


#' Get metadata for SRA accessions using a query
#'
#' Get metadata for SRA accessions associated with a query to NCBI.
#'
#' @param query The query to use
#' @param max_count The maximum number or proportion of results to download
#' @param prefer_unique Give preference to diverse taxa, but still return
#'   duplicates if not enough unique taxa are found to satisfy `max_count`.
#' @param prefer_binomial Prefer genomes with a standard looking genus/species
#'   name (no numbers)
#'
#' @keywords internal
get_ncbi_sra_runs <- function(query, max_count = 10, prefer_unique = TRUE, prefer_binomial = TRUE) {
  # Search for SRA accession but don't download metadata yet
  search_result <- rentrez::entrez_search(db = 'sra', query, retmax = 10000, use_history = TRUE)

  # Convert max_count if a proportion is supplied instead of a count
  if (max_count < 1) {
    max_count <- ceiling(max_count * search_result$count)
  }

  # Parse 500 results at a time
  starts <- seq(from = 0 , to = length(search_result$ids) - 1, by = 500)
  results <- NULL
  preferred_results <- NULL
  for (start in starts) {
    summary_result <- rentrez::entrez_summary(db = 'sra', retmax = 500, retstart = start, web_history = search_result$web_history)
    if (length(search_result$ids) == 1) {
      summary_result <- list(summary_result)
    }
    run_data <- unlist(lapply(summary_result, function(x) {
      if (length(x$runs) > 1) {
        warning('Accession with multiple runs found. Only using first run.')
      }
      x$runs[1]
    }))
    expxml_data <- unlist(lapply(summary_result, function(x) {
      if (length(x$expxml) > 1) {
        warning('Accession with multiple runs found. Only using first run.')
      }
      x$expxml[1]
    }))
    batch_data <- data.frame(
      ncbi_acc = gsub(run_data, pattern = '.+ acc="(.+?)" .+', replacement = '\\1'),
      total_spots = gsub(run_data, pattern = '.+ total_spots="(.+?)" .+', replacement = '\\1'),
      total_bases = gsub(run_data, pattern = '.+ total_bases="(.+?)" .+', replacement = '\\1'),
      title = gsub(expxml_data, pattern = '.+<Title>(.+?)</Title>.+', replacement = '\\1'),
      instrument_model = gsub(expxml_data, pattern = '.+<Platform instrument_model="(.+?)">(.+?)</Platform>.+', replacement = '\\1'),
      platform = gsub(expxml_data, pattern = '.+<Platform instrument_model="(.+?)">(.+?)</Platform>.+', replacement = '\\2'),
      total_runs = gsub(expxml_data, pattern = '.+ total_runs="(.+?)" .+', replacement = '\\1'),
      total_size = gsub(expxml_data, pattern = '.+ total_size="(.+?)" .+', replacement = '\\1'),
      experiment_acc = gsub(expxml_data, pattern = '.+<Experiment acc="(.+?)" .+', replacement = '\\1'),
      study_acc = gsub(expxml_data, pattern = '.+<Study acc="(.+?)" .+', replacement = '\\1'),
      sample_acc = gsub(expxml_data, pattern = '.+<Sample acc="(.+?)" .+', replacement = '\\1'),
      taxid = gsub(expxml_data, pattern = '.+ taxid="(.+?)" .+', replacement = '\\1'),
      scientific_name = gsub(expxml_data, pattern = '.+ ScientificName="(.+?)"/>.+', replacement = '\\1'),
      library_strategy =  gsub(expxml_data, pattern = '.+<LIBRARY_STRATEGY>(.+?)</LIBRARY_STRATEGY>.+', replacement = '\\1'),
      library_source =  gsub(expxml_data, pattern = '.+<LIBRARY_SOURCE>(.+?)</LIBRARY_SOURCE>.+', replacement = '\\1'),
      library_selection =  gsub(expxml_data, pattern = '.+<LIBRARY_SELECTION>(.+?)</LIBRARY_SELECTION>.+', replacement = '\\1'),
      bioproject =  gsub(expxml_data, pattern = '.+<Bioproject>(.+?)</Bioproject>.+', replacement = '\\1'),
      biosample =  gsub(expxml_data, pattern = '.+<Biosample>(.+?)</Biosample>.+', replacement = '\\1')
    )
    rownames(batch_data) <- NULL
    batch_data <- batch_data[batch_data$library_strategy == 'WGS' & batch_data$library_source == 'GENOMIC' & batch_data$total_bases > 100000, , drop = FALSE]
    results <- rbind(results, batch_data)
    preferred_batch_data <- batch_data
    if (prefer_binomial) {
      preferred_batch_data <- preferred_batch_data[is_latin_binomial(preferred_batch_data$scientific_name), , drop = FALSE]
    }
    if (prefer_unique) {
      is_new <- ! duplicated(preferred_batch_data$taxid)
      if (! is.null(preferred_results)) {
        is_new <- is_new & ! preferred_batch_data$taxid %in% preferred_results$taxid
      }
      preferred_batch_data <- preferred_batch_data[is_new, , drop = FALSE]
    }
    preferred_results <- rbind(preferred_results, preferred_batch_data)
    if (nrow(preferred_results) >= max_count) {
      break
    }
  }

  # Subset results to query maximum
  if (nrow(preferred_results) >= max_count) {
    output <- preferred_results
  } else {
    output <- rbind(preferred_results, results[! results$ncbi_acc %in% preferred_results$ncbi_acc, , drop = FALSE])
  }
  output <- output[seq_len(min(c(nrow(output), max_count))), , drop = FALSE]
  rownames(output) <- NULL
  return(output)
}


#' @keywords internal
validate_source_assembly_query <- function(metadata, prefer_unique = TRUE, prefer_binomial = TRUE) {
  is_query <- ! is.na(metadata$type) &  ! is.na(metadata$source) & metadata$type == 'NCBI assembly query' & metadata$enabled
  unique_query_data <- unique(metadata[is_query, c('source', 'query_max'), drop = FALSE])
  ncbi_result <- lapply(seq_len(nrow(unique_query_data)), function(i) {
    get_ncbi_assemblies(unique_query_data$source[i], max_count = unique_query_data$query_max[i])
  })
  new_sample_data <- do.call(rbind, lapply(which(is_query), function(i) {
    query_data <- ncbi_result[unique_query_data$source == metadata$source[i] & unique_query_data$query_max == metadata$query_max[i]][[1]]
    output <- metadata[rep(i, nrow(query_data)), , drop = FALSE]
    if (is.na(metadata$data_id[i])) {
      output$data_id <- query_data$assemblyaccession
    } else {
      output$data_id <- paste0(metadata$data_id[i], query_data$assemblyaccession)
    }
    if (is.na(metadata$bio_id[i])) {
      output$bio_id <- query_data$biosampleaccn
    } else {
      output$bio_id <- paste0(metadata$bio_id[i], query_data$biosampleaccne)
    }
    if (is.na(metadata$name[i])) {
      output$name <- query_data$speciesname
    } else {
      output$name <- paste0(metadata$name[i], query_data$speciesname)
    }
    if (is.na(metadata$description[i])) {
      output$description <- query_data$organism
    } else {
      output$description <- paste0(metadata$description[i], query_data$organism)
    }
    output$type <- "NCBI accession"
    output$source <- query_data$assemblyaccession
    rownames(output) <- NULL
    output
  }))
  metadata$enabled[is_query] <- FALSE
  metadata <- rbind(
    metadata,
    new_sample_data
  )
  list(metadata = metadata, messages = NULL)
}


#' Get metadata for NCBI assemblies using a query
#'
#' Get metadata for NCBI assemblies associated with a query to NCBI.
#'
#' @param query The query to use
#' @param max_count The maximum number or proportion of results to download
#' @param prefer_unique Give preference to diverse taxa, but still return
#'   duplicates if not enough unique taxa are found to satisfy `max_count`.
#' @param prefer_binomial Prefer genomes with a standard looking genus/species
#'   name (no numbers)
#' @param prefer_complete Prefer genomes with "Complete Genome" or "Chromosome"
#'   assembly status
#' @param prefer_refseq Prefer genomes in RefSeq.
#' @param only_latest Only include the most recent version of a genome assembly
#'   if there are multiple versions in the result.
#'
#' @return A `data.frame`
#'
#' @keywords internal
get_ncbi_assemblies <- function(query, max_count = 10, prefer_unique = TRUE,
                                prefer_binomial = TRUE, prefer_complete = TRUE,
                                prefer_refseq = TRUE, only_latest = TRUE) {
  # Search for SRA accession but don't download metadata yet
  search_result <- rentrez::entrez_search(db = 'assembly', query, retmax = 10000, use_history = TRUE)

  # Convert max_count if a proportion is supplied instead of a count
  if (max_count < 1) {
    max_count <- ceiling(max_count * search_result$count)
  }

  # Parse 500 results at a time
  to_auto_extract <- c(
    "uid", "rsuid", "gbuid", "assemblyaccession", "lastmajorreleaseaccession",
    "latestaccession", "chainid", "assemblyname", "ucscname", "ensemblname",
    "taxid", "organism", "speciestaxid", "speciesname", "assemblytype",
    "assemblystatus", "wgs", "biosampleaccn", "biosampleid", "coverage",
    "partialgenomerepresentation", "primary", "assemblydescription",
    "releaselevel", "releasetype", "asmreleasedate_genbank", "asmreleasedate_refseq",
    "seqreleasedate", "asmupdatedate", "submissiondate", "lastupdatedate",
    "submitterorganization", "refseq_category", "fromtype", "annotrpturl",
    "ftppath_genbank", "ftppath_refseq", "ftppath_assembly_rpt",
    "ftppath_stats_rpt", "ftppath_regions_rpt", "sortorder", "meta"
  )
  starts <- seq(from = 0 , to = length(search_result$ids) - 1, by = 500)
  results <- NULL
  preferred_results <- NULL
  for (start in starts) {
    # Download a batch of assembly metadata
    summary_result <- rentrez::entrez_summary(db = 'assembly', retmax = 500, retstart = start, web_history = search_result$web_history)

    # Convert nested lists and XML to a table
    batch_data <- do.call(rbind, lapply(summary_result, function(x) {
      out <- as.data.frame(x[to_auto_extract])
      out$is_full_genome <- "full-genome-representation" %in% x$propertylist
      out$is_latest <- "latest" %in% x$propertylist
      out$bioprojects <- paste0(x$gb_bioprojects$bioprojectaccn, collapse = ';')
      out$exclfromrefseq <- paste0(unlist(x$exclfromrefseq), collapse = ';')
      out$sex <- x$biosource$sex
      out$isolate <- x$biosource$isolate
      out$infraspecies <- paste0(x$biosource$infraspecieslist$sub_value, collapse = ';')
      meta <- xml2::read_xml(paste0('<root>', out$meta, '</root>'))
      stats <- paste0(
        xml2::xml_attr(xml2::xml_find_all(meta, './/Stat'), 'category'),
        '_',
        xml2::xml_attr(xml2::xml_find_all(meta, './/Stat'), 'sequence_tag')
      )
      values <- xml2::xml_contents(xml2::xml_find_all(meta, './/Stat'))
      out[stats] <- as.list(as.character(values))
      out$assembly_level <- as.character(xml2::xml_contents(xml2::xml_find_all(meta, './assembly-level')))
      out$assembly_status <- as.character(xml2::xml_contents(xml2::xml_find_all(meta, './assembly-status')))
      out$representative_status <- as.character(xml2::xml_contents(xml2::xml_find_all(meta, './representative-status')))
      out$submitter_organization <- as.character(xml2::xml_contents(xml2::xml_find_all(meta, './submitter-organization')))
      out$meta <- NULL
      return(out)
    }))
    rownames(batch_data) <- NULL

    # Save all results for this batch
    results <- rbind(results, batch_data)

    # Remove assemblies that have a more recent accession in the results
    filter_for_latest <- function(table, all_results) {
      grouped <- split(all_results$assemblyaccession, sub(all_results$assemblyaccession, pattern = '\\.[0-9]+$', replacement = ''))
      latest <- vapply(grouped, FUN.VALUE = character(1), function(x) {
        versions <- as.numeric(sub(x, pattern = '^.+\\.([0-9]+)$', replacement = '\\1'))
        x[which.max(versions)]
      })
      table[table$assemblyaccession %in% latest, , drop = FALSE]
    }
    if (only_latest) {
      results <- filter_for_latest(results, results)
      batch_data <- filter_for_latest(batch_data, results)
    }

    # Create subset containing only preferred types of assemblies to see if there are enough yet
    preferred_batch_data <- batch_data
    if (prefer_complete) {
      is_complete <- function(table) {
        table$assemblystatus %in% c('Complete Genome', 'Chromosome') & table$partialgenomerepresentation == 'false'
      }
      preferred_batch_data <- preferred_batch_data[is_complete(preferred_batch_data), , drop = FALSE]
    }
    if (prefer_refseq) {
      is_refseq <- function(table) {
        startsWith(table$assemblyaccession, 'GCF_')
      }
      preferred_batch_data <- preferred_batch_data[is_refseq(preferred_batch_data), , drop = FALSE]
    }
    if (prefer_binomial) {
      preferred_batch_data <- preferred_batch_data[is_latin_binomial(preferred_batch_data$speciesname), , drop = FALSE]
    }
    if (prefer_unique) {
      is_new <- ! duplicated(preferred_batch_data$speciestaxid)
      if (! is.null(preferred_results)) {
        is_new <- is_new & ! preferred_batch_data$speciestaxid %in% preferred_results$speciestaxid
      }
      preferred_batch_data <- preferred_batch_data[is_new, , drop = FALSE]
    }
    preferred_results <- rbind(preferred_results, preferred_batch_data)
    if (nrow(preferred_results) >= max_count) {
      break
    }
  }

  # Subset results to query maximum
  if (nrow(preferred_results) >= max_count) {
    output <- preferred_results
  } else {
    remaining <- results[! results$assemblyaccession %in% preferred_results$assemblyaccession, , drop = FALSE]
    if (prefer_binomial || prefer_complete || prefer_refseq) {
      priority <- list(decreasing = TRUE)
      if (prefer_complete) {
        priority <- c(priority, list(is_complete(remaining)))
      }
      if (prefer_refseq) {
        priority <- c(priority, list(is_refseq(remaining)))
      }
      if (prefer_binomial) {
        priority <- c(priority, list(is_latin_binomial(remaining$speciesname)))
      }
      remaining <- remaining[do.call(order, priority), ]
    }
    output <- rbind(preferred_results, remaining)
  }
  output <- output[seq_len(min(c(nrow(output), max_count))), , drop = FALSE]
  rownames(output) <- NULL
  return(output)
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
    ends_with = valid_fasta_extentions()
  )
}

#' @keywords internal
validate_source_pacbio <- function(metadata) {
  validate_source_generic(
    metadata,
    type = 'Pacbio',
    must_exist = TRUE,
    max_count = 1,
    ends_with = valid_fastq_extentions()
  )
}

#' @keywords internal
validate_source_nanopore <- function(metadata) {
  validate_source_generic(
    metadata,
    type = 'Nanopore',
    must_exist = TRUE,
    max_count = 1,
    ends_with = valid_fastq_extentions()
  )
}

#' @keywords internal
validate_source_illumina <- function(metadata) {
  validate_source_generic(
    metadata,
    type = 'Illumina',
    must_exist = TRUE,
    max_count = 2,
    ends_with = valid_fastq_extentions()
  )
}


#' @keywords internal
valid_fastq_extentions <- function() {
  c('.fastq', '.fastq.gz', '.fq', '.fq.gz')
}

#' @keywords internal
valid_fasta_extentions <- function() {
  c('.fasta', '.fasta.gz', '.fa', '.fa.gz', '.fas', '.fas.gz')
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
      paste0('Invalid prefix for value in the `source` column. For `', type, '` input, one of the following prefixs are required: ',
             paste0(starts_with, collapse = ', '))
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
      paste0('Invalid suffix for value in the `source` column. For `', type, '` input, one of the following suffixes are required: ',
             paste0(ends_with, collapse = ', '))
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
  metadata$enabled <- is_enabled
  metadata$enabled[is_invalid] <- FALSE
  messages <- make_warning_message(
    metadata,
    which(is_invalid),
    'Invalid value in `enabled` column. Must be `TRUE`, `FALSE`, or left empty.'
  )
  list(metadata = metadata, messages = messages)
}


#' @keywords internal
validate_color_by <- function(metadata) {
  defualt_color_by <- c('usage')
  metadata$color_by <- vapply(strsplit(metadata$color_by, split = ';'), FUN.VALUE = character(1), function(parts) {
    paste0(unique(c(defualt_color_by, parts[! is.na(parts)])), collapse = ';')
  })
  validate_categorical(
    metadata,
    column = 'color_by',
    choices = colnames(metadata),
    clean = FALSE
  )
}


#' @keywords internal
validate_type <- function(metadata) {
  validate_categorical(
    metadata,
    column = 'type',
    choices = valid_types()
  )
}


#' @keywords internal
validate_usage <- function(metadata) {
  validate_categorical(
    metadata,
    column = 'usage',
    choices = valid_usages()
  )
}


#' @keywords internal
validate_ref_primary_usage <- function(metadata) {
  validate_categorical(
    metadata,
    column = 'ref_primary_usage',
    choices = valid_ref_usages()
  )
}


#' @keywords internal
validate_ref_contextual_usage <- function(metadata) {
  validate_categorical(
    metadata,
    column = 'ref_contextual_usage',
    choices = valid_ref_usages()
  )
}


#' @keywords internal
validate_categorical <- function(metadata, column, choices, clean = TRUE) {
  metadata[[column]] <- gsub(metadata[[column]], pattern = ' +', replacement = ' ')
  if (column %in% multi_input_columns()) {
    if (clean) {
      cleaned_parts <- lapply(strsplit(metadata[[column]], split = ';'), function(parts) {
        choices[match(tolower(parts), tolower(choices))]
      })
    } else {
      cleaned_parts <- strsplit(metadata[[column]], split = ';')
    }
    is_invalid <- vapply(cleaned_parts, FUN.VALUE = logical(1), function(parts) {
      any(is.na(parts))
    })
    cleaned_values <- vapply(cleaned_parts, FUN.VALUE = character(1), paste0, collapse = ';')
  } else {
    if (clean) {
      cleaned_values <- choices[match(tolower(metadata[[column]]), tolower(choices))]
    } else {
      cleaned_values <- metadata[[column]]
    }
    is_invalid <- is.na(cleaned_values)
  }
  messages <- make_warning_message(
    metadata,
    which(is_invalid),
    paste0('Invalid value in `', column, '` column. Must be one of: ', paste0(choices, collapse = ', '))
  )
  metadata$enabled[is_invalid] <- FALSE
  metadata[[column]] <- cleaned_values
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
#' @keywords internal
read_input_tables <- function(input_paths, read_all_sheets = TRUE, add_input_metadata = TRUE) {

  # Read data from whatever format it is stored in
  input_paths <- unique(input_paths)
  read_one <- function(path) {
    if (endsWith(path, '.csv')) {
      output <- list(tibble::as_tibble(utils::read.csv(path, check.names = FALSE)))
      names(output) <- path
      sheets <- NA
    } else if (endsWith(path, '.tsv')) {
      output <- list(tibble::as_tibble(utils::read.csv(path, check.names = FALSE, sep = '\t')))
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
  id_cols <- c('data_id', 'bio_id', 'report_id')
  replace_id_chars <- function(values) {
    out <- gsub(values, pattern = invalid_id_char_pattern(), replacement = '_')
    gsub(out, pattern = '_+', replacement = '_') # Some programs replace runs of underscores
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
      'The following columns in an input CSV have values but no header: ',
      paste0(which(is_headerless & ! is_empty), collapse = ', ')
    ))
  }
  table <- table[, ! is_headerless, drop = FALSE]

  # Check for duplicated columns
  present_known_cols <- colnames(table)[colnames(table) %in% known_columns]
  duplicated_cols <- unique(present_known_cols[duplicated(present_known_cols)])
  if (length(duplicated_cols) > 0) {
    stop(call. = FALSE,
         'The following columns occur more than once in an input CSV: ',
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
    # 'time_zone',
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
#' @export
default_column_values <- function() {
  c(
    usage = 'sample',
    report_id = 'miscellaneous',
    enabled = TRUE,
    query_max = 10,
    ref_primary_usage = 'optional',
    ref_contextual_usage = 'optional'
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


#' Values for the `ref_*_usage` columns
#'
#' Valid values for the `ref_primary_usage` and `ref_contextual_usage` columns
#' of `pathogensurveillance` input tables.
#'
#' @return A character `vector`
#'
#' @examples
#' valid_ref_usages()
#'
#' @export
valid_ref_usages <- function() {
  c(
    'optional',
    'required',
    'excluded',
    'exclusive'
  )
}
