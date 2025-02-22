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


#' @param read_all_sheets If `TRUE`, read all sheets in spreadsheet files like
#'   .ods and .xlsx
#'
#' @keywords internal
read_input_table <- function(input_paths, read_all_sheets = TRUE) {
  input_paths <- unique(input_paths)
  read_one <- function(path) {
    if (endsWith(path, '.csv')) {
      output <- list(read.csv(path, check.names = FALSE))
      names(output) <- path
    } else if (endsWith(path, '.tsv')) {
      output <- list(read.csv(path, check.names = FALSE, sep = '\t'))
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
