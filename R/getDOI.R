
#' Get BibTeX Citation from DOI
#'
#' @author Ricardo I. Alcala
#' @title Get BibTeX Citation from DOI
#' @description This function retrieves BibTeX citation from a given DOI using the rcrossref package.
#'
#' @param doi The DOI (Digital Object Identifier) for the publication.
#' @return A character string containing the BibTeX citation.
#'
#' @examples
#' \dontrun{
#' doi <- "https://doi.org/10.1002/cpbi.102"
#' bibtex_citation <- get_bibtex_from_doi(doi)
#' cat(bibtex_citation)
#' }
#'
#' @import rcrossref
#'
#' @keywords rcrossref DOI BibTeX
#'
#' @export

# function
get_bibtex_from_doi <- function(doi, ...) {
  # x <- cr_works(doi = doi)
  if (length(doi) > 0) {
  bibtex_citation <- cat(rcrossref::cr_cn(dois = doi, ...))
    return(bibtex_citation)
  } else {
    return(NULL)
  }
}
