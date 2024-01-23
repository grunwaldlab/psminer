## `Psminer`, an R package for analysis of the nf-core/PathogenSurveillance pipeline

`psminer` is an R package with functions that can read, summarize, plot, and manipulate data produced by the nextflow pipeline [`pathogensurveillance`](https://github.com/nf-core/pathogensurveillance). 

## Installation

You can also install the development version for the newest features,
bugs, and bug fixes in R as follows:

``` r
install.packages("devtools")
devtools::install_github("grunwaldlab/psminer")
```

``` R
# Replace "your_doi_here" with the actual DOI
doi <- "https://doi.org/10.1002/cpbi.102"
bibtex_citation <- get_bibtex_from_doi(doi, format = "bibtex")
```

## License

This work is subject to the [MIT
License](https://github.com/grunwaldlab/metacoder/blob/master/LICENSE).
