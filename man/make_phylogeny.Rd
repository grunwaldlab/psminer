% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_core_gene_phylogeny.R, R/make_trees.R
\name{make_phylogeny}
\alias{make_phylogeny}
\title{Make core gene phylogeny}
\usage{
make_phylogeny(
  input,
  sample_data,
  ref_data,
  interactive = knitr::is_html_output()
)

make_phylogeny(
  input,
  sample_data,
  ref_data,
  interactive = knitr::is_html_output()
)
}
\arguments{
\item{input}{Either a file path to a .treefile or a "phylo" object}

\item{interactive}{}
}
\value{
Core gene phylogeny

Core gene phylogeny
}
\description{
The tree will either be made using ggtree, if output is a PDF, or phylo

The tree will either be made using ggtree, if output is a PDF, or phylo
}
