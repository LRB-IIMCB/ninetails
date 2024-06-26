% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ninetails_data_postprocessing_functions.R
\name{count_nonA_abundance}
\alias{count_nonA_abundance}
\title{Counts reads with certain amount of nonA occurrences (instances)}
\usage{
count_nonA_abundance(residue_data, grouping_factor = NA)
}
\arguments{
\item{residue_data}{A dataframe or tibble containig non-A residue predictions
made by ninetails pipeline}

\item{grouping_factor}{grouping variable (e.g. "sample_name")}
}
\value{
A tibble with counts for each amount of nonA occurrences (single, two, more)
}
\description{
The function counts reads with single, two or more separate instances of nonA
occurrences as reported by ninetails main pipeline.
}
\examples{
\dontrun{
nonA_abundance <- ninetails::count_nonA_abundance(residue_data=residue_data,
                                                  grouping_factor="sample_name")
}
}
