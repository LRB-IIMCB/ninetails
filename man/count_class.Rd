% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ninetails_data_postprocessing_functions.R
\name{count_class}
\alias{count_class}
\title{Counts read classes found in read_classes dataframe produced
by ninetails pipeline.}
\usage{
count_class(class_data, grouping_factor = NA, detailed = TRUE)
}
\arguments{
\item{class_data}{A dataframe or tibble containig read_classes predictions
made by ninetails pipeline}

\item{grouping_factor}{character string. A grouping variable
(e.g. "sample_name")}

\item{detailed}{logical [TRUE/FALSE]. If TRUE, the counts will be provided
based on the "comments" column, which contains detailed information on the
assigned class. If FALSE, the counts will be provided based on "class" column
which gives more crude glimpse on the classification - i.e. provides an info
whether the reads were considered as "decorated", "blank" and
"unclassified" only. By default, the TRUE option is set.}
}
\value{
A tibble with counts for each non-A residue
}
\description{
Process the information returned by ninetails found in the read_classes
dataframe in the prediction column.
}
\examples{
\dontrun{
class_counted <- ninetails::count_class(class_data=out[[1]],
                                        grouping_factor=NA,
                                        detailed=TRUE)
}
}
