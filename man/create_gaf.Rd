% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ninetails_core_functions.R
\name{create_gaf}
\alias{create_gaf}
\title{Converts ONT signal to Gramian Angular Field.}
\usage{
create_gaf(tail_chunk, method = "s")
}
\arguments{
\item{tail_chunk}{numeric. A numeric vector representing signal chunk
within the analyzed dataset.}

\item{method}{character string specifying the type of Gramian Angular Field:
"s" can be used to produce summation field (GASF) and "d" to produce
difference field (GADF). Defaults to summation ["s"].}
}
\value{
an array (100,100,1) with values (coordinates) representing ONT signal.
}
\description{
This function represents time series data (ont squiggle) in a polar coordinate
system instead of the typical Cartesian coordinates. This is a Pythonic
pyts.image equivalent written in R language.
}
\details{
Two methods of such transformation are available: Gramian angular summation
field (GASF) and Gramian angular difference field (GADF).
}
\examples{
\dontrun{

ninetails::create_gaf(tail_chunk = tail_chunk,
                      method="s")

}
}
