% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ninetails_core_functions.R
\name{combine_gafs}
\alias{combine_gafs}
\title{Creates a two-dimensional array containing GASF and GADF resulting from
the transformation of a given ONT tail chunk.}
\usage{
combine_gafs(tail_chunk)
}
\arguments{
\item{tail_chunk}{numeric. A numeric vector representing signal chunk
within the analyzed dataset.}
}
\value{
an array (100,100,2) with values (coordinates) representing GASF
(first dimension) and GADF (second dimension) produced by the
\code{\link{create_gaf}} function applied to given fragment (tail chunk)
of analyzed ONT signal.
}
\description{
This function allows for the classification of signal fragments based
on angular coordinates generated by two methods (summation & difference)
simultaneously.
}
\details{
Using this approach increases the sensitivity of the classification.
It overcomes the limitations of each method.
}
\examples{
\dontrun{

ninetails::combine_gafs(tail_chunk = tail_chunk)

}
}
