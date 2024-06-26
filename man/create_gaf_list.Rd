% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ninetails_core_functions.R
\name{create_gaf_list}
\alias{create_gaf_list}
\title{Creates list of gramian angular matrices produced based on
list of splitted tails (tail chunks).}
\usage{
create_gaf_list(tail_chunk_list, num_cores)
}
\arguments{
\item{tail_chunk_list}{character string. The list object produced
by create_chunk_list function.}

\item{num_cores}{numeric [1]. Number of physical cores to use in processing
the data. Do not exceed 1 less than the number of cores at your disposal.}
}
\value{
A list of gaf matrices organized by the read ID_index
is returned. Always assign this returned list to a variable, otherwise
the long list will be printed to the console, which may crash your R session.
}
\description{
Creates list of gramian angular matrices produced based on
list of splitted tails (tail chunks).
}
\examples{
\dontrun{

gl <- ninetails::create_gaf_list(tail_chunk_list = tcl,
                                 num_cores = 2)

}
}
