% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ninetails_plotting_functions.R
\name{plot_nanopolish_qc}
\alias{plot_nanopolish_qc}
\title{Plots qc data (qc_tag) inherited from nanopolish polya function.}
\usage{
plot_nanopolish_qc(processing_info, frequency = TRUE)
}
\arguments{
\item{processing_info}{the output of nanopolish_qc function}

\item{frequency}{logical [TRUE/FALSE].}
}
\value{
a ggplot object
}
\description{
Provides insight into the proportion of each quality tag category per sample
or experiment condition defined by the user.
}
\details{
This is the ninetails' implementation of the
\code{\link[nanotail:plot_nanopolish_qc]{name}} function originally written
by P. Krawczyk (smeagol) and incorporated within the NanoTail package.

For original source code, see:
https://github.com/LRB-IIMCB/nanotail/blob/master/R/polya_plots.R

The variable names were adjusted according to the naming convention within
ninetails to avoid confusion.

Many thanks to the author of original source code for help and advice.
}
