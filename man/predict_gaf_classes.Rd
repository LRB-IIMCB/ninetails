% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ninetails_core_functions.R
\name{predict_gaf_classes}
\alias{predict_gaf_classes}
\title{Performs classification of matrices produced from read chunks with CNN.}
\usage{
predict_gaf_classes(gaf_list)
}
\arguments{
\item{gaf_list}{[list] A list of gaf matrices organized by the read ID_index.}
}
\value{
a list of gafs predictions based on used model.
}
\description{
This function allows to assign gafs corresponding to the given signals
into one of 4 categories (A, C, G, U, respectively). This function in its
current implementation allows to load the pretrained CNN model. It uses
tensorflow backend to produce the classification output.
}
\examples{
\dontrun{

pl <- ninetails::predict_gaf_classes(gl)

}
}
