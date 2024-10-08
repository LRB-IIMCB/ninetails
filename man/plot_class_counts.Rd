% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ninetails_plotting_functions.R
\name{plot_class_counts}
\alias{plot_class_counts}
\title{Plotting read classes data per category assigned to the analyzed reads.}
\usage{
plot_class_counts(
  class_data,
  grouping_factor = NA,
  frequency = TRUE,
  type = "R"
)
}
\arguments{
\item{class_data}{A dataframe or tibble containing read_classes predictions
made by ninetails pipeline}

\item{grouping_factor}{character string. A grouping variable (e.g. "sample_name")}

\item{frequency}{logical [TRUE/FALSE]. If TRUE, the frequency will be plotted.
If FALSE, raw counts will be shown. This parameter is set to TRUE by default.}

\item{type}{character string ["R"/"N"/"A"]. This variable controls the level
of detail of the resulting plot:
\itemize{
\item "R" - detailed classification (based on column "comments")
\item "N" - crude classification (based on column "class")
\item "A" - reads decorated with non-As exclusively
}
By default, the "R" option is set.}
}
\value{
ggplot object with read class prediction
}
\description{
This function requires as input the dataframe with read_classes provided by
ninetails pipeline. Function works in 3 flavours, by plotting either: \itemize{
\item detailed classification (based on column "comments")
\item crude classification (based on column "class")
\item reads decorated with non-As exclusively
}
}
\details{
Function  based on the Nanotail equivalent: https://github.com/LRB-IIMCB/nanotail
Many thanks to smeagol (Pawel Krawczyk) for advice & support!
}
\examples{
\dontrun{

ninetails::plot_class_counts(class_data=read_classes_dataframe,
                             grouping_factor="sample_name",
                             frequency=TRUE,
                             type="R")
}


}
