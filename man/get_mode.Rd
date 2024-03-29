% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ninetails_misc_helper_functions.R
\name{get_mode}
\alias{get_mode}
\title{Calculates statistical mode of given vector.}
\usage{
get_mode(x, method = "density", na.rm = FALSE)
}
\arguments{
\item{x}{vector with values [numeric] to compute the mode for}

\item{method}{[character] string; which method of computing statistical mode
is meant to be used by the function. 2 options available: "value"/"density".
The latter is set by default.}

\item{na.rm}{logical [TRUE/FALSE] whether or not to remove NA values from the
calculation. Set to FALSE by default.}
}
\value{
statistical mode of given vector of values [numeric]
}
\description{
This function operates in 2 flavours. Either of them can be defined
by the "method" parameter. There are 2 options available: "density"
and "value".
}
\details{
If the "density" (default) option was chosen, then the
statistical mode (most frequent value) would be computed from the
normalized data density distribution. In this mode, the functon returns
a single value.

If the "value" option was selected, then the function return the mode based
on actual value. This means that, if the dataset was bi- or multimodal, the
returned vector would contain more than one most frequent values.

This function was written based on the following SO thread:
https://stackoverflow.com/questions/2547402/how-to-find-the-statistical-mode
Special thanks to Chris and hugovdberg!
}
\examples{
\dontrun{

test1 <- c(rep(2,5), rep(3,4), rep(1,4), rep(8,2), rep(7,3), rep(5,3)) # 1 most freq val
test2 <- c(rep(2,5), rep(3,4), rep(1,4), rep(8,5), rep(7,3), rep(5,3)) # 2 most freq val
test3 <- c(rep(2,5), rep(3,4), rep(1,4), rep(8,5), rep(7,5), rep(5,3)) # 3 most freq val
test4 <- c("Lorem", "ipsum", "dolor", "sit", "amet")

result <- ninetails::get_mode(x=test1, method= "value")
print(result)
class(result)

}

}
