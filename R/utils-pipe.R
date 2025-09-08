# Package imports and utility functions
#
# This file contains package imports and the pipe operator

# Tidyverse imports
#' @import dplyr
#' @import forcats
#' @import purrr
#' @import stringr
#' @import tibble
#' @import tidyr
#' @import tidyselect

# Data handling
#' @importFrom data.table fread as.data.table
#' @import vroom

# Parallel processing
#' @import doSNOW
#' @importFrom foreach %dopar%
#' @import parallel

# Statistics and utilities
#' @importFrom stats median sd
#' @importFrom utils head tail
#' @import rlang

# Bioconductor
#' @import Rsamtools
#' @import S4Vectors

# Input validation
#' @importFrom assertthat assert_that
#' @importFrom checkmate test_string

# CLI and logging
#' @import cli

# Machine learning
#' @import keras

# Additional utilities
#' @import magrittr
#' @import rrapply
#' @importFrom reshape2 melt

#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL
