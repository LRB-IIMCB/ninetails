#' Winsorizes nanopore signal - removes high cliffs (extending above & below
#' signal range). Slightly affects signal extremes.
#'
#' Based on the A. Signorell's DescTools https://github.com/AndriSignorell/DescTools/blob/master/R/DescTools.r
#'
#'
#' @param signal character string. Name of the given ONT signal.
#'
#' @return winsorized signal.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' winsorize_signal <- function(signal='12fdcb3-ewfd543-34552-1234ddta345')
#'
#' }
#'
winsorize_signal <- function(signal){

  #assertions

  if (missing(signal)) {
    stop("Signal vector is missing. Please provide a valid signal argument.", .call = FALSE)
  }

  assertthat::assert_that(assertive::is_numeric(signal), msg=paste("Signal vector must be numeric. Please provide a valid argument."))
  assertthat::assert_that(assertive::is_atomic(signal), msg=paste("Signal vector must be atomic. Please provide a valid argument."))

  signal_q <- stats::quantile(x=signal, probs=c(0.005, 0.995), na.rm=TRUE, type=7)
  minimal_val <- signal_q[1L]
  maximal_val <- signal_q[2L]

  signal[signal<minimal_val] <- minimal_val
  signal[signal>maximal_val] <- maximal_val

  winsorized_signal <- as.integer(signal)

  return(winsorized_signal)

}
