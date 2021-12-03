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


  signal_q <- stats::quantile(x=signal, probs=c(0.005, 0.995), na.rm=TRUE, type=7)
  minval <- signal_q[1L]
  maxval <- signal_q[2L]

  signal[signal<minval] <- minval
  signal[signal>maxval] <- maxval

  winsorized_signal <- as.integer(signal)

  return(winsorized_signal)

}
