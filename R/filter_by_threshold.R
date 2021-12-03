#' Peak signal detection in realtime timeseries data (nanopore signal).
#' An R implementation of code originally designed for Matlab (see link below).
#' Herein adjusted for nanopore data.
#'
#' @param chunkname character string. Name of the given signal chunk within the
#' analyzed dataset.
#'
#' @param tail_chunk_list character string. Full path of the list object produced
#' by create_chunk_list function.
#'
#' @param lag numeric [1] the lag of the moving window (amount of observations
#' to be taken into consideration for filtering)
#'
#' @param threshold numeric [1] the z-score at which the algorithm signals (how many sd
#' away from the actual signal the threshold is set).
#'
#' @param influence numeric [1] the influence (between 0 and 1) of new signals on the mean and standard deviation
#'
#' @return a list of signals with regions falling outside the threshold.
#' @export
#'
#' @examples
#'\dontrun{
#' filter_by_threshold(chunkname ='1234-trew-t33451-123455_1',
#'                     tail_chunk_list=tail_chunk_list,lag = 50,
#'                     threshold = 3, influence=0)
#'}

# thresholding algorithm based on this thread:
#https://stackoverflow.com/questions/22583391/peak-signal-detection-in-realtime-timeseries-data/54507329#54507329


filter_by_threshold <- function(chunkname, tail_chunk_list, lag, threshold, influence){

  signal <- tail_chunk_list[[chunkname]]
  signals <- rep(0,length(signal))
  filtered_signal <- signal[0:lag]
  avg_filter <- NULL
  std_filter <- NULL
  avg_filter[lag] <- mean(signal[0:lag], na.rm=TRUE)
  std_filter[lag] <- sd(signal[0:lag], na.rm=TRUE)
  for (i in (lag+1):length(signal)){
    if (abs(signal[i]-avg_filter[i-1]) > threshold*std_filter[i-1]) {
      if (signal[i] > avg_filter[i-1]) {
        signals[i] <- 1;
      } else {
        signals[i] <- -1;
      }
      filtered_signal[i] <- influence*signal[i] + (1-influence)*filtered_signal[i-1]
    } else {
      signals[i] <- 0
      filtered_signal[i] <- signal[i]
    }
    avg_filter[i] <- mean(filtered_signal[(i-lag):i], na.rm=TRUE)
    std_filter[i] <- sd(filtered_signal[(i-lag):i], na.rm=TRUE)
  }
  return(signals)

}
