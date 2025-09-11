#' Extract poly(A) tail signal segments from POD5 files
#'
#' The function uses \pkg{reticulate} to import the Python library `pod5`. For each pod5 file referenced in
#' `polya_data` it opens a `pod5.Reader`, enumerates reads, matches `read_id` and extracts the numeric
#' slice `signal[polya_start:polya_end]` whenever indices are valid. If a read is missing or indices are
#' invalid, `polya_signal` for that row will be `numeric(0)` and `signal_length = 0`.
#'
#' The function performs several input checks and will stop early with informative messages when
#' required arguments or columns are missing or have incorrect types.
#'
#' @param polya_data A data frame with at least the following columns:
#' \describe{
#'  \item{read_id}{Character: ONT read identifiers.}
#'  \item{pod5_file}{Character: corresponding POD5 filename (relative to `pod5_dir`).}
#'  \item{polya_start}{Integer/numeric: start index of the poly(A) region in the signal.}
#'  \item{polya_end}{Integer/numeric: end index of the poly(A) region in the signal.}
#' }
#'
#' @param pod5_dir Character string. Path to the directory containing POD5 files referenced in
#' `polya_data$pod5_file`. Equivalent to guppy workspace.
#'
#' @returns A tibble combining the original rows from `polya_data` with:
#' \describe{
#'  \item{polya_signal}{list-column: each element is an atomic numeric vector
#'  with the extracted signal (length may be zero if extraction failed or read missing).}
#'  \item{signal_length}{integer: length of the extracted `polya_signal`.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' res <- extract_tails_from_pod5(polya_data, "/path/to/pod5_dir")
#' }
#'
#'
extract_tails_from_pod5 <- function(polya_data, pod5_dir) {

  ## variable bindings to satisfy R CMD check
  i <- reader <- NULL

  # Assertions
  assertthat::assert_that(is.character(pod5_dir), dir.exists(pod5_dir),
                          msg = "pod5_dir must be a valid directory path (character string and directory must exist)")

  assertthat::assert_that(is.data.frame(polya_data) && nrow(polya_data) > 0,
                          msg = "polya_data must be a non-empty data frame")

  required_cols <- c("read_id", "pod5_file", "polya_start", "polya_end")
  missing_cols <- setdiff(required_cols, colnames(polya_data))
  if (length(missing_cols) > 0) {
    stop(sprintf("polya_data is missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  assertthat::assert_that(is.character(polya_data$read_id),
                          msg = "polya_data$read_id must be a character vector")
  assertthat::assert_that(is.character(polya_data$pod5_file),
                          msg = "polya_data$pod5_file must be a character vector")
  assertthat::assert_that(is.numeric(polya_data$polya_start),
                          msg = "polya_data$polya_start must be numeric")
  assertthat::assert_that(is.numeric(polya_data$polya_end),
                          msg = "polya_data$polya_end must be numeric")

  # ensure pod5 python module is available before import (clear message)
  if (!reticulate::py_module_available("pod5")) {
    stop("Python module 'pod5' is not available in the current environment.")
  }

  # expose pod5 to R via reticulate
  pod5 <- reticulate::import("pod5")


  polya_by_file <- split(polya_data, polya_data$pod5_file)

  results <- lapply(names(polya_by_file), function(pod5_file) {
    current_data <- polya_by_file[[pod5_file]]
    pod5_path <- file.path(pod5_dir, pod5_file)

    if (!file.exists(pod5_path)) {
      warning(sprintf("POD5 file not found: %s", pod5_path))
      return(NULL)
    }

    tryCatch({
      reader <- pod5$Reader(pod5_path)
      all_reads <- reticulate::iterate(reader$reads())
      available_ids <- sapply(all_reads, function(r) trimws(r$read_id))
      reads_by_id <- stats::setNames(all_reads, available_ids)

      result <- lapply(seq_len(nrow(current_data)), function(i) {
        read_id <- trimws(current_data$read_id[i])
        if (!read_id %in% names(reads_by_id)) {
          message(sprintf("Read ID not found: %s in file %s", read_id, pod5_file))
          return(list(
            read_id = read_id,
            pod5_file = pod5_file,
            polya_signal = numeric(0),
            signal_length = 0
          ))
        }
        signal <- reticulate::py_to_r(reads_by_id[[read_id]]$signal)
        start_idx <- current_data$polya_start[i]
        end_idx <- current_data$polya_end[i]
        # Diagnostic print
        message(sprintf("Extracting: %s, signal length: %d, start: %d, end: %d", read_id, length(signal), start_idx, end_idx))
        if (is.numeric(start_idx) && is.numeric(end_idx) && start_idx > 0 && end_idx > 0 && start_idx < end_idx && end_idx <= length(signal)) {
          polya_signal <- signal[start_idx:end_idx]
        } else {
          message(sprintf("Invalid indices for read %s: start=%d, end=%d, signal_length=%d", read_id, start_idx, end_idx, length(signal)))
          polya_signal <- numeric(0)
        }
        list(
          read_id = read_id,
          pod5_file = pod5_file,
          polya_signal = polya_signal,
          signal_length = length(polya_signal)
        )
      })
      reader$close()
      result_df <- tibble::tibble(
        read_id = vapply(result, function(x) x$read_id, character(1)),
        pod5_file = vapply(result, function(x) x$pod5_file, character(1)),
        polya_signal = lapply(result, function(x) x$polya_signal),
        signal_length = vapply(result, function(x) x$signal_length, numeric(1))
      )
      dplyr::left_join(current_data, result_df, by = c("read_id", "pod5_file"))
    }, error = function(e) {
      warning(sprintf("Error processing %s: %s", pod5_path, e$message))
      NULL
    })
  })

  results <- dplyr::bind_rows(results)
  return(results)
}



