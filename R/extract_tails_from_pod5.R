#' Extract poly(A) tail signal segments from POD5 files
#'
#' Uses reticulate to import the Python library `pod5`. For each pod5 file referenced in
#' `polya_data`, opens a `pod5.Reader`, matches `read_id`, and extracts the numeric
#' slice `signal[poly_tail_start:poly_tail_end]` if indices are valid. Applies winsorization and interpolation.
#'
#' Only keeps reads if:
#' \itemize{
#'   \item `poly_tail_length > 10` (if the column exists), and
#'   \item `poly_tail_start > 0` (to exclude artifacts, e.g. caused by clogged pores).
#' }
#'
#' @param polya_data Data frame with columns:
#'   \describe{
#'     \item{read_id}{Character: ONT read identifiers}
#'     \item{filename}{Character: corresponding POD5 filename}
#'     \item{poly_tail_start}{Integer/numeric: start index of the poly(A) region}
#'     \item{poly_tail_end}{Integer/numeric: end index of the poly(A) region}
#'   }
#' @param pod5_dir Character. Path to the directory containing POD5 files.
#'
#' @return List of numeric vectors, each containing the extracted signal for a read.
#' @export
#'
#' @examples
#' \dontrun{
#' res <- extract_tails_from_pod5(polya_data, "/path/to/pod5_dir")
#' }
extract_tails_from_pod5 <- function(polya_data, pod5_dir) {
  i <- reader <- NULL

  # Check for required columns
  required_cols <- c("read_id", "filename", "poly_tail_start", "poly_tail_end")
  missing_cols <- setdiff(required_cols, colnames(polya_data))
  if (length(missing_cols) > 0) {
    stop(sprintf("polya_data is missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  if (!reticulate::py_module_available("pod5")) {
    stop("Python module 'pod5' is not available in the current environment.")
  }
  pod5 <- reticulate::import("pod5")

  # Keep only valid reads:
  # - poly_tail_length > 10 (if column exists)
  # - poly_tail_start > 0
  if ("poly_tail_length" %in% colnames(polya_data)) {
    polya_data <- polya_data[polya_data$poly_tail_length > 10 & polya_data$poly_tail_start > 0, ]
    if (nrow(polya_data) == 0) {
      stop("No valid reads found: require poly_tail_length > 10 and poly_tail_start > 0.")
    }
  } else {
    polya_data <- polya_data[polya_data$poly_tail_start > 0, ]
    if (nrow(polya_data) == 0) {
      stop("No valid reads found: require poly_tail_start > 0.")
    }
  }

  polya_by_file <- split(polya_data, polya_data$filename)

  signals_list <- list()
  for (filename in names(polya_by_file)) {
    current_data <- polya_by_file[[filename]]
    pod5_path <- file.path(pod5_dir, filename)
    if (!file.exists(pod5_path)) {
      warning(sprintf("POD5 file not found: %s", pod5_path))
      next
    }
    tryCatch({
      reader <- pod5$Reader(pod5_path)
      all_reads <- reticulate::iterate(reader$reads())
      available_ids <- sapply(all_reads, function(r) trimws(r$read_id))
      reads_by_id <- stats::setNames(all_reads, available_ids)
      for (i in seq_len(nrow(current_data))) {
        read_id <- trimws(current_data$read_id[i])
        if (!read_id %in% names(reads_by_id)) next
        signal <- reticulate::py_to_r(reads_by_id[[read_id]]$signal)
        start_idx <- current_data$poly_tail_start[i]
        end_idx <- current_data$poly_tail_end[i]
        if (is.numeric(start_idx) && is.numeric(end_idx) &&
            start_idx > 0 && end_idx > 0 &&
            start_idx < end_idx && end_idx <= length(signal)) {
          polya_signal <- signal[start_idx:end_idx]
          # Winsorize and interpolate
          polya_signal <- ninetails::winsorize_signal(polya_signal)
          polya_signal <- round(stats::approx(polya_signal, method = "linear", n = ceiling(0.2 * length(polya_signal)))[[2]], digits = 0)
        } else {
          polya_signal <- numeric(0)
        }
        signals_list[[read_id]] <- as.numeric(polya_signal)
      }
      reader$close()
    }, error = function(e) {
      warning(sprintf("Error processing %s: %s", pod5_path, e$message))
    })
  }
  return(signals_list)
}
