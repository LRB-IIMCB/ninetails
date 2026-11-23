################################################################################
# cDNA-SPECIFIC SIGNAL PROCESSING AND CLASSIFICATION (LACHESIS / DORADO)
################################################################################
# cDNA tails are classified with two dedicated models, one per read
# orientation (poly(A) reads and poly(T) reads).
#
# In contrast to the DRS pipelines, the raw (5 kHz) signal is used without
# downsampling and without pseudomove prefiltering: 100-sample windows slide along
# the entire tail and every window is classified. Class labels always refer to
# the residue in the original RNA tail (A, C, G, U), in both orientations, so no
# reverse complementing of the predictions is needed.
################################################################################

# Internal: frozen settings of the cDNA models and of the calling rule.
# Thresholds were tuned on validation reads and must not be changed without
# retraining/re-evaluating the models.
.get_cdna_model_config <- function(tail_type) {

  if (tail_type == "polyA") {
    model_file <- "gasf_gadf_cdna_polya_model_20261122.h5"
    theta_mid <- 0.95
    theta_edge <- 0.98
  } else if (tail_type == "polyT") {
    model_file <- "gasf_gadf_cdna_polyt_model_20261122.h5"
    theta_mid <- 0.80
    theta_edge <- 0.70
  } else {
    stop("Wrong tail_type provided. Please provide either 'polyA' or 'polyT'",
         call. = FALSE)
  }

  config <- list(model_file = model_file,
                 theta_mid = theta_mid,   # threshold for windows inside the tail
                 theta_edge = theta_edge, # threshold for windows near a tail boundary
                 min_run = 2,             # min. number of consecutive windows above threshold
                 flank = 100,             # samples extracted on each side of the tail
                 batch_size = 256)        # windows per prediction batch (memory bound)

  return(config)
} # .get_cdna_model_config


################################################################################
# TAIL COORDINATES (LACHESIS / DORADO ADAPTERS)
################################################################################
#' Read poly(A)/poly(T) tail coordinates for the cDNA pipeline
#'
#' Converts the tail coordinates produced by Lachesis or Dorado into one
#' common table used by the cDNA pipeline, so that the downstream functions
#' do not depend on the tool that delimited the tail.
#'
#' @details
#' Coordinates are positions in the raw signal (samples), with the tail
#' spanning \code{signal[(tail_start + 1):tail_end]}.
#'
#' @section Lachesis input:
#' The tail is delimited by \code{adapter_end} and \code{polya_end} in both
#' orientations. The orientation is taken from the \code{strand} column:
#' \code{"+"} is a poly(A) read, \code{"-"} is a poly(T) read.
#'
#' @section Dorado input:
#' The tail is delimited by \code{poly_tail_start} and \code{poly_tail_end}.
#' The summary does not carry the read orientation, so \code{tail_type} is
#' \code{NA} and has to be filled in by the caller (see
#' \code{\link{detect_orientation_multiple}}). Both \code{alignment_mapq} and
#' \code{alignment_mapping_quality}, as well as \code{filename} and
#' \code{input_filename}, are accepted.
#'
#' @param tail_coordinates Character string or data frame. Path to the
#'   Lachesis output / Dorado summary file, or the same table in memory.
#'
#' @param source Character string ["lachesis"]. Tool that produced the
#'   coordinates: \code{"lachesis"} or \code{"dorado"}.
#'
#' @return A data frame with one row per read:
#' \describe{
#'   \item{read_id}{Character. Read identifier}
#'   \item{filename}{Character. POD5 file name (\code{NA} if not available)}
#'   \item{contig}{Character. Reference name (\code{NA} if not available)}
#'   \item{tail_start}{Numeric. First sample of the tail (0-based)}
#'   \item{tail_end}{Numeric. End of the tail (exclusive)}
#'   \item{polya_length}{Numeric. Estimated tail length in nucleotides}
#'   \item{qc_tag}{Character. Lachesis QC tag or Dorado mapping quality}
#'   \item{tail_type}{Character. \code{"polyA"}, \code{"polyT"},
#'     \code{"unidentified"} or \code{NA}}
#' }
#'
#' @seealso \code{\link{extract_tail_signals_cdna}},
#'   \code{\link{check_tails_dorado_cDNA}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' coordinates <- ninetails::read_tail_coordinates_cdna(
#'   tail_coordinates = '/path/to/lachesis_polya.tsv',
#'   source = "lachesis")
#'
#' }
read_tail_coordinates_cdna <- function(tail_coordinates,
                                       source = "lachesis") {

  #assertions
  if (missing(tail_coordinates)) {
    stop("Tail coordinates are missing. Please provide a valid tail_coordinates argument.",
         call. = FALSE)
  }

  if (is_string(tail_coordinates)) {
    assert_file_exists(tail_coordinates, "Tail coordinates")
    tail_coordinates <- vroom::vroom(tail_coordinates, show_col_types = FALSE)
  }

  if (!is.data.frame(tail_coordinates) || nrow(tail_coordinates) == 0) {
    stop("Empty data frame provided as an input (tail_coordinates). Please provide valid input",
         call. = FALSE)
  }

  #coerce tibble to df
  tail_coordinates <- as.data.frame(tail_coordinates)

  if (source == "lachesis") {
    required_cols <- c("readname", "contig", "adapter_end", "polya_end",
                       "polya_length", "qc_tag", "strand")
  } else if (source == "dorado") {
    # From dorado 1.4.0, column name for pod5 files is "input_filename"
    if (!"filename" %in% colnames(tail_coordinates) && "input_filename" %in% colnames(tail_coordinates)) {
      colnames(tail_coordinates)[colnames(tail_coordinates) == "input_filename"] <- "filename"
    }
    # newer dorado summaries name the mapping quality "alignment_mapping_quality"
    if (!"alignment_mapq" %in% colnames(tail_coordinates) && "alignment_mapping_quality" %in% colnames(tail_coordinates)) {
      colnames(tail_coordinates)[colnames(tail_coordinates) == "alignment_mapping_quality"] <- "alignment_mapq"
    }
    required_cols <- c("read_id", "poly_tail_length", "poly_tail_start", "poly_tail_end")
  } else {
    stop("Wrong source provided. Please provide either 'lachesis' or 'dorado'",
         call. = FALSE)
  }

  missing_cols <- setdiff(required_cols, colnames(tail_coordinates))
  if (length(missing_cols) > 0) {
    stop(sprintf("Required columns missing: %s", paste(missing_cols, collapse = ", ")),
         call. = FALSE)
  }

  # columns absent from the input are reported as NA
  .optional_column <- function(column_name) {
    if (column_name %in% colnames(tail_coordinates)) {
      as.character(tail_coordinates[[column_name]])
    } else {
      rep(NA_character_, nrow(tail_coordinates))
    }
  } # .optional_column

  if (source == "lachesis") {
    # "+" strand: poly(A) in the signal; "-" strand: poly(T) in the signal
    tail_type <- ifelse(tail_coordinates$strand == "+", "polyA",
                        ifelse(tail_coordinates$strand == "-", "polyT", "unidentified"))
    tail_type[is.na(tail_type)] <- "unidentified"

    result <- data.frame(read_id = as.character(tail_coordinates$readname),
                         filename = .optional_column("filename"),
                         contig = as.character(tail_coordinates$contig),
                         tail_start = as.numeric(tail_coordinates$adapter_end),
                         tail_end = as.numeric(tail_coordinates$polya_end),
                         polya_length = as.numeric(tail_coordinates$polya_length),
                         qc_tag = as.character(tail_coordinates$qc_tag),
                         tail_type = tail_type,
                         stringsAsFactors = FALSE)
  } else {
    result <- data.frame(read_id = as.character(tail_coordinates$read_id),
                         filename = .optional_column("filename"),
                         contig = .optional_column("alignment_genome"),
                         tail_start = as.numeric(tail_coordinates$poly_tail_start),
                         tail_end = as.numeric(tail_coordinates$poly_tail_end),
                         polya_length = as.numeric(tail_coordinates$poly_tail_length),
                         qc_tag = .optional_column("alignment_mapq"),
                         tail_type = NA_character_,
                         stringsAsFactors = FALSE)
  }

  # one row per read
  result <- result[!duplicated(result$read_id), ]

  return(result)
}



################################################################################
# RAW TAIL SIGNAL EXTRACTION
################################################################################
#' Extract raw poly(A)/poly(T) tail signals with flanks from POD5 files (cDNA)
#'
#' Extracts the raw signal of each tail together with 100 samples on each
#' side (clipped at the read ends). In contrast to
#' \code{\link{extract_tails_from_pod5}}, the signal is neither winsorized
#' nor downsampled, because the cDNA models work on the raw 5 kHz signal.
#'
#' @details
#' If the \code{filename} column is not available (Lachesis output), every
#' POD5 file in \code{pod5_dir} is scanned for the requested reads.
#'
#' @param tail_coordinates Data frame. Output of
#'   \code{\link{read_tail_coordinates_cdna}}.
#'
#' @param pod5_dir Character string. Path to the directory with POD5 files.
#'
#' @param num_cores Numeric [1]. Number of physical cores to use.
#'
#' @return A named list of numeric vectors (raw signal), one per read with
#'   valid coordinates, named with read IDs. The tail starts at sample
#'   \code{min(tail_start, 100)} (0-based) of each vector.
#'
#' @seealso \code{\link{read_tail_coordinates_cdna}},
#'   \code{\link{predict_nonA_cdna}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' signal_list <- ninetails::extract_tail_signals_cdna(
#'   tail_coordinates = coordinates,
#'   pod5_dir = '/path/to/pod5/',
#'   num_cores = 2)
#'
#' }
extract_tail_signals_cdna <- function(tail_coordinates,
                                      pod5_dir,
                                      num_cores = 1) {

  #assertions
  if (missing(tail_coordinates)) {
    stop("Tail coordinates are missing. Please provide a valid tail_coordinates argument.",
         call. = FALSE)
  }

  if (missing(pod5_dir)) {
    stop("Directory with POD5 files is missing. Please provide a valid pod5_dir argument.",
         call. = FALSE)
  }

  if (!is.data.frame(tail_coordinates) || nrow(tail_coordinates) == 0) {
    stop("Empty data frame provided as an input (tail_coordinates). Please provide valid input",
         call. = FALSE)
  }

  assert_condition(is.numeric(num_cores),
                   "Declared core number must be numeric. Please provide a valid argument.")
  assert_dir_exists(pod5_dir, "POD5")

  required_cols <- c("read_id", "tail_start", "tail_end")
  missing_cols <- setdiff(required_cols, colnames(tail_coordinates))
  if (length(missing_cols) > 0) {
    stop(sprintf("Required columns missing: %s", paste(missing_cols, collapse = ", ")),
         call. = FALSE)
  }

  # Filter valid reads
  valid <- !is.na(tail_coordinates$tail_start) & !is.na(tail_coordinates$tail_end) &
    tail_coordinates$tail_start > 0 & tail_coordinates$tail_end > tail_coordinates$tail_start
  tail_coordinates <- tail_coordinates[valid, ]

  if (nrow(tail_coordinates) == 0) {
    return(list()) # No valid reads, return empty list
  }

  python_script <- system.file("extdata", "extract_pod5_signals.py", package = "ninetails")
  if (!file.exists(python_script)) {
    stop("POD5 extraction script not found. Please reinstall ninetails package.",
         call. = FALSE)
  }

  # the python helper expects dorado-style column names
  python_input <- data.frame(read_id = tail_coordinates$read_id,
                             poly_tail_start = tail_coordinates$tail_start,
                             poly_tail_end = tail_coordinates$tail_end,
                             stringsAsFactors = FALSE)
  if ("filename" %in% colnames(tail_coordinates) && !any(is.na(tail_coordinates$filename))) {
    python_input$filename <- tail_coordinates$filename
  }

  input_file <- tempfile(pattern = "cdna_input_", fileext = ".csv")
  output_file <- tempfile(pattern = "cdna_output_", fileext = ".pkl")
  on.exit(unlink(c(input_file, output_file)), add = TRUE)

  utils::write.csv(python_input, input_file, row.names = FALSE)

  python_cmd <- sprintf('"%s" "%s" --input "%s" --pod5_dir "%s" --output "%s" --num_cores %d --raw --flank %d',
                        reticulate::py_config()$python,
                        python_script,
                        input_file,
                        pod5_dir,
                        output_file,
                        as.integer(num_cores),
                        .get_cdna_model_config("polyA")$flank)

  exit_code <- system(python_cmd, intern = FALSE)

  if (exit_code != 0 || !file.exists(output_file)) {
    stop("POD5 extraction failed. Ensure 'pod5' Python module is installed: pip install pod5",
         call. = FALSE)
  }

  pkl <- reticulate::import("pickle", convert = FALSE)
  py_builtin <- reticulate::import_builtins()
  f <- py_builtin$open(output_file, "rb")
  signal_list <- reticulate::py_to_r(pkl$load(f))
  f$close()

  # Ensure numeric vectors
  signal_list <- lapply(signal_list, function(x) as.numeric(unlist(x)))

  return(signal_list)
}



################################################################################
# NORMALISATION, WINDOWS AND GRAMIAN ANGULAR FIELDS
################################################################################
#' Normalise a raw cDNA tail signal
#'
#' Centres the signal on the tail median, scales it by the sample-to-sample
#' noise of the tail and squashes it to the (-1, 1) range:
#' \code{tanh(((signal - median(tail)) / noise) / 8)}, where
#' \code{noise = 1.4826 * mad(diff(tail)) / sqrt(2)}.
#'
#' @details
#' The statistics come from the tail only, but the whole provided signal
#' (tail with flanks) is transformed. The result does not depend on the
#' signal units (raw ADC values or picoamperes give the same output).
#'
#' @param signal Numeric vector. Raw signal (tail with flanks).
#'
#' @param tail_start Numeric. First sample of the tail within \code{signal}
#'   (0-based).
#'
#' @param tail_end Numeric. End of the tail within \code{signal} (exclusive).
#'
#' @return Numeric vector of the same length as \code{signal}.
#'
#' @seealso \code{\link{create_tail_windows_cdna}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' normalised <- ninetails::normalise_tail_cdna(signal = signal,
#'                                              tail_start = 100,
#'                                              tail_end = 800)
#'
#' }
normalise_tail_cdna <- function(signal, tail_start, tail_end) {

  #assertions
  if (missing(signal)) {
    stop("Signal is missing. Please provide a valid signal argument.", call. = FALSE)
  }
  if (missing(tail_start)) {
    stop("Tail start is missing. Please provide a valid tail_start argument.", call. = FALSE)
  }
  if (missing(tail_end)) {
    stop("Tail end is missing. Please provide a valid tail_end argument.", call. = FALSE)
  }

  assert_condition(is.numeric(signal), "Signal must be numeric. Please provide a valid argument.")
  assert_condition(tail_start >= 0 && tail_end > tail_start && tail_end <= length(signal),
                   "Tail coordinates must lie within the signal. Please provide valid tail_start and tail_end arguments.")

  tail_signal <- signal[(tail_start + 1):tail_end]
  tail_median <- stats::median(tail_signal)

  # robust noise estimate from first differences (insensitive to level shifts)
  tail_diff <- diff(tail_signal)
  noise <- 1.4826 * stats::median(abs(tail_diff - stats::median(tail_diff))) / sqrt(2)

  #prevent bugs: flat or single-sample tails
  if (is.na(noise) || noise <= 0) {
    noise <- 1
  }

  normalised_signal <- tanh(((signal - tail_median) / noise) / 8)

  return(normalised_signal)
}



#' Split a normalized cDNA tail signal into sliding windows
#'
#' Cuts 100-sample windows every 25 samples, with window centers running
#' from 50 samples before the tail to 50 samples behind it. Windows that
#' would leave the provided signal are skipped (no padding).
#'
#' @inheritParams normalise_tail_cdna
#'
#' @return A list with two elements:
#' \describe{
#'   \item{centres}{Numeric vector. Window centres within \code{signal}
#'     (0-based)}
#'   \item{windows}{Numeric matrix with one normalised window per row and
#'     100 columns}
#' }
#'
#' @seealso \code{\link{normalise_tail_cdna}}, \code{\link{create_gaf_cdna}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' tail_windows <- ninetails::create_tail_windows_cdna(signal = signal,
#'                                                     tail_start = 100,
#'                                                     tail_end = 800)
#'
#' }
create_tail_windows_cdna <- function(signal, tail_start, tail_end) {

  normalised_signal <- ninetails::normalise_tail_cdna(signal, tail_start, tail_end)

  centres <- seq(tail_start - 50, tail_end + 50, by = 25)
  centres <- centres[centres - 50 >= 0 & centres + 50 <= length(normalised_signal)]

  if (length(centres) == 0) {
    return(list(centres = numeric(0), windows = matrix(numeric(0), nrow = 0, ncol = 100))) # Signal too short, return no windows
  }

  # row i holds samples (centre - 50 + 1):(centre + 50) in R (1-based) indexing
  window_index <- outer(centres - 50, seq_len(100), "+")
  windows <- matrix(normalised_signal[window_index], nrow = length(centres))

  result <- list(centres = centres, windows = windows)

  return(result)
}



#' Create gramian angular fields from cDNA signal windows
#'
#' Converts normalised signal windows into two-channel images: the gramian
#' angular summation field \code{cos(phi_i + phi_j)} with
#' \code{phi = acos(x)} and the field \code{sin(psi_i + psi_j)} with
#' \code{psi = asin(x)}.
#'
#' @details
#' In contrast to \code{\link{create_gaf}}, the windows are \strong{not}
#' rescaled to their own minimum and maximum: the normalisation of
#' \code{\link{normalise_tail_cdna}} already maps the signal to (-1, 1) and
#' keeps the amplitude relative to the whole tail, which carries most of the
#' information in cDNA signals.
#'
#' @param windows Numeric matrix. One normalised window per row
#'   (as returned by \code{\link{create_tail_windows_cdna}}).
#'
#' @return A numeric array with dimensions (number of windows, window
#'   length, window length, 2).
#'
#' @seealso \code{\link{create_tail_windows_cdna}},
#'   \code{\link{predict_nonA_cdna}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' gafs <- ninetails::create_gaf_cdna(windows = tail_windows$windows)
#'
#' }
create_gaf_cdna <- function(windows) {

  #assertions
  if (missing(windows)) {
    stop("Windows are missing. Please provide a valid windows argument.", call. = FALSE)
  }

  assert_condition(is.matrix(windows) && is.numeric(windows),
                   "Windows must be a numeric matrix. Please provide a valid argument.")

  n_windows <- nrow(windows)
  n_points <- ncol(windows)

  x <- pmin(pmax(windows, -1), 1)
  s <- sqrt(1 - x^2)

  # [b, i, j] = value at position i of window b; the transposed copy gives position j
  x_i <- array(x, c(n_windows, n_points, n_points))
  x_j <- aperm(x_i, c(1, 3, 2))
  s_i <- array(s, c(n_windows, n_points, n_points))
  s_j <- aperm(s_i, c(1, 3, 2))

  # cos(a + b) = cos(a)cos(b) - sin(a)sin(b); sin(a + b) with psi = asin(x) analogously
  gafs <- array(c(x_i * x_j - s_i * s_j,
                  x_i * s_j + x_j * s_i),
                c(n_windows, n_points, n_points, 2))

  return(gafs)
}



################################################################################
# CALLING NON-ADENOSINES
################################################################################
#' Call non-adenosine residues from window-wise cDNA predictions
#'
#' Turns the class probabilities of consecutive signal windows into
#' non-adenosine calls. Windows whose non-A score (\code{1 - p(A)}) exceeds
#' the threshold are merged into runs; each sufficiently long run gives one
#' call.
#'
#' @details
#' Windows centred within 100 samples of a tail boundary are compared with
#' \code{theta_edge}, the remaining ones with \code{theta_mid}. The call is
#' positioned at the centre of the highest-scoring window of the run and its
#' class is the one with the highest summed probability over the run. Calls
#' positioned outside the tail are dropped.
#'
#' @param centres Numeric vector. Window centres (0-based, same coordinate
#'   system as \code{tail_start} and \code{tail_end}).
#'
#' @param probabilities Numeric matrix. Class probabilities with one row
#'   per window and four columns (A, C, G, U).
#'
#' @param tail_start Numeric. First sample of the tail (0-based).
#'
#' @param tail_end Numeric. End of the tail (exclusive).
#'
#' @param theta_mid Numeric. Threshold for windows inside the tail.
#'
#' @param theta_edge Numeric. Threshold for windows near a tail boundary.
#'
#' @param min_run Numeric [2]. Minimal number of consecutive windows above
#'   the threshold required for a call.
#'
#' @return A data frame with one row per call:
#' \describe{
#'   \item{prediction}{Character. Residue in the RNA tail (C, G or U)}
#'   \item{centr_signal_pos}{Numeric. Position of the call in the signal}
#'   \item{score}{Numeric. Highest non-A score of the run}
#' }
#'
#' @seealso \code{\link{predict_nonA_cdna}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' calls <- ninetails::call_nonA_cdna(centres = tail_windows$centres,
#'                                    probabilities = probabilities,
#'                                    tail_start = 100,
#'                                    tail_end = 800,
#'                                    theta_mid = 0.95,
#'                                    theta_edge = 0.98)
#'
#' }
call_nonA_cdna <- function(centres,
                           probabilities,
                           tail_start,
                           tail_end,
                           theta_mid,
                           theta_edge,
                           min_run = 2) {

  #assertions
  if (missing(centres)) {
    stop("Centres are missing. Please provide a valid centres argument.", call. = FALSE)
  }
  if (missing(probabilities)) {
    stop("Probabilities are missing. Please provide a valid probabilities argument.", call. = FALSE)
  }
  if (missing(theta_mid) || missing(theta_edge)) {
    stop("Thresholds are missing. Please provide valid theta_mid and theta_edge arguments.", call. = FALSE)
  }

  assert_condition(is.matrix(probabilities) && ncol(probabilities) == 4,
                   "Probabilities must be a matrix with 4 columns (A, C, G, U). Please provide a valid argument.")
  assert_condition(length(centres) == nrow(probabilities),
                   "Number of centres must match the number of rows in probabilities.")

  prediction_dict <- c("C", "G", "U")

  calls <- data.frame(prediction = character(0),
                      centr_signal_pos = numeric(0),
                      score = numeric(0),
                      stringsAsFactors = FALSE)

  if (length(centres) == 0) {
    return(calls) # No windows, return no calls
  }

  nonA_score <- 1 - probabilities[, 1]
  near_boundary <- abs(centres - tail_start) <= 100 | abs(centres - tail_end) <= 100
  above_threshold <- nonA_score > ifelse(near_boundary, theta_edge, theta_mid)

  runs <- rle(above_threshold)
  run_ends <- cumsum(runs$lengths)
  run_starts <- run_ends - runs$lengths + 1

  for (i in which(runs$values & runs$lengths >= min_run)) {
    run_index <- run_starts[i]:run_ends[i]
    best_window <- run_index[which.max(nonA_score[run_index])]
    position <- centres[best_window]

    # calls positioned outside the tail are dropped
    if (position < tail_start || position > tail_end) {
      next
    }

    class_index <- which.max(colSums(probabilities[run_index, 2:4, drop = FALSE]))

    calls[nrow(calls) + 1, ] <- list(prediction_dict[class_index],
                                     position,
                                     nonA_score[best_window])
  } # for run

  return(calls)
}



#' Predict non-adenosine residues in cDNA tails
#'
#' Classifies sliding windows of raw tail signals with the orientation-
#' specific cDNA model and calls non-adenosine residues.
#'
#' @details
#' For every read the tail signal is normalised
#' (\code{\link{normalise_tail_cdna}}), cut into windows
#' (\code{\link{create_tail_windows_cdna}}), converted to gramian angular
#' fields (\code{\link{create_gaf_cdna}}) and classified. Calls are made
#' with \code{\link{call_nonA_cdna}} using thresholds frozen together with
#' the models.
#'
#' The position of each call is reported as the distance from the 3' end of
#' the RNA tail in nucleotides. Poly(A) reads end with the 3' end of the
#' tail, poly(T) reads begin with it:
#' \itemize{
#'   \item poly(A): \code{(tail_end - position) / samples_per_nt}
#'   \item poly(T): \code{(position - tail_start) / samples_per_nt}
#' }
#' where \code{samples_per_nt = (tail_end - tail_start) / polya_length}.
#'
#' @param signal_list Named list of numeric vectors. Raw tail signals with
#'   flanks, as returned by \code{\link{extract_tail_signals_cdna}}.
#'
#' @param tail_coordinates Data frame. Output of
#'   \code{\link{read_tail_coordinates_cdna}} (columns \code{read_id},
#'   \code{tail_start}, \code{tail_end} and \code{polya_length} are used).
#'
#' @param tail_type Character string. Orientation of the reads:
#'   \code{"polyA"} or \code{"polyT"}.
#'
#' @param keras_model_path Character string or missing. Full path of the
#'   \code{.h5} file with a Keras model. If missing (default), the
#'   pretrained cDNA model for the given \code{tail_type} is loaded.
#'
#' @return A data frame with one row per called non-adenosine residue:
#' \describe{
#'   \item{read_id}{Character. Read identifier}
#'   \item{prediction}{Character. Residue in the RNA tail (C, G or U)}
#'   \item{est_nonA_pos}{Numeric. Estimated position from the 3' end of the
#'     tail (nucleotides)}
#'   \item{score}{Numeric. Highest non-A score of the call}
#' }
#'
#' @seealso \code{\link{create_outputs_dorado_cdna}},
#'   \code{\link{load_keras_model}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' nonA_calls <- ninetails::predict_nonA_cdna(
#'   signal_list = signal_list,
#'   tail_coordinates = coordinates,
#'   tail_type = "polyA")
#'
#' }
predict_nonA_cdna <- function(signal_list,
                              tail_coordinates,
                              tail_type,
                              keras_model_path) {

  #assertions
  if (missing(signal_list)) {
    stop("Signal list is missing. Please provide a valid signal_list argument.", call. = FALSE)
  }
  if (missing(tail_coordinates)) {
    stop("Tail coordinates are missing. Please provide a valid tail_coordinates argument.", call. = FALSE)
  }
  if (missing(tail_type)) {
    stop("Tail type is missing. Please provide a valid tail_type argument.", call. = FALSE)
  }

  assert_condition(is.list(signal_list), "Signal list must be a list. Please provide a valid argument.")
  assert_condition(is.data.frame(tail_coordinates),
                   "Tail coordinates must be a data frame. Please provide a valid argument.")

  required_cols <- c("read_id", "tail_start", "tail_end", "polya_length")
  missing_cols <- setdiff(required_cols, colnames(tail_coordinates))
  if (length(missing_cols) > 0) {
    stop(sprintf("Required columns missing: %s", paste(missing_cols, collapse = ", ")),
         call. = FALSE)
  }

  config <- .get_cdna_model_config(tail_type)

  nonA_calls <- data.frame(read_id = character(0),
                           prediction = character(0),
                           est_nonA_pos = numeric(0),
                           score = numeric(0),
                           stringsAsFactors = FALSE)

  tail_coordinates <- tail_coordinates[tail_coordinates$read_id %in% names(signal_list), ]

  if (nrow(tail_coordinates) == 0) {
    return(nonA_calls) # No reads with signal, return no calls
  }

  if (rlang::is_missing(keras_model_path)) {
    keras_model_path <- system.file("extdata", "cnn_model", config$model_file, package = "ninetails")
  }
  keras_model <- ninetails::load_keras_model(keras_model_path)

  # Header comm
  cat(paste0('[', as.character(Sys.time()), '] ', 'Classifying ', tail_type, ' tail signals...', '\n', sep = ''))

  # reads are processed in small groups to bound the memory used by the images
  read_groups <- split(seq_len(nrow(tail_coordinates)), ceiling(seq_len(nrow(tail_coordinates)) / 25))

  pb <- utils::txtProgressBar(min = 0, max = length(read_groups), style = 3, width = 50, char = "=", file = stderr())

  calls_list = list()

  for (g in seq_along(read_groups)) {
    group_coordinates <- tail_coordinates[read_groups[[g]], ]

    group_windows <- lapply(seq_len(nrow(group_coordinates)), function(i) {
      signal <- signal_list[[group_coordinates$read_id[i]]]
      # position of the tail within the extracted (flanked) signal
      local_start <- min(group_coordinates$tail_start[i], config$flank)
      local_end <- local_start + (group_coordinates$tail_end[i] - group_coordinates$tail_start[i])
      if (length(signal) < local_end) {
        return(NULL) # Truncated signal, skip read
      }
      tail_windows <- ninetails::create_tail_windows_cdna(signal, local_start, local_end)
      tail_windows[["local_start"]] <- local_start
      tail_windows[["local_end"]] <- local_end
      tail_windows
    })

    all_windows <- do.call(rbind, lapply(group_windows, function(x) x$windows))

    if (is.null(all_windows) || nrow(all_windows) == 0) {
      utils::setTxtProgressBar(pb, g)
      next
    }

    probabilities <- matrix(NA_real_, nrow = nrow(all_windows), ncol = 4)
    batches <- split(seq_len(nrow(all_windows)), ceiling(seq_len(nrow(all_windows)) / config$batch_size))
    for (batch in batches) {
      probabilities[batch, ] <- keras_model %>%
        stats::predict(ninetails::create_gaf_cdna(all_windows[batch, , drop = FALSE]), verbose = 0)
    }

    first_row <- 0
    for (i in seq_len(nrow(group_coordinates))) {
      tail_windows <- group_windows[[i]]
      n_windows <- length(tail_windows$centres)
      if (is.null(tail_windows) || n_windows == 0) {
        next
      }
      read_rows <- first_row + seq_len(n_windows)
      first_row <- first_row + n_windows

      calls <- ninetails::call_nonA_cdna(centres = tail_windows$centres,
                                         probabilities = probabilities[read_rows, , drop = FALSE],
                                         tail_start = tail_windows$local_start,
                                         tail_end = tail_windows$local_end,
                                         theta_mid = config$theta_mid,
                                         theta_edge = config$theta_edge,
                                         min_run = config$min_run)

      if (nrow(calls) > 0) {
        samples_per_nt <- (tail_windows$local_end - tail_windows$local_start) / group_coordinates$polya_length[i]
        # distance from the 3' end of the RNA tail:
        # poly(A) reads end with the 3' end, poly(T) reads begin with it
        if (tail_type == "polyA") {
          est_nonA_pos <- (tail_windows$local_end - calls$centr_signal_pos) / samples_per_nt
        } else {
          est_nonA_pos <- (calls$centr_signal_pos - tail_windows$local_start) / samples_per_nt
        }
        calls_list[[length(calls_list) + 1]] <- data.frame(read_id = group_coordinates$read_id[i],
                                                           prediction = calls$prediction,
                                                           est_nonA_pos = round(est_nonA_pos, 2),
                                                           score = calls$score,
                                                           stringsAsFactors = FALSE)
      }
    } # for read

    utils::setTxtProgressBar(pb, g)
  } # for read group

  close(pb)

  if (length(calls_list) > 0) {
    nonA_calls <- do.call(rbind, calls_list)
  }

  # Done comm
  cat(paste0('[', as.character(Sys.time()), '] ', 'Done!', '\n', sep = ''))

  return(nonA_calls)
}
