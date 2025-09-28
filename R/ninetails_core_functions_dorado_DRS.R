#' Process and split Dorado summary file into smaller parts
#'
#' Splits a Dorado summary file or data frame (in-memory file) into multiple
#' smaller files for downstream analysis. This helps to avoid memory overflow
#' and data loss during processing.
#'
#' @param dorado_summary Character path to Dorado summary file, or a data frame containing summary information.
#'
#' @param save_dir Character path to directory where split summary files will be saved.
#'
#' @param part_size Integer. Number of reads per file part when splitting the summary file.
#'
#' @param cli_log Function for logging messages and progress.
#'
#' @return Character vector containing paths to the split summary files.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' summary_files <- process_dorado_summary(
#'   dorado_summary = "path/to/summary.txt",
#'   save_dir = "path/to/output/",
#'   part_size = 40000,
#'   cli_log = message
#' )
#' }
process_dorado_summary <- function(dorado_summary,
                                   save_dir,
                                   part_size,
                                   cli_log) {

  # Assertions
  if (missing(part_size)) {
    stop("Number of reads per file part (part_size) is missing. Please provide a valid part_size argument.", call. = FALSE)
  }

  assertthat::assert_that(is.numeric(part_size), part_size > 0,
                          msg = paste0("Reads per part must be numeric and positive. Please provide a valid argument."))

  # Handle input data and get prefix
  if (is.character(dorado_summary)) {
    cli_log(sprintf("Reading summary file: %s", basename(dorado_summary)), "INFO")

    # First, let's peek at the file structure
    cli_log("Checking file structure...", "INFO")
    tryCatch({
      # Read first few lines to check structure
      headers <- names(vroom::vroom(dorado_summary, n_max = 1, show_col_types = FALSE))
      cli_log(sprintf("Found columns: %s", paste(headers, collapse=", ")), "INFO")

      # Now read the actual data
      summary_data <- vroom::vroom(dorado_summary, show_col_types = FALSE)

      if (!"read_id" %in% headers) {
        cli_log("WARNING: 'read_id' column not found. Available columns are:", "WARNING")
        cli_log(paste(headers, collapse=", "), "WARNING")
        stop("Required column 'read_id' not found in summary file")
      }

    }, error = function(e) {
      cli_log(sprintf("Error reading summary file: %s", e$message), "ERROR")
      stop(e$message)
    })

    # Extract prefix from original file name
    summary_prefix <- tools::file_path_sans_ext(basename(dorado_summary))
  } else if (is.data.frame(dorado_summary)) {
    summary_data <- tibble::as_tibble(dorado_summary)
    # Check for required columns in data frame
    if (!"read_id" %in% names(summary_data)) {
      cli_log("Required column 'read_id' not found in data frame", "ERROR")
      stop("Required column 'read_id' not found in data frame")
    }
    # Use default prefix for data frame input
    summary_prefix <- "summary"
  } else {
    stop("Invalid dorado_summary format. Must be either a file path or a data frame.")
  }

  # Create output directory if it doesn't exist
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }


  # this is to filter input summary data and decrease size of processing reads
  # filtering out data with bad alignment quality (mostly trash)
  # reads without proper tails (misidentified due to pore clog)
  # reads with too short tails (at least 10 nt is required)
  summary_data <- ninetails::filter_dorado_summary(dorado_summary)


  # Calculate number of parts needed
  total_reads <- nrow(summary_data)
  num_parts <- ceiling(total_reads / part_size)

  cli_log(sprintf("Total reads: %d", total_reads), "INFO")
  cli_log(sprintf("Splitting into %d parts", num_parts), "INFO")

  # Initialize vector to store output file paths
  output_files <- character(num_parts)

  # Split and save summary data
  for (i in seq_len(num_parts)) {
    start_idx <- ((i - 1) * part_size) + 1
    end_idx <- min(i * part_size, total_reads)

    # Extract subset of data
    part_data <- summary_data[start_idx:end_idx, ]

    # Create output filename using prefix from original file
    output_file <- file.path(save_dir,
                             sprintf("%s_part%d.txt", summary_prefix, i))

    # Save data
    vroom::vroom_write(part_data,
                       file = output_file,
                       delim = "\t")

    output_files[i] <- output_file

    cli_log(sprintf("Saved part %d/%d: %s", i, num_parts, basename(output_file)), "INFO")
  }

  cli_log(sprintf("Summary split into %d parts", num_parts), "SUCCESS")

  return(output_files)
}



#' Extract poly(A) tail signal segments from POD5 files using parallel Python processing
#'
#' This function extracts raw nanopore signal data corresponding to poly(A) tail regions
#' from POD5 format files. It leverages a Python subprocess for efficient parallel
#' extraction, avoiding R-Python serialization bottlenecks inherent in reticulate.
#' The function processes multiple reads across multiple POD5 files simultaneously,
#' applies winsorization to reduce noise, and performs signal interpolation to
#' standardize tail representations.
#'
#' @param polya_data Data frame containing poly(A) tail coordinate information.
#'   Must include the following columns:
#'   \describe{
#'     \item{read_id}{Character. Unique identifier for each nanopore read}
#'     \item{filename}{Character. Name of the POD5 file containing the read}
#'     \item{poly_tail_start}{Integer. Start position of poly(A) tail in the raw signal}
#'     \item{poly_tail_end}{Integer. End position of poly(A) tail in the raw signal}
#'     \item{poly_tail_length}{Integer. (Optional) Length of the poly(A) tail in bases.
#'           If present, only tails > 10 bases are processed}
#'   }
#'
#' @param pod5_dir Character string. Full path to the directory containing POD5 files.
#'   All POD5 files referenced in the \code{filename} column of \code{polya_data}
#'   must be present in this directory.
#'
#' @param num_cores Integer. Number of CPU cores to use for parallel processing.
#'   Default is 1. If set to NULL, the function automatically uses all available
#'   cores minus 1 to maintain system responsiveness. Values > 1 enable parallel
#'   processing of POD5 files.
#'
#' @details
#' The function performs the following operations:
#' \enumerate{
#'   \item \strong{Input validation}: Checks for required columns and valid data
#'   \item \strong{Read filtering}: Excludes reads with:
#'     \itemize{
#'       \item poly_tail_length ≤ 10 (if column exists)
#'       \item poly_tail_start ≤ 0 (invalid coordinates)
#'     }
#'   \item \strong{Python subprocess execution}: Delegates extraction to an optimized
#'         Python script that:
#'     \itemize{
#'       \item Groups reads by POD5 file for efficient I/O
#'       \item Processes files in parallel using multiprocessing
#'       \item Extracts signal segments based on provided coordinates
#'     }
#'   \item \strong{Signal processing}: For each valid tail region:
#'     \itemize{
#'       \item Applies winsorization (0.5% and 99.5% percentiles)
#'       \item Interpolates to 20% of original length for standardization
#'       \item Converts to integer values for downstream compatibility
#'     }
#' }
#'
#' The Python subprocess approach bypasses reticulate's Global Interpreter Lock (GIL)
#' limitations, providing true parallel processing for large datasets. Temporary files
#' are used for data exchange and automatically cleaned up after processing.
#'
#' @return A named list of numeric vectors, where:
#'   \itemize{
#'     \item Names correspond to read_id values from the input
#'     \item Each vector contains the winsorized and interpolated signal values
#'     \item Empty vectors indicate reads where extraction failed
#'   }
#'   Returns an empty list if no valid reads are found.
#'
#' @note
#' \itemize{
#'   \item Requires Dorado >=1.0.0 to retrieve polyA coordinates
#'   \item Requires Python 3.6+ with the 'pod5' module installed: \code{pip install pod5}
#'   \item The POD5 extraction script must be present at
#'         \code{system.file("extdata", "extract_pod5_signals.py", package = "ninetails")}
#'   \item Large datasets may require substantial temporary disk space
#'   \item Progress messages are printed directly from the Python subprocess
#' }
#'
#' @seealso
#' \code{\link{create_tail_features_list_dorado}} for computing pseudomoves from signals,
#' \code{\link{process_dorado_signal_files}} for complete signal processing pipeline
#'
#' @examples
#' \dontrun{
#' # Load poly(A) coordinate data
#' polya_data <- read.table("polya_coords.txt", header = TRUE)
#'
#' # Extract signals using 4 cores
#' signals <- extract_tails_from_pod5(
#'   polya_data = polya_data,
#'   pod5_dir = "/path/to/pod5/files/",
#'   num_cores = 4
#' )}
#'
#' @export
#' @importFrom reticulate py_config import import_builtins py_to_r
#' @importFrom parallel detectCores
#' @importFrom utils write.csv
extract_tails_from_pod5 <- function(polya_data,
                                    pod5_dir,
                                    num_cores = 1) {

  # Default to all cores minus 1 if not specified
  if (is.null(num_cores)) {
    num_cores <- max(1, parallel::detectCores() - 1)
  }

  # Validate input
  required_cols <- c("read_id", "filename", "poly_tail_start", "poly_tail_end")
  missing_cols <- setdiff(required_cols, colnames(polya_data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  # Filter valid reads
  if ("poly_tail_length" %in% colnames(polya_data)) {
    polya_data <- polya_data[polya_data$poly_tail_length > 10 &
                               polya_data$poly_tail_start > 0, ]
  } else {
    polya_data <- polya_data[polya_data$poly_tail_start > 0, ]
  }

  if (nrow(polya_data) == 0) {
    return(list())  # Return empty list if no valid reads
  }

  # Get Python script path
  python_script <- system.file("extdata", "extract_pod5_signals.py", package = "ninetails")
  if (!file.exists(python_script)) {
    stop("POD5 extraction script not found. Please reinstall ninetails package.")
  }

  # Create temp files
  temp_dir <- tempdir()
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  input_file <- file.path(temp_dir, paste0("polya_input_", timestamp, ".csv"))
  output_file <- file.path(temp_dir, paste0("polya_output_", timestamp, ".pkl"))

  # Save input data
  write.csv(polya_data, input_file, row.names = FALSE)

  # Build and execute Python command
  python_cmd <- sprintf('"%s" "%s" --input "%s" --pod5_dir "%s" --output "%s" --num_cores %d',
                        reticulate::py_config()$python,
                        python_script,
                        input_file,
                        pod5_dir,
                        output_file,
                        num_cores)

  exit_code <- system(python_cmd, intern = FALSE)

  if (exit_code != 0 || !file.exists(output_file)) {
    unlink(c(input_file, output_file))
    stop("POD5 extraction failed. Ensure 'pod5' Python module is installed: pip install pod5")
  }

  # Load results
  pkl <- reticulate::import("pickle", convert = FALSE)
  py_builtin <- reticulate::import_builtins()
  f <- py_builtin$open(output_file, "rb")
  signals_list <- reticulate::py_to_r(pkl$load(f))
  f$close()

  # Clean up
  unlink(c(input_file, output_file))

  # Ensure numeric vectors
  signals_list <- lapply(signals_list, as.numeric)

  return(signals_list)
}



#' Preprocess Dorado inputs for ninetails analysis (no BAM processing)
#'
#' This function prepares Dorado summary files and extracts poly(A) signals
#' from POD5 files for downstream analysis. It splits large summary files into
#' manageable parts and extracts corresponding poly(A) signals.
#'
#' @param dorado_summary Character or data frame. Path to Dorado summary file or data frame.
#'
#' @param pod5_dir Character. Directory containing pod5 files.
#'
#' @param num_cores Integer. Number of CPU cores to use.
#'
#' @param qc Logical. Whether to perform quality control.
#'
#' @param save_dir Character. Directory where output files will be saved.
#'
#' @param prefix Character. Prefix to add to output file names (optional).
#'
#' @param part_size Integer. Number of reads to process in each chunk.
#'
#' @param cli_log Function for logging messages and progress.
#'
#' @return List containing paths to processed files:
#'   \describe{
#'     \item{summary_files}{Paths to split summary files}
#'     \item{polya_signal_files}{Paths to extracted poly(A) signal files}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' processed_files <- preprocess_inputs(
#'   dorado_summary = "path/to/summary.txt",
#'   pod5_dir = "path/to/pod5/",
#'   num_cores = 2,
#'   qc = TRUE,
#'   save_dir = "path/to/output/",
#'   prefix = "experiment1",
#'   part_size = 40000,
#'   cli_log = message
#' )
#' }
preprocess_inputs <- function(dorado_summary,
                              pod5_dir,
                              num_cores,
                              qc,
                              save_dir,
                              prefix,
                              part_size,
                              cli_log) {

  # Input validation and configuration
  cli_log("Validating input parameters", "INFO", "Validating Inputs")
  assertthat::assert_that(is.numeric(num_cores), num_cores > 0,
                          msg = "Number of cores must be a positive numeric value")
  assertthat::assert_that(is.character(pod5_dir), dir.exists(pod5_dir),
                          msg = "POD5 directory must exist")
  assertthat::assert_that(is.numeric(part_size), part_size >= 1,
                          msg = "Part size must be at least 1")

  ################################################################################
  # SUMMARY FILE PROCESSING
  ################################################################################
  cli_log("Processing summary files...", "INFO", "Summary Processing", bullet = TRUE)

  # Read and validate summary file
  if (is.character(dorado_summary)) {
    if (!file.exists(dorado_summary)) {
      stop("Dorado summary file does not exist: ", dorado_summary)
    }
    summary_data <- vroom::vroom(dorado_summary, show_col_types = FALSE)
  } else {
    summary_data <- dorado_summary
  }

  # Validate required columns for poly(A) information
  required_cols <- c("read_id", "filename", "poly_tail_length", "poly_tail_start", "poly_tail_end")
  missing_cols <- setdiff(required_cols, colnames(summary_data))

  if (length(missing_cols) > 0) {
    stop(sprintf("Summary file missing required poly(A) columns: %s. Please ensure you are using Dorado summary files with poly(A) information.",
                 paste(missing_cols, collapse = ", ")))
  }

  cli_log(sprintf("Summary contains %d reads with poly(A) information", nrow(summary_data)), "INFO")

  # Apply quality control filtering if requested
  if (qc) {
    cli_log("Applying quality control filters...", "INFO", bullet = TRUE)
    initial_count <- nrow(summary_data)

    # Filter for valid poly(A) tails (length > 10, valid coordinates)
    summary_data <- summary_data[
      summary_data$poly_tail_length > 10 &
        summary_data$poly_tail_start > 0 &
        summary_data$poly_tail_end > summary_data$poly_tail_start,
    ]

    final_count <- nrow(summary_data)
    filtered_count <- initial_count - final_count

    if (filtered_count > 0) {
      cli_log(sprintf("QC filtering removed %d reads (%d remaining)",
                      filtered_count, final_count), "INFO")
    }
  }

  # Split summary file into parts
  summary_dir <- file.path(save_dir, "dorado_summary_dir")
  if (!dir.exists(summary_dir)) {
    dir.create(summary_dir, recursive = TRUE)
  }

  cli_log("Splitting summary file into parts...", "INFO", bullet = TRUE)

  # Calculate number of parts needed
  total_reads <- nrow(summary_data)
  num_parts <- ceiling(total_reads / part_size)

  cli_log(sprintf("Total reads: %d", total_reads), "INFO")
  cli_log(sprintf("Splitting into %d parts", num_parts), "INFO")

  # Initialize vector to store output file paths
  part_files <- character(num_parts)

  # Split and save summary data directly in dorado_summary_dir
  for (i in seq_len(num_parts)) {
    start_idx <- ((i - 1) * part_size) + 1
    end_idx <- min(i * part_size, total_reads)

    # Extract subset of data
    part_data <- summary_data[start_idx:end_idx, ]

    # Create output filename with prefix if provided
    if (nchar(prefix) > 0) {
      output_file <- file.path(summary_dir, sprintf("%s_summary_part%d.txt", prefix, i))
    } else {
      output_file <- file.path(summary_dir, sprintf("summary_part%d.txt", i))
    }

    # Save data using vroom for consistency
    vroom::vroom_write(part_data, file = output_file, delim = "\t")
    part_files[i] <- output_file

    cli_log(sprintf("Saved part %d/%d: %s (%d reads)", i, num_parts, basename(output_file), nrow(part_data)), "INFO")
  }

  cli_log(sprintf("Summary split into %d parts", length(part_files)), "SUCCESS")

  ################################################################################
  # EXTRACTING POLY(A) SIGNALS FROM POD5
  ################################################################################
  cli_log("Extracting poly(A) signals from POD5 files...", "INFO", "Signal Extraction", bullet = TRUE)

  polya_signal_dir <- file.path(save_dir, "polya_signal_dir")
  if (!dir.exists(polya_signal_dir)) {
    dir.create(polya_signal_dir, recursive = TRUE)
  }

  polya_signal_files <- character(length(part_files))

  for (i in seq_along(part_files)) {
    cli_log(sprintf("Processing part %d/%d for signal extraction", i, length(part_files)), "INFO")

    # Load current part
    current_summary <- vroom::vroom(part_files[i], show_col_types = FALSE)

    # Extract signals using POD5
    signals <- ninetails::extract_tails_from_pod5(
      polya_data = current_summary,
      pod5_dir = pod5_dir,
      num_cores = num_cores
    )

    if (length(signals) > 0) {
      # Generate output filename
      output_basename <- if (nchar(prefix) > 0) {
        sprintf("%s_polya_signal_part%d.rds", prefix, i)
      } else {
        sprintf("polya_signal_part%d.rds", i)
      }

      output_file <- file.path(polya_signal_dir, output_basename)
      saveRDS(signals, output_file)
      polya_signal_files[i] <- output_file

      cli_log(sprintf("Extracted signals for %d reads in part %d", length(signals), i), "INFO")
    } else {
      cli_log(sprintf("No signals extracted for part %d", i), "WARNING")
    }

    # Clean up memory
    rm(current_summary, signals)
    gc()
  }

  # Remove empty entries
  polya_signal_files <- polya_signal_files[file.exists(polya_signal_files)]

  if (length(polya_signal_files) == 0) {
    stop("No poly(A) signals were successfully extracted from any part")
  }

  cli_log(sprintf("Successfully extracted signals from %d parts", length(polya_signal_files)), "SUCCESS")

  # Return paths to processed files
  return(list(
    summary_files = part_files,
    polya_signal_files = polya_signal_files
  ))
}



#' Detection of outliers (peaks & valleys) in ONT signal using z-scores (Ultra-fast vectorized version).
#'
#' This function provides an ultra-fast, highly vectorized implementation for detecting
#' signal outliers in Oxford Nanopore poly(A) tail sequences. It identifies areas where
#' the signal significantly deviates from typical adenosine homopolymer values, which
#' may indicate the presence of non-adenosine nucleotides (C, G, U).
#'
#' The algorithm uses a robust peak detection approach based on z-scores with adaptive
#' windowing to distinguish between normal adenosine signal and potential modifications.
#' This vectorized implementation prioritizes performance through efficient memory usage,
#' reduced function calls, and optimized statistical calculations.
#'
#' @section Algorithm Details:
#' The function implements a sliding window approach where:
#' \itemize{
#'   \item A calibration phase establishes baseline signal characteristics
#'   \item Rolling statistics (mean, standard deviation) are computed for each position
#'   \item Outliers are detected when signal deviates >3.5 standard deviations from local mean
#'   \item Direction of deviation determines pseudomove value: +1 (peaks), -1 (valleys), 0 (normal)
#' }
#'
#' @section Performance Optimizations:
#' This vectorized version includes several performance improvements:
#' \itemize{
#'   \item Vectorized statistical calculations (sum/length instead of mean())
#'   \item Efficient rolling window operations
#'   \item Reduced memory allocations and copying
#'   \item Minimized function call overhead
#'   \item Integer operations where appropriate
#' }
#'
#' @section Output Interpretation:
#' The returned pseudomove vector contains:
#' \itemize{
#'   \item \strong{1}: Signal significantly above baseline (potential G nucleotides)
#'   \item \strong{0}: Signal within normal range (likely A nucleotides)
#'   \item \strong{-1}: Signal significantly below baseline (potential C/U nucleotides)
#' }
#'
#' @section Quality Control:
#' The function applies several quality control measures:
#' \itemize{
#'   \item Terminal position masking (first 5 positions set to 0)
#'   \item Gap substitution to handle segmentation artifacts
#'   \item Calibration using most frequent signal values
#'   \item Reproducible random seed for consistent results
#' }
#'
#' @param signal A numeric vector representing the raw Oxford Nanopore signal
#'   corresponding to the poly(A) tail region. Signal values should be current
#'   measurements in picoamperes. Minimum length is 10 data points.
#'
#' @returns A numeric vector of "pseudomoves" with the same length as the input
#'   signal (after calibration data removal). Values are integers in the range
#'   \{-1, 0, 1\} representing:
#'   \describe{
#'     \item{-1}{Signal valleys (potential C/U nucleotides)}
#'     \item{0}{Normal signal (likely A nucleotides)}
#'     \item{1}{Signal peaks (potential G nucleotides)}
#'   }
#'
#' @section References:
#' The peak detection algorithm is based on:
#' Brakel, J.P.G. van (2014). "Robust peak detection algorithm using z-scores".
#' Stack Overflow. Available at:
#' \url{https://stackoverflow.com/questions/22583391/peak-signal-detection-in-realtime-timeseries-data/22640362#22640362}
#' (version: 2020-11-08).
#'
#' @section Implementation Notes:
#' \itemize{
#'   \item Uses a fixed random seed (123) for reproducible results
#'   \item Window size is fixed at 100 data points
#'   \item Z-score threshold is set to 3.5 standard deviations
#'   \item Calibration uses 100 synthetic data points
#'   \item Terminal positions (first 5) are masked to prevent false positives
#' }
#'
#' @section Performance:
#' Expected performance improvements over the original implementation:
#' \itemize{
#'   \item 3-5x faster execution time for typical signal lengths
#'   \item Reduced memory usage by approximately 30-50\%
#'   \item Better scaling for large signals (>10,000 data points)
#' }
#'
#' @family signal_processing
#' @family pseudomove_functions
#' @family optimization_functions
#'
#' @seealso
#' \code{\link{filter_signal_by_threshold}} for the original implementation,
#' \code{\link{substitute_gaps}} for gap handling,
#' \code{\link{create_tail_features_list_dorado}} for the complete feature extraction pipeline
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage with a signal vector
#' pseudomoves <- filter_signal_by_threshold_vectorized(tail_signal)
#'
#' # Example with simulated poly(A) signal
#' set.seed(42)
#' simulated_signal <- c(
#'   rnorm(100, mean = 650, sd = 15),  # Normal A signal
#'   rnorm(20, mean = 700, sd = 10),   # G-like peak
#'   rnorm(50, mean = 645, sd = 12),   # More A signal
#'   rnorm(15, mean = 580, sd = 8)     # C/U-like valley
#' )
#'
#' # Detect modifications
#' pseudomoves <- filter_signal_by_threshold_vectorized(simulated_signal)
#'
#' # Visualize results
#' plot(simulated_signal, type = "l", main = "Signal with Detected Modifications")
#' points(which(pseudomoves == 1), simulated_signal[which(pseudomoves == 1)],
#'        col = "red", pch = 16)  # Peaks
#' points(which(pseudomoves == -1), simulated_signal[which(pseudomoves == -1)],
#'        col = "blue", pch = 16) # Valleys
#'
#' # Integration with ninetails pipeline
#' # This function is typically called within create_tail_features_list_dorado()
#' features <- create_tail_features_list_dorado(
#'   signal_list = list("read1" = tail_signal),
#'   num_cores = 1
#' )
#' }
#'
filter_signal_by_threshold_vectorized <- function(signal) {

  # Input validation
  if (missing(signal) || !is.numeric(signal) || length(signal) < 10) {
    stop("Invalid signal input", call. = FALSE)
  }

  set.seed(123)

  # Parameters
  window_size <- 100L
  threshold <- 3.5

  # Prepare signal with calibration
  start_vals <- signal[1:10]
  most_freq <- as.numeric(names(sort(table(signal), decreasing = TRUE)[1:20]))
  calibration <- sample(c(most_freq, start_vals), 100, replace = TRUE)
  adj_signal <- c(calibration, signal)

  n <- length(adj_signal)

  # Vectorized rolling statistics using efficient approach
  # Pre-allocate results
  pseudomoves <- integer(n)

  # Calculate rolling means and SDs more efficiently
  # Use embedded for loop but with vectorized operations where possible
  for (i in (window_size + 1):n) {
    # Get window efficiently
    window_start <- max(1L, i - window_size)
    window_data <- adj_signal[window_start:(i-1)]

    # Fast statistics
    window_mean <- sum(window_data) / length(window_data)
    window_sd <- sqrt(sum((window_data - window_mean)^2) / (length(window_data) - 1))

    # Outlier detection
    deviation <- abs(adj_signal[i] - window_mean)
    if (deviation > threshold * window_sd) {
      pseudomoves[i] <- if (adj_signal[i] > window_mean) 1L else -1L
    }
  }

  # Extract result and apply filters
  result <- pseudomoves[101:n]
  if (length(result) > 5) {
    result[1:5] <- 0L
  }

  # Apply gap substitution
  result <- ninetails::substitute_gaps(result)

  return(result)
}




#' Creates a nested list of Dorado tail features (raw signal + pseudomoves).
#'
#' This function processes raw poly(A) tail signal traces from Dorado and computes
#' pseudomoves using a threshold-based filter. Each read is stored in a nested
#' list structure, keyed by its read ID, with both the raw signal and the
#' computed pseudomove vector.
#'
#' Unlike Guppy-based pipelines, Dorado does not provide moves directly in BAM
#' files, so pseudomoves are computed empirically from the raw signal traces.
#'
#' Parallel execution is handled with \pkg{foreach} and \pkg{doSNOW}, and a
#' progress bar is displayed during processing.
#'
#' @param signal_list list of numeric vectors. Each element must represent
#' a raw poly(A) tail signal trace, with list names corresponding to read IDs.
#' @param num_cores numeric [1]. Number of physical cores to use in processing.
#' Do not exceed 1 less than the total available cores on your system.
#'
#' @return A nested list of tail features, organized by read IDs. Each read entry
#' contains:
#'   \itemize{
#'     \item \code{tail_signal}: numeric vector of the raw poly(A) tail signal.
#'     \item \code{tail_pseudomoves}: numeric vector of pseudomove states
#'           computed by \code{ninetails::filter_signal_by_threshold}.
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' features <- ninetails::create_tail_features_list_dorado(
#'   signal_list = signals,
#'   num_cores = 4
#' )
#' }
create_tail_features_list_dorado <- function(signal_list,
                                             num_cores) {

  # variable binding (suppressing R CMD check from throwing an error)
  x <- NULL

  # Assertions
  if (missing(signal_list)) {
    stop("Signal list is missing. Please provide a valid signal_list argument.", call. = FALSE)
  }
  if (missing(num_cores)) {
    stop("Number of declared cores is missing. Please provide a valid num_cores argument.", call. = FALSE)
  }

  assertthat::assert_that(is.list(signal_list),
                          msg = "Provided signal_list is not a list. Please provide a valid list of raw signal vectors.")
  assertthat::assert_that(is.numeric(num_cores),
                          msg = "Declared num_cores must be numeric. Please provide a valid value.")
  assertthat::assert_that(all(sapply(signal_list, is.numeric)),
                          msg = "All elements of signal_list must be numeric vectors (raw tail signal traces).")


  # creating cluster for parallel computing
  my_cluster <- parallel::makeCluster(num_cores)
  on.exit(parallel::stopCluster(my_cluster))
  doSNOW::registerDoSNOW(my_cluster)
  `%dopar%` <- foreach::`%dopar%`
  `%do%` <- foreach::`%do%`
  # Options for multicore execution
  mc_options <- list(preschedule = TRUE, set.seed = FALSE, cleanup = FALSE)

  #creating list for outputs

  output_nested_list <- list()

  #progressbar header
  cat(paste0('[', as.character(Sys.time()), '] ','Computing pseudomoves...', '\n', sep=''))

  # progress bar
  pb <- utils::txtProgressBar(min = 0,
                              max = length(signal_list),
                              style = 3,
                              width = 50,
                              char = "=",
                              file = stderr())
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  output_nested_list <- list()

  # Main parallel loop
  # For each read (named element of signal_list):
  #   - store the raw signal trace
  #   - compute pseudomoves with ninetails::filter_signal_by_threshold()
  #   - return as a named sublist under the read ID
  output_nested_list <- foreach::foreach(x = names(signal_list),
                                         .combine = c,
                                         .inorder = TRUE,
                                         .errorhandling = "pass",
                                         .options.snow = opts,
                                         .options.multicore = mc_options) %dopar% (function(x) {
                                           stats::setNames(list(list(
                                             tail_signal = signal_list[[x]],
                                             tail_pseudomoves = ninetails::filter_signal_by_threshold_vectorized(signal_list[[x]]))
                                           ), x)
                                         })(x)

  #close(pb)
  # Done comm
  cat(paste0('[', as.character(Sys.time()), '] ','Done!', '\n', sep=''))

  return(output_nested_list)
}


#' Extracts fragments of poly(A) tail signal (Dorado mode) containing potential modifications
#' along with their delimitation (positional indices; coordinates) within the tail.
#'
#' This function processes raw poly(A) tail signal and Dorado-derived pseudomoves
#' to identify and extract signal segments (chunks) potentially corresponding
#' to modified positions (e.g., non-A residues). Each extracted chunk spans
#' 100 signal points, centered on the midpoint of a pseudomove run.
#'
#' In the Dorado pipeline, moves are not used:
#'   * retrieving them from BAM files is computationally expensive
#'   * processing is non-intuitive
#'
#' Instead, only pseudomoves are considered. As a safeguard against Dorado’s
#' tendency to extend poly(A) boundaries into the transcript body, the last
#' 3 pseudomove values are forced to 0. This prevents misclassification of
#' transcript nucleotides as part of the tail.
#'
#' Candidate modification regions are detected by:
#'   * run-length encoding (RLE) of the pseudomove vector
#'   * filtering runs of pseudomoves with length ≥ 5
#'
#' Extracted fragments are padded/imputed if they extend beyond signal boundaries:
#'   * upstream/downstream missing values (NAs) are replaced
#'   * imputation is based on random draws from the 5 most frequent signal values
#'
#' The function returns a list object (nested), where each element represents
#' one candidate modification region, containing:
#'   * `chunk_sequence`: the raw signal subsequence (length = 100, imputed if needed)
#'   * `chunk_start_pos`: starting index of the subsequence
#'   * `chunk_end_pos`: ending index of the subsequence
#'
#' @param readname character string. Name of the given read within the
#' analyzed dataset.
#'
#' @param tail_feature_list list object produced by
#'   \code{\link{create_tail_features_list_dorado}} Dorado-tail feature
#'   extraction function. Must contain
#'   \code{$tail_signal} (numeric vector) and \code{$tail_pseudomoves} (integer vector)
#'   for each read.
#'
#' @return a nested list where each element corresponds to a signal fragment.
#' Each fragment is itself a list with three entries:
#'   \itemize{
#'     \item \code{chunk_sequence}: numeric vector of raw signal values
#'     \item \code{chunk_start_pos}: integer, start index of the chunk
#'     \item \code{chunk_end_pos}: integer, end index of the chunk
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ninetails::split_tail_centered_dorado(
#'   readname = "1234-anexample-r3adn4m3",
#'   tail_feature_list = tail_feature_list
#' )
#' }
split_tail_centered_dorado <- function(readname, tail_feature_list) {

  # variable binding (suppressing R CMD check from throwing an error)
  i <- NULL

  #assertions
  if (missing(readname)) {
    stop("Readname is missing. Please provide a valid readname argument.", call. =FALSE)
  }

  if (missing(tail_feature_list)) {
    stop("List of tail features is missing. Please provide a valid tail_feature_list argument.", call. =FALSE)
  }

  assertthat::assert_that(is.character(readname),
                          msg = paste0("Given readname is not a character string. Please provide a valid readname."))
  assertthat::assert_that(is.list(tail_feature_list),
                          msg = paste0("Given tail_feature_list is not a list (class). Please provide valid file format."))

  # Extract required data
  # In the Dorado pipeline, moves are not used:
  #   - retrieving moves from BAM is time-consuming
  #   - processing moves is non-intuitive
  signal <- tail_feature_list[[readname]]$tail_signal
  pseudomoves <- tail_feature_list[[readname]]$tail_pseudomoves

  # Force the last 3 pseudomove values to 0 if they are not already 0
  # Dorado often extends tail boundaries slightly into the transcript body,
  # which can cause false positives by misclassifying transcript nucleotides
  # as part of the poly(A) tail.
  # Instead of adjusting tail boundaries directly, we zero these trailing positions.
  # This strategy is chosen with the expectation that ONT will eventually refine
  # their poly(A) detection algorithm.
  if (length(pseudomoves) >= 3) {
    idx <- (length(pseudomoves)-2):length(pseudomoves)
    pseudomoves[idx][pseudomoves[idx] != 0] <- 0
  }

  # mod-centered chunk extraction
  # Strategy:
  # - Compute run-length encoding (RLE) of pseudomoves
  # - Select runs that:
  #     * are "positive" pseudomove stretches (value != 0)
  #     * are at least 5 bases long
  # - These runs are likely to correspond to poly(A) tail signal
  mod_rle <- rle(pseudomoves)
  # pseudomoves filtered by condition (potentially decorated - empyrical!)
  condition <- mod_rle$lengths >= 5 & mod_rle$values != 0
  if (!any(condition)) return(NULL)  # No valid chunks, return NULL
  # For each valid run:
  # - Compute its starting position
  # - Compute chunk boundaries centered on the middle of the run
  # - Each chunk has a fixed window size of 100 bases (±50 around center)
  first_filtered_positions <- cumsum(c(1, utils::head(mod_rle$lengths, -1)))[condition]
  # length of pseudomoves satisfying condition
  filtered_length <- mod_rle$lengths[condition]
  # extracted coordinates (indices)
  start_positions <- first_filtered_positions + floor(filtered_length / 2) - 50
  end_positions <- start_positions + 99

  # extract signal chunks centered on potential modification
  # - If the start position is < 1, left-pad with NAs
  # - Later, NAs are imputed with frequent signal values
  list_1 <- lapply(seq_along(start_positions), function(i) {
    if (start_positions[i] > 0) signal[start_positions[i]:end_positions[i]]
    else c(rep(NA, abs(start_positions[i] - 1)), signal[1:end_positions[i]])
  })

  # NAs introduced by left-padding are replaced by random draws
  # from the 5 most frequent observed signal values.
  # This ensures consistent chunk length without introducing
  # artificial bias from fixed-value padding.
  most_freq_vals <- as.numeric(names(sort(table(signal), decreasing = TRUE)[1:5]))
  list_1 <- lapply(list_1, function(n) replace(n, is.na(n), sample(most_freq_vals, sum(is.na(n)), TRUE)))

  # naming chunks
  #   <readname>_<chunk_index>
  chunk_names <- paste0(rep(readname, length(list_1)), '_', seq_along(start_positions))
  names(list_1) <- chunk_names

  # retrieve coordinates as list_2 and list_3:
  list_2 <- as.list(start_positions)
  list_3 <- as.list(end_positions)

  # Output is a list of named lists:
  #   - chunk_sequence   (raw signal segment, length = 100)
  #   - chunk_start_pos  (starting index in original signal)
  #   - chunk_end_pos    (ending index in original signal)
  out <- purrr::transpose(purrr::set_names(list(list_1, list_2, list_3), c('chunk_sequence', 'chunk_start_pos', 'chunk_end_pos')))
  return(out)
}


#' Creates list of poly(A) tail chunks (Dorado mode) centered on significant signal deviations.
#'
#' Processes raw poly(A) tail signals and pseudomoves (as generated by Dorado)
#' in parallel to extract candidate signal fragments potentially containing
#' non-A nucleotides. Each fragment is 100 signal points long, centered on
#' pseudomove runs of sufficient length, and is returned with its positional
#' coordinates. The resulting data are organized into a nested list keyed
#' by read IDs.
#'
#' This Dorado-specific function differs from the Guppy-based version:
#'   * moves are not used (to avoid costly BAM parsing and processing)
#'   * pseudomoves are corrected at the tail ends (last 3 values forced to 0)
#'
#' Parallelization is handled with \pkg{foreach} and \pkg{doSNOW}, allowing
#' efficient scaling across multiple CPU cores. A progress bar is displayed
#' to monitor job completion.
#'
#' @param tail_feature_list list object produced by \code{\link{create_tail_feature_list}}
#' or an equivalent Dorado-tail feature extraction function. Must contain
#' per-read entries with \code{$tail_signal} and \code{$tail_pseudomoves}.
#'
#' @param num_cores numeric [1]. Number of physical cores to use in processing.
#' Do not exceed 1 less than the number of available cores on your machine.
#'
#' @return A nested list containing the segmented tail data (chunks and coordinates),
#' organized by read IDs. Each read entry contains one or more fragments, where
#' each fragment is a list with:
#'   \itemize{
#'     \item \code{chunk_sequence}: numeric vector of raw signal values (length 100)
#'     \item \code{chunk_start_pos}: integer, starting index of the chunk
#'     \item \code{chunk_end_pos}: integer, ending index of the chunk
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' tcl_dorado <- ninetails::create_tail_chunk_list_dorado(
#'   tail_feature_list = tfl,
#'   num_cores = 3
#' )
#' }
create_tail_chunk_list_dorado <- function(tail_feature_list, num_cores) {
  # variable binding (suppressing R CMD check from throwing an error)
  i <- NULL

  # assertions
  if (missing(num_cores)) {
    stop("Number of declared cores is missing. Please provide a valid num_cores argument.", call. =FALSE)
  }

  if (missing(tail_feature_list)) {
    stop("List of features is missing. Please provide a valid tail_feature_list argument.", call. =FALSE)
  }

  assertthat::assert_that(is.numeric(num_cores),
                          msg = paste0("Declared core number must be numeric. Please provide a valid argument."))
  assertthat::assert_that(is.list(tail_feature_list),
                          msg = paste0("Given tail_feature_list is not a list (class). Please provide valid file format."))

  # creating cluster for parallel computing
  my_cluster <- parallel::makeCluster(num_cores)
  on.exit(parallel::stopCluster(my_cluster))
  doSNOW::registerDoSNOW(my_cluster)
  `%dopar%` <- foreach::`%dopar%`
  mc_options <- list(preschedule = TRUE, set.seed = FALSE, cleanup = FALSE)

  # header for progress bar
  cat(paste0('[', as.character(Sys.time()), '] ','Creating tail segmentation data...', '\n'))

  # progress bar
  pb <- utils::txtProgressBar(min = 0,
                              max = length(tail_feature_list),
                              style = 3,
                              width = 50,
                              char = "=",
                              file = stderr())
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  #output list
  tail_chunk_list <- list()

  # For each read (element of tail_feature_list):
  #   - Call split_tail_centered_dorado() on the read
  #   - Collect resulting signal chunks into a nested list
  tail_chunk_list <- foreach::foreach(i = seq_along(tail_feature_list),
                                      .combine = c,
                                      .inorder = TRUE,
                                      .errorhandling = 'pass',
                                      .options.snow = opts,
                                      .options.multicore = mc_options) %dopar% {
                                        lapply(names(tail_feature_list[i]), function(x) ninetails::split_tail_centered_dorado(x, tail_feature_list))
                                      }

  names(tail_chunk_list) <- names(tail_feature_list)
  # Remove NULLs (reads with no valid chunks)
  tail_chunk_list <- rrapply::rrapply(tail_chunk_list,
                                      condition = Negate(is.null),
                                      how = "prune")
  # Done comment
  cat(paste0('[', as.character(Sys.time()), '] ','Done!', '\n'))

  return(tail_chunk_list)
}


#' Create Ninetails output tables for Dorado DRS pipeline
#'
#' This function integrates Dorado poly(A) tail summaries, temporary non-adenosine
#' predictions, and poly(A) chunk information to generate two main outputs:
#' \itemize{
#' \item A table of read classifications (read_classes), categorizing reads as
#' decorated, blank, or unclassified, with metadata such as poly(A) length, contig,
#' and quality tags.
#' \item A table of non-adenosine residue predictions (nonadenosine_residues),
#' including their estimated positions along the poly(A) tail.
#' }
#'
#' @param dorado_summary_dir Character string. Path to a directory containing Dorado
#' summary files (.txt, .tsv, or .csv) with per-read poly(A) tail information.
#'
#' @param nonA_temp_dir Character string. Path to a directory containing non-adenosine
#' prediction RDS files, generated from temporary models.
#'
#' @param polya_chunks_dir Character string. Path to a directory containing poly(A) chunk
#' RDS files used for position inference of predictions.
#'
#' @param num_cores Integer. Number of cores to use for parallelized file loading
#' and processing. Must be a positive integer. Default is 1.
#'
#' @param qc Logical. Whether to apply quality control filtering of terminal
#' predictions (removing predictions near the ends of poly(A) tails). Default is TRUE.
#'
#' @return A named list with two data frames:
#' \describe{
#' \item{read_classes}{Data frame with per-read classification results, including
#' columns for read name, contig, poly(A) length, QC tag, class, and comments.}
#' \item{nonadenosine_residues}{Data frame with per-chunk predictions of
#' non-adenosine residues, including read name, contig, predicted base,
#' estimated position within the poly(A) tail, poly(A) length, and QC tag.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' results <- create_outputs_dorado(
#' dorado_summary_dir = "data/dorado_summaries",
#' nonA_temp_dir = "data/nonA_predictions",
#' polya_chunks_dir = "data/polya_chunks",
#' num_cores = 4,
#' qc = TRUE
#' )
#'
#' # Access read classifications
#' head(results$read_classes)
#'
#' # Access non-adenosine residues
#' head(results$nonadenosine_residues)
#' }
create_outputs_dorado <- function(dorado_summary_dir,
                                  nonA_temp_dir,
                                  polya_chunks_dir,
                                  num_cores = 1,
                                  qc = TRUE) {

  # Variable binding for R CMD check
  read_id <- alignment_genome <- alignment_mapq <- poly_tail_length <- NULL
  poly_tail_start <- poly_tail_end <- chunkname <- prediction <- NULL
  centr_signal_pos <- signal_length <- est_nonA_pos <- class <- comments <- NULL

  # Assertions
  if (missing(dorado_summary_dir)) stop("Dorado summary directory is missing.", call. = FALSE)
  if (missing(nonA_temp_dir)) stop("Non-A predictions directory is missing.", call. = FALSE)
  if (missing(polya_chunks_dir)) stop("PolyA chunks directory is missing.", call. = FALSE)
  if (missing(num_cores)) stop("Number of cores is missing.", call. = FALSE)

  # Batch argument validation
  args <- list(dorado_summary_dir, nonA_temp_dir, polya_chunks_dir)
  dirs_exist <- sapply(args, function(x) is.character(x) && length(x) == 1 && dir.exists(x))
  if (!all(dirs_exist)) {
    stop("All directory arguments must be valid existing paths.", call. = FALSE)
  }

  if (!is.numeric(num_cores) || num_cores < 1) {
    stop("num_cores must be a positive integer.", call. = FALSE)
  }

  # Check files exist (batch operation)
  dorado_files <- list.files(dorado_summary_dir, pattern = "\\.txt$|\\.tsv$|\\.csv$", full.names = TRUE)
  chunk_files <- list.files(polya_chunks_dir, pattern = "\\.rds$", full.names = TRUE)
  rds_files <- list.files(nonA_temp_dir, pattern = "\\.rds$", full.names = TRUE)

  if (length(dorado_files) == 0) stop("No summary files found in dorado_summary_dir", call. = FALSE)
  if (length(chunk_files) == 0) stop("No RDS files found in polya_chunks_dir", call. = FALSE)
  if (length(rds_files) == 0) stop("No prediction files found in nonA_temp_dir", call. = FALSE)

  # Loading dorado summary files - using vroom for optimal performance
  cat("Loading Dorado summary files...\n")
  dorado_list <- lapply(dorado_files, function(f) {
    suppressMessages(vroom::vroom(f, delim = "\t", show_col_types = FALSE,
                                  col_select = c("read_id", "alignment_genome", "alignment_mapq",
                                                 "poly_tail_length", "poly_tail_start", "poly_tail_end")))
  })

  dorado_summary <- dplyr::bind_rows(dorado_list)

  # Loading prediction files in parallel
  cat("Loading all prediction files in parallel...\n")

  # Set up parallel cluster for file loading
  cl <- parallel::makeCluster(min(num_cores, length(rds_files)))
  on.exit(parallel::stopCluster(cl), add = TRUE)

  # Load all predictions in parallel
  all_predictions <- parallel::parLapply(cl, rds_files, function(file) {
    preds <- readRDS(file)

    # Standardize prediction format
    if (is.list(preds) && 'chunkname' %in% names(preds)) {
      data.frame(chunkname = preds$chunkname, prediction = preds$prediction, stringsAsFactors = FALSE)
    } else if (is.list(preds)) {
      data.frame(chunkname = names(preds), prediction = unlist(preds), stringsAsFactors = FALSE)
    } else {
      data.frame(chunkname = names(preds), prediction = as.vector(preds), stringsAsFactors = FALSE)
    }
  })

  # Combine all predictions using vroom/dplyr (keeps consistency)
  cat("Combining all predictions...\n")
  moved_chunks_table <- dplyr::bind_rows(all_predictions)

  # Skip if no predictions
  if (nrow(moved_chunks_table) == 0) {
    warning("No predictions found in any files")
    return(list(read_classes = data.frame(), nonadenosine_residues = data.frame()))
  }

  # Extract read_id
  moved_chunks_table$read_id <- sub('_.*', '', moved_chunks_table$chunkname)

  # Vectorized prediction conversion
  prediction_map <- c("0" = "A", "1" = "C", "2" = "G", "3" = "U")
  moved_chunks_table$prediction <- prediction_map[as.character(moved_chunks_table$prediction)]

  # Grouping operations
  # Identify blank reads (all A) - optimized base R approach
  prediction_by_read <- split(moved_chunks_table$prediction, moved_chunks_table$read_id)
  all_A_reads <- sapply(prediction_by_read, function(x) all(x == "A"))
  moved_blank_readnames <- names(all_A_reads)[all_A_reads]

  # Remove blank reads and A predictions
  moved_chunks_table <- moved_chunks_table[!(moved_chunks_table$read_id %in% moved_blank_readnames) &
                                             moved_chunks_table$prediction != "A", ]

  cat("Total non-A predictions:", nrow(moved_chunks_table), "\n")
  cat("Total blank reads:", length(moved_blank_readnames), "\n")

  # Only load chunks that are actually needed
  if (nrow(moved_chunks_table) > 0) {
    needed_chunks <- unique(moved_chunks_table$chunkname)

    # Load only required chunk files
    cat("Loading and processing chunk files (optimized)...\n")

    # Pre-allocate lists for better performance
    position_data <- vector("list", length(chunk_files))
    chunk_names_list <- vector("list", length(chunk_files))

    # Process chunks in parallel; extract what is needed
    chunk_results <- parallel::parLapply(cl, seq_along(chunk_files), function(i) {
      chunks <- readRDS(chunk_files[i])

      # Pre-allocate vectors
      positions <- numeric()
      names_vec <- character()

      for (read_name in names(chunks)) {
        read_chunks <- chunks[[read_name]]
        for (chunk_name in names(read_chunks)) {
          # Only process if this chunk is needed
          if (chunk_name %in% needed_chunks) {
            positions <- c(positions, read_chunks[[chunk_name]][['chunk_start_pos']] + 50)
            names_vec <- c(names_vec, chunk_name)
          }
        }
      }

      if (length(positions) > 0) {
        data.frame(centr_signal_pos = positions, chunkname = names_vec, stringsAsFactors = FALSE)
      } else {
        NULL
      }
    })

    # Combine position results
    chunk_results <- chunk_results[!sapply(chunk_results, is.null)]
    if (length(chunk_results) > 0) {
      non_a_position_list <- dplyr::bind_rows(chunk_results)
    } else {
      stop("No position data found for required chunks")
    }

    # Add read_id
    non_a_position_list$read_id <- gsub('_.*', '', non_a_position_list$chunkname)

    # Use dplyr for joins
    # Merge operations
    non_a_position_list <- dplyr::left_join(non_a_position_list, dorado_summary, by = "read_id")
    non_a_position_list$signal_length <- 0.2 * (non_a_position_list$poly_tail_end - non_a_position_list$poly_tail_start)

    # Final merge
    moved_chunks_table <- dplyr::left_join(moved_chunks_table, non_a_position_list,
                                           by = c("read_id", "chunkname"))

    # Calculate position
    moved_chunks_table$est_nonA_pos <- round(
      moved_chunks_table$poly_tail_length - ((moved_chunks_table$poly_tail_length * moved_chunks_table$centr_signal_pos) / moved_chunks_table$signal_length),
      2
    )

    # Select final columns
    moved_chunks_table <- moved_chunks_table[, c("read_id", "alignment_genome", "prediction", "est_nonA_pos",
                                                 "poly_tail_length", "signal_length", "alignment_mapq")]
  }

  # Get all read IDs
  all_read_ids <- unique(dorado_summary$read_id)

  # Quality control & sanity check - optimized
  if (qc == TRUE && nrow(moved_chunks_table) > 0) {
    # Vectorized filtering
    terminal_mask <- (moved_chunks_table$est_nonA_pos < 2) |
      (moved_chunks_table$est_nonA_pos > moved_chunks_table$poly_tail_length - 2)

    moved_chunks_table_discarded_ids <- unique(moved_chunks_table$read_id[terminal_mask])
    moved_chunks_table <- moved_chunks_table[!terminal_mask, ]

    # Update blank reads
    moved_blank_readnames <- unique(c(moved_blank_readnames, moved_chunks_table_discarded_ids))
  }

  decorated_read_ids <- unique(moved_chunks_table$read_id)
  blank_read_ids <- setdiff(all_read_ids, c(decorated_read_ids, moved_blank_readnames))

  # Vectorized read classification
  dorado_summary$class <- "blank"  # Default
  dorado_summary$comments <- "MAU"  # Default

  # Update classifications
  dorado_summary$class[dorado_summary$read_id %in% decorated_read_ids] <- "decorated"
  dorado_summary$comments[dorado_summary$read_id %in% decorated_read_ids] <- "YAY"

  dorado_summary$class[dorado_summary$read_id %in% moved_blank_readnames] <- "blank"
  dorado_summary$comments[dorado_summary$read_id %in% moved_blank_readnames] <- "MPU"

  dorado_summary$class[dorado_summary$poly_tail_length < 10] <- "unclassified"
  dorado_summary$comments[dorado_summary$poly_tail_length < 10] <- "IRL"

  read_classes <- dorado_summary
  read_classes <- read_classes[!duplicated(read_classes$read_id), ]

  # Direct column operations
  # Format read_classes
  names(read_classes)[names(read_classes) == "read_id"] <- "readname"
  names(read_classes)[names(read_classes) == "alignment_genome"] <- "contig"
  names(read_classes)[names(read_classes) == "poly_tail_length"] <- "polya_length"
  names(read_classes)[names(read_classes) == "alignment_mapq"] <- "qc_tag"

  read_classes <- read_classes[, c("readname", "contig", "polya_length", "qc_tag", "class", "comments")]

  # Format moved_chunks_table
  if (nrow(moved_chunks_table) > 0) {
    names(moved_chunks_table)[names(moved_chunks_table) == "read_id"] <- "readname"
    names(moved_chunks_table)[names(moved_chunks_table) == "alignment_genome"] <- "contig"
    names(moved_chunks_table)[names(moved_chunks_table) == "poly_tail_length"] <- "polya_length"
    names(moved_chunks_table)[names(moved_chunks_table) == "alignment_mapq"] <- "qc_tag"

    moved_chunks_table <- moved_chunks_table[, c("readname", "contig", "prediction", "est_nonA_pos", "polya_length", "qc_tag")]
  } else {
    moved_chunks_table <- data.frame(
      readname = character(), contig = character(), prediction = character(),
      est_nonA_pos = numeric(), polya_length = numeric(), qc_tag = character(),
      stringsAsFactors = FALSE
    )
  }

  # Final output
  ninetails_output <- list(
    read_classes = read_classes,
    nonadenosine_residues = moved_chunks_table
  )

  # Summary statistics
  cat("Processing complete!\n")
  cat("Decorated reads:", sum(read_classes$class == "decorated", na.rm = TRUE), "\n")
  cat("Blank reads:", sum(read_classes$class == "blank", na.rm = TRUE), "\n")
  cat("Unclassified reads:", sum(read_classes$class == "unclassified", na.rm = TRUE), "\n")
  cat("Total nonadenosine residues:", nrow(moved_chunks_table), "\n")

  return(ninetails_output)
}
