#' Extract poly(A) tail signals from POD5 files using parallel processing
#'
#' @param polya_data Data frame with columns: read_id, filename, poly_tail_start, poly_tail_end
#' @param pod5_dir Character. Path to the directory containing POD5 files
#' @param num_cores Integer. Number of cores to use (default: 1)
#'
#' @return List of numeric vectors containing extracted signals
#' @export
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





