#' Process and split Dorado summary file into smaller parts
#'
#' @param dorado_summary Path to Dorado summary file or a data frame containing summary information
#' @param save_dir Directory where split summary files will be saved
#' @param part_size Number of reads per file part when splitting the summary file
#' @param cli_log Function for logging messages and progress
#'
#' @returns Character vector containing paths to the split summary files
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
    stop("Number of reads per file chunk (part_size) is missing. Please provide a valid part_size argument.", call. = FALSE)
  }

  assertthat::assert_that(is.numeric(part_size), part_size > 0,
                          msg = paste0("Reads per chunk must be numeric and positive. Please provide a valid argument."))

  # Handle input data and get prefix
  if (is.character(dorado_summary)) {
    cli_log(sprintf("Reading summary file: %s", basename(dorado_summary)), "INFO")

    # First, let's peek at the file structure
    cli_log("Checking file structure...", "INFO")
    tryCatch({
      # Read first few lines to check structure
      headers <- names(data.table::fread(dorado_summary, nrows = 1))
      cli_log(sprintf("Found columns: %s", paste(headers, collapse=", ")), "INFO")

      # Now read the actual data
      summary_data <- data.table::fread(dorado_summary)

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
    summary_data <- data.table::as.data.table(dorado_summary)
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
    part_data <- summary_data[start_idx:end_idx]

    # Create output filename using prefix from original file
    output_file <- file.path(save_dir,
                             sprintf("%s_part%d.txt", summary_prefix, i))

    # Save data
    data.table::fwrite(part_data,
                       file = output_file,
                       sep = "\t")

    output_files[i] <- output_file

    cli_log(sprintf("Saved part %d/%d: %s", i, num_parts, basename(output_file)), "INFO")
  }

  cli_log(sprintf("Summary split into %d parts", num_parts), "SUCCESS")

  return(output_files)
}
