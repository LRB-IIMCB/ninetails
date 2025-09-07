#' Process and optionally split dorado summary file
#'
#' This function processes a dorado alignment summary file. If the file contains more
#' reads than the specified part_size, it splits it into multiple smaller files.
#' Otherwise, it creates a copy in the summary directory. All files are saved in
#' 'dorado_summary_dir' within the specified save directory.
#'
#' The function performs the following steps:
#' 1. Validates input parameters and checks file existence
#' 2. Reads the summary file if provided as path, or uses the provided data frame
#' 3. Determines the number of parts needed based on total reads and part_size
#' 4. Creates output files with naming pattern: original_filename_partn.txt
#' 5. Splits and saves the data maintaining the original format
#'
#' If input is a data frame, uses "summary" as the default prefix for output files.
#'
#' The 'dorado_summary_dir' folder structure:
#' - If input size <= part_size: Contains a single file named 'dorado_summary.txt'
#' - If input size > part_size: Contains multiple files named 'part_XX.txt' where XX
#'   is a zero-padded number indicating the part number
#'
#'  Important notes:
#' - Output files are tab-separated text files
#' - The original file name (without extension) is used as prefix for part files
#' - Directory structure will be created if it doesn't exist
#'
#' @param dorado_summary character string or data.frame. Either full path of the .txt file
#' with dorado summary or an in-memory data.frame containing dorado summary data.
#'
#' @param save_dir character string. Full path of the directory where the processed
#' files should be stored.
#'
#' @param part_size numeric [100000] (optional). If provided, defines maximum
#' number of reads per output part.
#'
#' @param cli_log Function for logging. This function is encoded in main
#' pipeline wrapper. Its purpose is to provide neatly formatted & informative log file.
#'
#' @returns A character vector containing paths to the generated files.
#' @export
#'
#' @examples
#' \dontrun{
#' results <- process_dorado_summary(
#'   dorado_summary = "path/to/summary.txt",
#'   save_dir = "output/directory",
#'   prefix = "experiment1",
#'   part_size = 1000
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
    summary_data <- data.table::fread(dorado_summary)
    # Extract prefix from original file name
    summary_prefix <- tools::file_path_sans_ext(basename(dorado_summary))
  } else if (is.data.frame(dorado_summary)) {
    summary_data <- data.table::as.data.table(dorado_summary)
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
