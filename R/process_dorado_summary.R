#' Process and optionally split dorado summary file
#'
#' @description
#' This function processes a dorado alignment summary file. If the file contains more
#' reads than the specified part_size, it splits it into multiple smaller files.
#' Otherwise, it creates a copy in the summary directory. All files are saved in
#' 'dorado_summary_dir' within the specified save directory.
#'
#' @details
#' The 'dorado_summary_dir' folder structure:
#' - If input size <= part_size: Contains a single file named 'dorado_summary.txt'
#' - If input size > part_size: Contains multiple files named 'part_XX.txt' where XX
#'   is a zero-padded number indicating the part number
#'
#' @param dorado_summary character string or data.frame. Either full path of the .txt file
#' with dorado summary or an in-memory data.frame containing dorado summary data.
#' @param save_dir character string. Full path of the directory where the processed
#' files should be stored.
#' @param part_size numeric [100000] (optional). If provided, defines maximum
#' number of reads per output part.
#' @param cli_log Function for logging. This function is encoded in main
#' pipeline wrapper. Its purpose is to provide neatly formatted & informative log file.
#'
#' @return A character vector containing paths to the generated files.
#'
#' @importFrom assertthat assert_that
#' @importFrom checkmate test_string assert_file_exists
#' @importFrom vroom vroom vroom_write
process_dorado_summary <- function(dorado_summary,
                                   save_dir,
                                   part_size = 100000,
                                   cli_log = message) {

  # Variable binding (suppressing R CMD check from throwing an error)
  filename <- read_id <- NULL

  # Assertions
  if (missing(dorado_summary)) {
    stop("Dorado summary is missing. Please provide a valid dorado_summary argument.",
         call. = FALSE)
  }

  if (missing(save_dir)) {
    stop("Output directory is missing. Please provide a valid save_dir argument.",
         call. = FALSE)
  }

  assertthat::assert_that(is.numeric(part_size),
                          msg = "part_size must be numeric")
  assertthat::assert_that(part_size > 0,
                          msg = "part_size must be greater than 0")

  cli_log("[INFO] Processing dorado summary")

  # Create output directory if it doesn't exist
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }

  # Create summary directory
  summary_dir <- file.path(save_dir, "dorado_summary_dir")
  if (!dir.exists(summary_dir)) {
    dir.create(summary_dir, recursive = TRUE)
  }

  # Load dorado summary data
  if (checkmate::test_string(dorado_summary)) {
    # If string provided as an argument, read from file
    checkmate::assert_file_exists(dorado_summary)
    dorado_data <- vroom::vroom(dorado_summary,
                                show_col_types = FALSE)
    cli_log(sprintf("[INFO] Read summary file: %s", basename(dorado_summary)))
  } else {
    # Make sure dorado_summary is an object with rows
    if (!is.data.frame(dorado_summary) || nrow(dorado_summary) == 0) {
      stop("Empty data frame provided as input (dorado_summary). Please provide valid input")
    }
    cli_log("[INFO] Processing in-memory summary data")
    dorado_data <- dorado_summary
  }

  # Get total number of reads
  total_reads <- nrow(dorado_data)
  cli_log(sprintf("[INFO] Summary contains %d reads", total_reads))

  # If total reads is less than or equal to part_size, just copy the file
  if (total_reads <= part_size) {
    output_file <- file.path(summary_dir, "dorado_summary.txt")
    vroom::vroom_write(dorado_data,
                       output_file,
                       delim = "\t")
    cli_log(sprintf("[INFO] Summary contains %d reads, no splitting required", total_reads))
    cli_log(sprintf("[INFO] Saved to: %s", basename(output_file)))
    return(output_file)
  } else {
    # Calculate number of parts needed and pad width for naming
    num_parts <- ceiling(total_reads / part_size)
    pad_width <- nchar(as.character(num_parts))

    cli_log(sprintf("[INFO] Splitting into %d parts (%d reads per part)", num_parts, part_size))

    # Split the data and write files
    output_files <- character(num_parts)
    for (i in seq_len(num_parts)) {
      start_idx <- ((i - 1) * part_size) + 1
      end_idx <- min(i * part_size, total_reads)

      # Create padded file number
      file_num <- sprintf(paste0("%0", pad_width, "d"), i)
      output_file <- file.path(summary_dir,
                               sprintf("part_%s.txt", file_num))

      # Write the part to file
      vroom::vroom_write(dorado_data[start_idx:end_idx, ],
                         output_file,
                         delim = "\t")

      output_files[i] <- output_file

      cli_log(sprintf("[INFO] Created part %d/%d: %s",
                      i, num_parts, basename(output_file)))
    }

    cli_log(sprintf("[SUCCESS] Generated %d summary parts in: %s",
                    num_parts, basename(summary_dir)))

    return(output_files)
  }
}
