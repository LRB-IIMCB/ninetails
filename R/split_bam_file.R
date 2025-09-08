#' Split BAM file into smaller parts based on read IDs from summary file
#'
#' @param bam_file Path to the BAM file to be split
#' @param dorado_summary Path to the Dorado summary file containing read IDs for filtering
#' @param part_size Number of reads per output file part (default: 100000)
#' @param save_dir Directory where split BAM files will be saved
#' @param part_number Numerical identifier for the current part being processed
#' @param cli_log Function for logging messages and progress (default: message)
#'
#' @returns Character vector containing paths to the split BAM files
#' @export
#'
#' @examples
#' \dontrun{
#' bam_files <- split_bam_file(
#'   bam_file = "path/to/aligned.bam",
#'   dorado_summary = "path/to/summary_part1.txt",
#'   save_dir = "path/to/output/",
#'   part_number = 1
#' )
#' }
split_bam_file <- function(bam_file,
                           dorado_summary,
                           part_size = 100000,
                           save_dir,
                           part_number,
                           cli_log = message) {

  # Variable binding (suppressing R CMD check from throwing an error)
  read_id <- NULL

  # Input validation
  if (missing(bam_file)) {
    stop("BAM file is missing. Please provide a valid bam_file argument.",
         call. = FALSE)
  }

  if (missing(dorado_summary)) {
    stop("Summary file is missing. Please provide a valid dorado_summary argument.",
         call. = FALSE)
  }

  if (missing(save_dir)) {
    stop("Output directory is missing. Please provide a valid save_dir argument.",
         call. = FALSE)
  }

  assertthat::assert_that(file.exists(bam_file),
                          msg = "BAM file does not exist")
  assertthat::assert_that(file.exists(dorado_summary),
                          msg = "Summary file does not exist")
  assertthat::assert_that(dir.exists(save_dir),
                          msg = "Output directory does not exist")

  input_basename <- base::basename(bam_file)
  input_noext <- base::sub("\\.bam$", "", input_basename)

  # Create output prefix using the provided part number
  output_prefix <- file.path(save_dir, sprintf("%s_part_%d", input_noext, part_number))

  # Read and process the summary file
  cli_log(sprintf("[INFO] Reading summary file: %s", basename(dorado_summary)))

  t1 <- Sys.time()
  # Get read IDs efficiently with vroom
  read_ids <- vroom::vroom(dorado_summary,
                           col_select = "read_id",
                           show_col_types = FALSE,
                           altrep = TRUE)

  total_alignments <- nrow(read_ids)
  unique_reads <- length(unique(read_ids$read_id))

  cli_log(sprintf("[INFO] Summary processed in %.2f seconds",
                  as.numeric(difftime(Sys.time(), t1, units = "secs"))))

  cli_log(sprintf("[INFO] Found %d alignments for %d unique reads",
                  total_alignments, unique_reads))

  if(total_alignments > unique_reads) {
    cli_log(sprintf("[INFO] Average alignments per read: %.2f",
                    total_alignments/unique_reads))
  }

  # Get unique read IDs and clean up
  read_ids <- unique(read_ids$read_id)
  total_reads <- length(read_ids)
  gc()

  # Calculate how many output files we'll create
  num_parts <- ceiling(total_reads / part_size)

  # Calculate number of digits needed for file numbering
  num_digits <- nchar(as.character(num_parts))

  cli_log(sprintf("[INFO] Splitting into %d files with ~%d reads each",
                  num_parts, part_size))

  # Prepare output file paths with dynamic padding
  output_files <- sprintf(paste0("%s_%0", num_digits, "d.bam"),
                          output_prefix, seq_len(num_parts))

  # Split read IDs into parts of desired size
  read_parts <- split(read_ids, ceiling(seq_along(read_ids) / part_size))

  # Set up BAM reading parameters - we only need read names for filtering
  param <- Rsamtools::ScanBamParam(what = "qname")

  # Process each part
  for (i in seq_along(read_parts)) {
    # Log progress
    cli_log(sprintf("[INFO] Processing part %d/%d (expecting %d reads)",
                    i, num_parts, length(read_parts[[i]])))

    # Prepare filter function for the current part of read IDs
    current_part <- read_parts[[i]]
    filter <- function(x) {
      matches <- x$qname %in% current_part
      return(matches)
    }

    # Filter BAM file to create part
    Rsamtools::filterBam(bam_file, output_files[i],
                         filter = S4Vectors::FilterRules(filter),
                         param = param)

    # Count actual number of reads in output file
    actual_reads <- Rsamtools::countBam(output_files[i])$records

    cli_log(sprintf("[INFO] Created: %s (contains %d alignments)",
                    basename(output_files[i]), actual_reads))

    if(actual_reads > length(current_part)) {
      cli_log(sprintf("[INFO] Note: %d additional alignments included for these reads",
                      actual_reads - length(current_part)))
    }
  }

  cli_log(sprintf("[SUCCESS] Splitting complete"))

  return(output_files)
}

