#' Initially preprocess dorado input files before passing them to the
#' main ninetails pipeline
#'
#' Preprocesses dorado input files by validating inputs and splitting large files
#' into manageable parts if needed, so they can be efficiently parsed with
#' the main ninetails pipeline without memory overload.
#'
#' @param bam_file character string. Full path to the input BAM file.
#'
#' @param dorado_summary character string or data.frame. Path to dorado summary file
#' or a data frame containing the summary data.
#'
#' @param pod5_dir character string. Directory containing pod5 files.
#'
#' @param num_cores numeric [1]. Number of cores to use for processing.
#'
#' @param qc logical [TRUE/FALSE]. Whether to perform quality control.
#'
#' @param save_dir character string. Directory where output files will be saved.
#'
#' @param prefix character string. Prefix for output files.
#'
#' @param part_size numeric [100000]. Maximum number of reads to process in each part.
#'
#' @param cli_log Function for logging. This function is encoded in main pipeline wrapper.
#'
#' @returns A list containing paths to processed files and status information
#'
#' @export
#'
#' @examples
#'\dontrun{
#'
#' # Split a BAM file based on dorado summary part
#' split_bam_file(
#'   bam_file = '/path/to/input.bam',
#'   dorado_summary = '/path/to/summary_part.txt',
#'   save_dir = '/path/to/output',
#'   part_size = 1000
#' )
#'
#'}
preprocess_dorado_input <- function(bam_file,
                                    dorado_summary,
                                    pod5_dir,
                                    num_cores,
                                    qc,
                                    save_dir,
                                    prefix,
                                    part_size,
                                    cli_log) {
  # Input validation and configuration
  cli_log("Validating input parameters", "INFO", "Validating Inputs")
  # Assertions
  ###################################################
  assertthat::assert_that(is.numeric(num_cores), num_cores > 0,
                          msg = "Number of cores must be a positive numeric value")

  assertthat::assert_that(is.character(pod5_dir), dir.exists(pod5_dir),
                          msg = "Pod5 files directory does not exist or path is invalid")

  assertthat::assert_that(is.character(save_dir),
                          msg = "Output directory path must be a character string")

  assertthat::assert_that(is.logical(qc),
                          msg = "qc must be logical [TRUE/FALSE]")

  # Validate prefix (if provided)
  if (nchar(prefix) > 0) {
    assertthat::assert_that(is.character(prefix),
                            msg = "File name prefix must be a character string")
  }

  # Validate BAM input
  if (checkmate::test_string(bam_file)) {
    assertthat::assert_that(file.exists(bam_file),
                            msg = "BAM file does not exist")
  }

  # Validate sequencing summary input
  if (checkmate::test_string(dorado_summary)) {
    assertthat::assert_that(file.exists(dorado_summary),
                            msg = "Dorado summary file does not exist")
  } else {
    assertthat::assert_that(is.data.frame(dorado_summary) && nrow(dorado_summary) > 0,
                            msg = "Dorado summary must be a non-empty data frame or valid file path")
  }
  cli_log("Provided input files/paths are in correct format", "SUCCESS")

  # Process dorado summary
  cli_log("Processing dorado summary...", "INFO", "Dorado Summary Processing", bullet = TRUE)

  # Process summary file - function handles size checking internally
  part_files <- process_dorado_summary(
    dorado_summary = dorado_summary,
    save_dir = save_dir,
    part_size = part_size,
    cli_log = cli_log
  )

  # Create BAM directory
  bam_dir <- file.path(save_dir, "bam_dir")
  if (!dir.exists(bam_dir)) {
    dir.create(bam_dir, recursive = TRUE)
  }

  # Process BAM file based on whether summary was split
  cli_log("Processing BAM file...", "INFO", "BAM Processing", bullet = TRUE)

  if (length(part_files) > 1) {
    cli_log(sprintf("Splitting BAM file to match %d summary parts", length(part_files)), "INFO")

    bam_files <- vector("character", length(part_files))

    for (i in seq_along(part_files)) {
      cli_log(sprintf("Processing BAM part %d/%d", i, length(part_files)), "INFO")

      bam_files[i] <- split_bam_file(
        bam_file = bam_file,
        dorado_summary = part_files[i],
        part_size = part_size,
        save_dir = bam_dir
      )
    }

    cli_log(sprintf("BAM file split into %d parts", length(bam_files)), "SUCCESS")

  } else {
    # Just copy the original BAM file
    new_bam_path <- file.path(bam_dir, basename(bam_file))
    file.copy(from = bam_file, to = new_bam_path)
    bam_files <- new_bam_path

    cli_log("BAM file copied as single file", "SUCCESS")
  }

  cli_log("All input files processed and ready", "SUCCESS")

  # Return processed file paths
  return(list(
    summary_files = part_files,
    bam_files = bam_files
  ))
}
