#' Preprocess Dorado input files for poly(A) analysis
#'
#' @param bam_file Path to the BAM file containing aligned reads
#' @param dorado_summary Path to Dorado summary file or a data frame containing summary information
#' @param pod5_dir Directory containing pod5 files
#' @param num_cores Number of CPU cores to use for processing
#' @param qc Logical indicating whether to perform quality control checks
#' @param save_dir Directory where output files will be saved
#' @param prefix Prefix to add to output file names (optional)
#' @param part_size Number of reads to process in each chunk
#' @param cli_log Function for logging messages and progress
#'
#' @returns A list containing paths to processed files:
#'   \itemize{
#'     \item summary_files - Paths to split summary files
#'     \item bam_files - Paths to split BAM files
#'     \item polya_files - Paths to extracted poly(A) data files
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' processed_files <- preprocess_dorado_input(
#'   bam_file = "path/to/aligned.bam",
#'   dorado_summary = "path/to/summary.txt",
#'   pod5_dir = "path/to/pod5/",
#'   num_cores = 2,
#'   qc = TRUE,
#'   save_dir = "path/to/output/",
#'   prefix = "experiment1",
#'   part_size = 40000
#' )
#' }
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

  ###################################################
  # DORADO SUMMARY PROCESSING
  ###################################################

  # Process dorado summary
  cli_log("Processing dorado summary...", "INFO", "Dorado Summary Processing", bullet = TRUE)

  # Create dorado summary directory (single level)
  summary_dir <- file.path(save_dir, "dorado_summary_dir")
  if (!dir.exists(summary_dir)) {
    dir.create(summary_dir, recursive = TRUE)
  }

  # Process summary file - handles size checking internally
  part_files <- ninetails::process_dorado_summary(
    dorado_summary = dorado_summary,
    save_dir = summary_dir,  # Now saving directly to summary_dir
    part_size = part_size,
    cli_log = cli_log
  )

  # Verify summary parts were created
  if (length(part_files) == 0) {
    stop("No summary parts were created")
  }

  ###################################################
  # BAM PRE-PROCESSING
  ###################################################
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

      current_summary <- part_files[i]
      # Verify summary file exists before proceeding
      if (!file.exists(current_summary)) {
        stop(sprintf("Summary file does not exist: %s", current_summary))
      }

      bam_files[i] <- split_bam_file(
        bam_file = bam_file,
        dorado_summary = current_summary,
        part_size = part_size,
        save_dir = bam_dir,
        part_number = i,  # Added part number
        cli_log = cli_log  # Added cli_log
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

  ###################################################
  # EXTRACTING POLYA FROM BAM & SUMMARY
  ###################################################

  # Inform that poly(A) will be extracted
  cli_log("Extracting poly(A) info from BAM file...", "INFO", "Poly(A) info Extracting", bullet = TRUE)
  # Create polya directory
  polya_dir <- file.path(save_dir, "polya_data_dir")
  if (!dir.exists(polya_dir)) {
    dir.create(polya_dir, recursive = TRUE)
  }

  # Process BAM files sequentially with their corresponding summary parts
  polya_files <- character(length(bam_files))

  for (i in seq_along(bam_files)) {
    current_bam <- bam_files[i]
    current_summary <- part_files[i]

    cli_log(sprintf("Processing BAM file %d/%d for poly(A) extraction", i, length(bam_files)), "INFO")

    # Extract poly(A) data from current BAM file with corresponding summary
    polya_data <- ninetails::extract_polya_from_bam(
      bam_file = current_bam,
      summary_file = current_summary,
      cli_log = cli_log  # Added cli_log
    )

    # Save results if we have any
    if (nrow(polya_data) > 0) {
      output_file <- file.path(polya_dir,
                               sprintf("%spolya_data_part%d.tsv",
                                       ifelse(nchar(prefix) > 0, paste0(prefix, "_"), ""),
                                       i))

      vroom::vroom_write(polya_data, output_file, delim = "\t")
      polya_files[i] <- output_file
      cli_log(sprintf("Saved poly(A) data for BAM %d to %s (%d reads)",
                      i, output_file, nrow(polya_data)), "INFO")
    } else {
      cli_log(sprintf("No poly(A) data extracted from BAM %d", i), "INFO")
    }

    # Explicitly trigger garbage collection after processing each BAM file
    gc()
  }

  # Filter out any empty results
  polya_files <- polya_files[file.exists(polya_files)]

  if (length(polya_files) > 0) {
    cli_log(sprintf("Successfully processed %d BAM files", length(polya_files)), "SUCCESS")
  } else {
    cli_log("No poly(A) data was extracted from any BAM file", "WARNING")
  }

  # Return processed file paths
  return(list(
    summary_files = part_files,
    bam_files = bam_files,
    polya_files = polya_files
  ))
}
