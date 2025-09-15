#' Preprocess Dorado input files for poly(A) analysis
#'
#' Validates and splits input files, extracts poly(A) information and signals, and prepares all files for downstream analysis.
#'
#' @param bam_file Character. Path to the BAM file containing aligned reads.
#' @param dorado_summary Character or data frame. Path to Dorado summary file or data frame.
#' @param pod5_dir Character. Directory containing pod5 files.
#' @param num_cores Integer. Number of CPU cores to use.
#' @param qc Logical. Whether to perform quality control.
#' @param save_dir Character. Directory where output files will be saved.
#' @param prefix Character. Prefix to add to output file names (optional).
#' @param part_size Integer. Number of reads to process in each chunk.
#' @param cli_log Function for logging messages and progress.
#'
#' @return List containing paths to processed files:
#'   \describe{
#'     \item{summary_files}{Paths to split summary files}
#'     \item{bam_files}{Paths to split BAM files (if BAM processing performed)}
#'     \item{polya_files}{Paths to extracted poly(A) data files}
#'     \item{polya_signal_files}{Paths to extracted poly(A) signal files}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' processed_files <- preprocess_inputs(
#'   bam_file = "path/to/aligned.bam",
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
preprocess_inputs <- function(bam_file,
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
  assertthat::assert_that(is.numeric(num_cores), num_cores > 0,
                          msg = "Number of cores must be a positive numeric value")
  assertthat::assert_that(is.character(pod5_dir), dir.exists(pod5_dir),
                          msg = "Pod5 files directory does not exist or path is invalid")
  assertthat::assert_that(is.character(save_dir),
                          msg = "Output directory path must be a character string")
  assertthat::assert_that(is.logical(qc),
                          msg = "qc must be logical [TRUE/FALSE]")
  if (nchar(prefix) > 0) {
    assertthat::assert_that(is.character(prefix),
                            msg = "File name prefix must be a character string")
  }
  if (checkmate::test_string(bam_file)) {
    assertthat::assert_that(file.exists(bam_file),
                            msg = "BAM file does not exist")
  }
  if (checkmate::test_string(dorado_summary)) {
    assertthat::assert_that(file.exists(dorado_summary),
                            msg = "Dorado summary file does not exist")
  } else {
    assertthat::assert_that(is.data.frame(dorado_summary) && nrow(dorado_summary) > 0,
                            msg = "Dorado summary must be a non-empty data frame or valid file path")
  }
  cli_log("Provided input files/paths are in correct format", "SUCCESS")

  ################################################################################
  # CHECK FOR POLYA COLUMNS IN SUMMARY
  ################################################################################
  if (is.character(dorado_summary)) {
    summary_preview <- vroom::vroom(dorado_summary, n_max = 1, show_col_types = FALSE)
    summary_cols <- colnames(summary_preview)
  } else if (is.data.frame(dorado_summary)) {
    summary_cols <- colnames(dorado_summary)
  } else {
    stop("dorado_summary must be a file path or a data frame")
  }
  has_polya_cols <- all(c("poly_tail_length", "poly_tail_start", "polya_end") %in% summary_cols)

  ################################################################################
  # DORADO SUMMARY PROCESSING
  ################################################################################
  cli_log("Processing dorado summary...", "INFO", "Dorado Summary Processing", bullet = TRUE)
  summary_dir <- file.path(save_dir, "dorado_summary_dir")
  if (!dir.exists(summary_dir)) {
    dir.create(summary_dir, recursive = TRUE)
  }
  part_files <- process_dorado_summary(
    dorado_summary = dorado_summary,
    save_dir = summary_dir,
    part_size = part_size,
    cli_log = cli_log
  )
  if (length(part_files) == 0) {
    stop("No summary parts were created")
  }

  ################################################################################
  # ROUTE SELECTION: DIRECT SUMMARY OR BAM PROCESSING
  ################################################################################
  if (has_polya_cols) {
    cli_log("Summary contains poly(A) columns. Skipping BAM processing.", "INFO", "Route Selection", bullet = TRUE)
    bam_files <- NULL
    polya_files <- part_files
  } else {
    cli_log("Summary does not contain poly(A) columns. Performing BAM processing.", "INFO", "Route Selection", bullet = TRUE)
    ################################################################################
    # BAM PRE-PROCESSING
    ################################################################################
    bam_dir <- file.path(save_dir, "bam_dir")
    if (!dir.exists(bam_dir)) {
      dir.create(bam_dir, recursive = TRUE)
    }
    cli_log("Processing BAM file...", "INFO", "BAM Processing", bullet = TRUE)
    if (length(part_files) > 1) {
      cli_log(sprintf("Splitting BAM file to match %d summary parts", length(part_files)), "INFO")
      bam_files <- vector("character", length(part_files))
      for (i in seq_along(part_files)) {
        cli_log(sprintf("Processing BAM part %d/%d", i, length(part_files)), "INFO")
        current_summary <- part_files[i]
        if (!file.exists(current_summary)) {
          stop(sprintf("Summary file does not exist: %s", current_summary))
        }
        bam_files[i] <- split_bam_file(
          bam_file = bam_file,
          dorado_summary = current_summary,
          part_size = part_size,
          save_dir = bam_dir,
          part_number = i,
          cli_log = cli_log
        )
      }
      cli_log(sprintf("BAM file split into %d parts", length(bam_files)), "SUCCESS")
    } else {
      new_bam_path <- file.path(bam_dir, basename(bam_file))
      file.copy(from = bam_file, to = new_bam_path)
      bam_files <- new_bam_path
      cli_log("BAM file copied as single file", "SUCCESS")
    }
    ################################################################################
    # EXTRACTING POLYA FROM BAM
    ################################################################################
    cli_log("Extracting poly(A) info from BAM file...", "INFO", "Poly(A) info Extracting", bullet = TRUE)
    polya_dir <- file.path(save_dir, "polya_data_dir")
    if (!dir.exists(polya_dir)) {
      dir.create(polya_dir, recursive = TRUE)
    }
    polya_files <- character(length(bam_files))
    for (i in seq_along(bam_files)) {
      current_bam <- bam_files[i]
      current_summary <- part_files[i]
      cli_log(sprintf("Processing BAM file %d/%d for poly(A) extraction", i, length(bam_files)), "INFO")
      polya_data <- extract_polya_from_bam(
        bam_file = current_bam,
        summary_file = current_summary,
        cli_log = cli_log
      )
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
      gc()
    }
    polya_files <- polya_files[file.exists(polya_files)]
    if (length(polya_files) > 0) {
      cli_log(sprintf("Successfully processed %d BAM files", length(polya_files)), "SUCCESS")
    } else {
      cli_log("No poly(A) data was extracted from any BAM file", "WARNING")
    }
  }

  ################################################################################
  # EXTRACTING POLYA SIGNALS FROM POD5
  ################################################################################
  cli_log("Signal extracting", "INFO", "Signal Extracting")
  cli_log("Extracting poly(A) signals from pod5 file...", "INFO", bullet = TRUE)
  polya_signal_dir <- file.path(save_dir, "polya_signal_dir")
  if (!dir.exists(polya_signal_dir)) {
    dir.create(polya_signal_dir, recursive = TRUE)
  }
  polya_signal_files <- character(length(polya_files))
  for (i in seq_along(polya_files)) {
    polya_data <- vroom::vroom(polya_files[i], delim = "\t", show_col_types = FALSE)
    signal_list <- extract_tails_from_pod5(polya_data, pod5_dir)
    output_file <- file.path(polya_signal_dir,
                             sprintf("%spolya_signal_part%d.rds",
                                     ifelse(nchar(prefix) > 0, paste0(prefix, "_"), ""),
                                     i))
    saveRDS(signal_list, output_file)
    polya_signal_files[i] <- output_file
    cli_log(sprintf("Saved extracted signals for poly(A) data part %d to %s", i, output_file), "INFO")
  }

  ################################################################################
  # RETURN OUTPUT FILE PATHS
  ################################################################################
  return(list(
    summary_files = part_files,
    bam_files = bam_files,
    polya_files = polya_files,
    polya_signal_files = polya_signal_files
  ))
}
