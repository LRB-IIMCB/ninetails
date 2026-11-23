################################################################################
# WRAPPER FUNCTION FOR DORADO cDNA OUTPUTS
################################################################################
#' Complete Oxford Nanopore poly(A)/poly(T) tail analysis pipeline for Dorado cDNA data.
#'
#' This function extends the check_tails_dorado_DRS pipeline to handle cDNA sequencing data
#' by adding BAM file processing for sequence extraction and Dorado-style read orientation
#' classification. The pipeline identifies and characterizes non-adenosine nucleotides within
#' poly(A) and poly(T) tails using the same signal processing and machine learning approach
#' as the DRS pipeline, but with automatic classification of read orientations.
#'
#' @section Pipeline Overview:
#' \enumerate{
#'   \item \strong{Tail coordinates and orientation}: taken from Lachesis (\code{lachesis} param),
#'   or from Dorado summary with orientation classified from the basecalled sequences in the BAM file
#'   (based on SSP/VNP primer matching)
#'   \item \strong{Signal extraction}: raw tail signals with flanks are extracted from POD5 files
#'   \item \strong{Classification}: sliding windows of each tail are classified with the poly(A)
#'   or poly(T) cDNA model (see \code{\link{predict_nonA_cdna}})
#'   \item \strong{Output with Tail Types}: Produces standard read_classes and nonadenosine_residues with tail_type information
#' }
#'
#' @section Input Requirements:
#' This pipeline requires specific input formats:
#' \itemize{
#'   \item \strong{Lachesis output}: tail coordinates and read orientation from Lachesis software (recommended) or
#'   \item \strong{Dorado Summary}: Must contain standard columns for read information (\code{--estimate-poly-a}; Dorado >=1.0.0)
#'   \item \strong{BAM File}: Aligned cDNA reads with basecalled, untrimmed sequences (\code{--no-trim})
#'   \item \strong{POD5 Files}: Raw signal files corresponding to reads in summary
#' }
#'
#' @section Key Differences from DRS Pipeline:
#' \itemize{
#'   \item \strong{BAM Input}: Additional BAM file input for sequence extraction since cDNA data requires basecalled sequences
#'   \item \strong{Dorado-Style Read Orientation Classification}: matching of SSP/VNP primers to classify reads as polyA, polyT, or unidentified before processing
#'   \item \strong{Separate models}: poly(A) and poly(T) reads are classified with dedicated cDNA models working on the raw (not downsampled, nor winsorized) signal
#'   \item \strong{RNA-level reporting}: predictions refer to the residues of the original RNA tail in both orientations
#' }
#'
#' @param bam_file Character string. Path to the BAM file containing aligned cDNA reads
#' with basecalled sequences. This file will be split into parts for memory management.
#'
#' @param dorado_summary Character string or data frame. Path to Dorado summary file or
#' data frame containing per-read summary information. Must include standard columns
#' such as read_id, filename, etc.
#'
#' @param pod5_dir Character string. Path to directory containing POD5 files with raw
#' nanopore signal data corresponding to the reads in the summary file.
#'
#' @param num_cores Integer [1]. Number of CPU cores to use for parallel processing.
#' Recommend using `parallel::detectCores() - 1` for optimal performance while
#' maintaining system responsiveness.
#'
#' @param qc Logical [TRUE]. Whether to apply quality control filtering during
#' analysis. When TRUE, applies standard ninetails QC filters including tail length
#' filtering and coordinate validation.
#'
#' @param save_dir Character string. Path to directory where all output files and
#' intermediate results will be saved. Directory will be created if it doesn't exist.
#'
#' @param prefix Character string [""]. Optional prefix to add to all output file names.
#' Useful for distinguishing between different experimental conditions or samples.
#'
#' @param part_size Integer [40000]. Number of reads to process in each chunk when
#' splitting large input files. Larger values use more memory but may be faster.
#' Adjust based on available system memory.
#'
#' @param cleanup Logical [FALSE]. Whether to remove intermediate files after successful
#' pipeline completion. When FALSE, all intermediate files are preserved for inspection.
#'
#' @param lachesis Character string or data frame [NULL]. Path to the Lachesis
#' output file (or the same table in memory). When provided, the tail coordinates
#' and read orientation (\code{strand} column) are based on Lachesis predictions.
#' If that is the case, \code{bam_file} and \code{dorado_summary} are not required.
#'
#' @return A named list containing the final analysis results:
#' \describe{
#'   \item{read_classes}{Data frame with per-read classification results including
#'     readname, contig, poly(A) length, QC tag, class, comments, and tail_type}
#'   \item{nonadenosine_residues}{Data frame with predicted non-adenosine positions
#'     within poly(A)/poly(T) tails including readname, contig, prediction,
#'     estimated position, poly(A) length, QC tag, and tail_type. Predictions
#'     refer to the residue present in the original RNA tail (C, G or U) and
#'     positions are counted from the 3' end of the tail, in both orientations}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic cDNA analysis
#' results <- ninetails::check_tails_dorado_cDNA(
#'   bam_file = "path/to/aligned_cdna.bam",
#'   dorado_summary = "path/to/dorado_summary.txt",
#'   pod5_dir = "path/to/pod5_files/",
#'   num_cores = 4,
#'   save_dir = "path/to/output/"
#' )
#'
#' # Access results
#' head(results$read_classes)
#' head(results$nonadenosine_residues)
#'
#' # Analysis with custom settings
#' results <- ninetails::check_tails_dorado_cDNA(
#'   bam_file = "large_dataset.bam",
#'   dorado_summary = summary_df,  # Can pass data frame
#'   pod5_dir = "/data/pod5/",
#'   num_cores = 8,
#'   qc = TRUE,
#'   save_dir = "/results/experiment1/",
#'   prefix = "exp1_sample_A",
#'   part_size = 20000,  # Smaller chunks for limited memory
#'   cleanup = TRUE,  # Remove intermediate files
#'   lachesis = NULL # use Dorado input
#' )
#' }
check_tails_dorado_cDNA <- function(bam_file,
                                    dorado_summary,
                                    pod5_dir,
                                    num_cores = 1,
                                    qc = TRUE,
                                    save_dir,
                                    prefix = "",
                                    part_size = 40000,
                                    cleanup = FALSE,
                                    lachesis = NULL) {

  # Initialize warning flag
  warn_message <- FALSE


  ################################################################################
  # SETUP AND VALIDATION
  ################################################################################

  # Create output directory if needed
  tryCatch({
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }
  }, error = function(e) {
    cli::cli_alert_danger("Failed to create output directory: {e$message}")
    stop(e$message, call. = FALSE)
  })

  # Initialize log collector and setup logging
  log_collector <- vector("character")
  log_start_time <- Sys.time()
  log_filename <- format(log_start_time,
                         paste0("%Y-%m-%d_%H-%M-%S",
                                if(nchar(prefix) > 0) paste0("_", prefix) else "",
                                "_ninetails_cDNA.log"))
  log_filepath <- file.path(save_dir, log_filename)

  # Helper functions for logging (same as DRS pipeline)
  log_message <- function(message, type = "INFO", section = NULL) {
    timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")

    if (!is.null(section)) {
      header <- paste0("\n-- ", section, " --\n")
      cat(header, file = log_filepath, append = TRUE)
      log_collector <<- c(log_collector, header)
    }

    if (grepl("\n", message)) {
      lines <- unlist(strsplit(message, "\n"))
      for (line in lines) {
        if (nchar(trimws(line)) > 0) {
          log_entry <- sprintf("%s %s: %s", timestamp, type, line)
          log_collector <<- c(log_collector, log_entry)
          cat(paste0(log_entry, "\n"), file = log_filepath, append = TRUE)
        }
      }
    } else {
      log_entry <- sprintf("%s %s: %s", timestamp, type, message)
      log_collector <<- c(log_collector, log_entry)
      cat(paste0(log_entry, "\n"), file = log_filepath, append = TRUE)
    }
  }

  cli_log <- function(message, type = "INFO", section = NULL, bullet = FALSE) {
    if (!is.null(section)) {
      cli::cli_h2(section)
    }

    if (bullet) {
      cli::cli_bullets(c("*" = message))
    } else {
      switch(type,
             "INFO" = cli::cli_alert_info(message),
             "SUCCESS" = cli::cli_alert_success(message),
             "WARNING" = cli::cli_alert_warning(message),
             "ERROR" = cli::cli_alert_danger(message))
    }
    log_message(message, type, section)
  }

  ################################################################################
  # PIPELINE START
  ################################################################################

  tryCatch({
    # Initialize pipeline
    ###################################################
    cli::cli_h1("Ninetails cDNA Pipeline {utils::packageVersion('ninetails')}")
    cli_log(paste0("Welcome to Ninetails cDNA Pipeline! v", utils::packageVersion('ninetails')), "INFO")
    cli::cli_alert_info("Analysis toolkit for poly(A)/poly(T) tail composition in ONT cDNA data")

    # Write log file header with welcome message
    cat(sprintf("Ninetails cDNA Analysis Log - %s\n", format(log_start_time, "%Y-%m-%d %H:%M:%S")),
        file = log_filepath, append = FALSE)
    cat(sprintf("Ninetails v%s - cDNA Analysis for poly(A)/poly(T) tail composition\n",
                utils::packageVersion('ninetails')),
        file = log_filepath, append = TRUE)
    cat("=====================================\n\n", file = log_filepath, append = TRUE)

    # Write any messages from directory checking to the log file
    if (length(log_collector) > 0) {
      for (entry in log_collector) {
        cat(paste0(entry, "\n"), file = log_filepath, append = TRUE)
      }
    }

    cli_log("Launching *check_tails_dorado_cDNA* pipeline", "INFO","Launching cDNA Pipeline")

    # Configuration logging display
    #####################################################
    cli_log("Configuration", "INFO", "Configuration")
    if (!is.null(lachesis)) {
      cli_log(sprintf("Lachesis output:             %s",
                      if(!is.object(lachesis)) lachesis else deparse(substitute(lachesis))), bullet = TRUE)
    } else {
      cli_log(sprintf("BAM file:                    %s", bam_file), bullet = TRUE)
      cli_log(sprintf("Dorado summary:              %s",
                      if(!is.object(dorado_summary)) dorado_summary else deparse(substitute(dorado_summary))), bullet = TRUE)
    }
    cli_log(sprintf("Pod5 files directory:        %s", pod5_dir), bullet = TRUE)
    cli_log(sprintf("Number of cores:             %d", num_cores), bullet = TRUE)
    cli_log(sprintf("Output quality control:      %s", qc), bullet = TRUE)
    cli_log(sprintf("Output directory:            %s", save_dir), bullet = TRUE)
    cli_log(sprintf("Reads processed at once:     %s", part_size), bullet = TRUE)
    cli_log(sprintf("Cleanup intermediate files:  %s", cleanup), bullet = TRUE)

    ################################################################################
    # INPUT PREPROCESSING AND BAM PROCESSING
    ################################################################################

    if (!is.null(lachesis)) {

      ################################################################################
      # LACHESIS INPUT: COORDINATES AND READ ORIENTATION FROM LACHESIS
      ################################################################################

      assert_condition(is.numeric(num_cores) && num_cores > 0,
                       "Number of cores must be a positive numeric value")
      assert_condition(is.numeric(part_size) && part_size > 0,
                       "Part size must be a positive numeric value")
      assert_condition(is.logical(qc),
                       "qc must be logical [TRUE/FALSE]")
      assert_dir_exists(pod5_dir, "Pod5 files")

      cli_log("Reading tail coordinates from Lachesis output...", "INFO", "Tail Coordinates")
      tail_coordinates <- ninetails::read_tail_coordinates_cdna(lachesis, source = "lachesis")

      # Extract chunked raw tail signals (memory management)
      polya_signal_dir <- file.path(save_dir, "polya_signal_dir")
      if (!dir.exists(polya_signal_dir)) {
        dir.create(polya_signal_dir, recursive = TRUE)
      }

      read_parts <- split(seq_len(nrow(tail_coordinates)),
                          ceiling(seq_len(nrow(tail_coordinates)) / part_size))
      signal_files <- character(0)

      for (i in seq_along(read_parts)) {
        cli_log(sprintf("Extracting signals from POD5 for part %d/%d", i, length(read_parts)), "INFO")

        signal_list <- ninetails::extract_tail_signals_cdna(
          tail_coordinates = tail_coordinates[read_parts[[i]], ],
          pod5_dir = pod5_dir,
          num_cores = num_cores
        )

        if (length(signal_list) > 0) {
          output_file <- file.path(polya_signal_dir,
                                   sprintf("%spolya_signal_part%d.rds",
                                           ifelse(nchar(prefix) > 0, paste0(prefix, "_"), ""),
                                           i))
          saveRDS(signal_list, output_file)
          signal_files <- c(signal_files, output_file)
        }

        rm(signal_list)
        gc()
      }

    } else {

      # Preprocess input files - BAM processing and sequence extraction
      processed_files <- ninetails::preprocess_inputs_cdna(
        bam_file = bam_file,
        dorado_summary = dorado_summary,
        pod5_dir = pod5_dir,
        num_cores = num_cores,
        qc = qc,
        save_dir = save_dir,
        prefix = prefix,
        part_size = part_size,
        cli_log = cli_log
      )

      ################################################################################
      # SEQUENCE CLASSIFICATION (polyA/polyT/unidentified)
      ################################################################################

      cli_log("Starting Dorado-style poly tail classification...", "INFO", "Sequence Classification")
      cli_log("Using SSP/VNP primer edit distance matching with score and separation validation", "INFO")

      # Classify read orientations and add tail_type column
      classified_sequences <- ninetails::detect_orientation_multiple(
        sequence_files = processed_files$sequence_files,
        num_cores = num_cores,
        cli_log = cli_log
      )

      cli_log(sprintf("Classification completed: %d total sequences with tail_type column",
                      nrow(classified_sequences)), "SUCCESS")

      # Tail coordinates from the Dorado summary, orientation from the classification
      tail_coordinates <- do.call(rbind, lapply(processed_files$summary_files, function(x) {
        ninetails::read_tail_coordinates_cdna(x, source = "dorado")
      }))
      tail_coordinates$tail_type <- classified_sequences$tail_type[
        match(tail_coordinates$read_id, classified_sequences$read_id)]
      tail_coordinates$tail_type[is.na(tail_coordinates$tail_type)] <- "unidentified"

      signal_files <- processed_files$polya_signal_files
    }

    ################################################################################
    # DUAL PROCESSING: POLYA AND POLYT (SEPARATE MODELS)
    ################################################################################

    # Split reads by tail type
    polya_reads <- tail_coordinates[tail_coordinates$tail_type == "polyA", ]
    polyt_reads <- tail_coordinates[tail_coordinates$tail_type == "polyT", ]
    unidentified_reads <- tail_coordinates[tail_coordinates$tail_type == "unidentified", ]

    cli_log(sprintf("Read distribution: %d polyA, %d polyT, %d unidentified",
                    nrow(polya_reads), nrow(polyt_reads), nrow(unidentified_reads)), "INFO")

    # Note: Both processing functions use the same signal files since signal extraction
    # preceeds split by orientation. Processing functions filter signals internally
    # based on read IDs in their input tables.

    # Process poly(A) reads with the cDNA poly(A) model
    polya_results <- NULL
    if (nrow(polya_reads) > 0 && length(signal_files) > 0) {
      cli_log("Processing polyA reads...", "INFO", "PolyA Processing", bullet = TRUE)
      polya_results <- ninetails::process_polya_reads_cdna(
        polya_sequences = polya_reads,
        signal_files = signal_files,  # All signal files (filtering in processing function)
        num_cores = num_cores,
        qc = qc,
        save_dir = save_dir,
        prefix = prefix,
        cli_log = cli_log
      )
    }

    # Process poly(T) reads with the cDNA poly(T) model
    polyt_results <- NULL
    if (nrow(polyt_reads) > 0 && length(signal_files) > 0) {
      cli_log("Processing polyT reads...", "INFO", "PolyT Processing", bullet = TRUE)
      polyt_results <- ninetails::process_polyt_reads_cdna(
        polyt_sequences = polyt_reads,
        signal_files = signal_files,  # All signal files (filtering in processing function)
        num_cores = num_cores,
        qc = qc,
        save_dir = save_dir,
        prefix = prefix,
        cli_log = cli_log
      )
    }

    ################################################################################
    # MERGE RESULTS AND CREATE FINAL OUTPUT
    ################################################################################

    # Merge polyA and polyT results into standard format (same as DRS pipeline)
    cli_log("Merging results and creating final output...", "INFO", "Creating Final Output", bullet = TRUE)
    outputs <- ninetails::merge_cdna_results(
      polya_results = polya_results,
      polyt_results = polyt_results,
      unidentified_reads = unidentified_reads,
      save_dir = save_dir,
      prefix = prefix,
      cli_log = cli_log
    )

    # Validate output format - a list with read_classes and nonadenosine_residues (same as DRS)
    if (!is.list(outputs)) {
      stop("merge_cdna_results must return a list", call. = FALSE)
    }

    # Check for expected components
    expected_names <- c("read_classes", "nonadenosine_residues")
    if (!all(expected_names %in% names(outputs))) {
      cli_log("Warning: Output does not contain expected components. Found components:", "WARNING")
      cli_log(paste(names(outputs), collapse = ", "), "WARNING")
    }

    # Save final outputs using standard Ninetails function (same as DRS pipeline)
    cli_log("Saving final outputs...", "INFO", "Saving Results", bullet = TRUE)
    ninetails::save_outputs(outputs, save_dir, prefix)

    ################################################################################
    # CLEANUP AND FINALIZATION
    ################################################################################

    # Cleanup intermediate files if requested
    if (cleanup) {
      cli_log("Cleaning up intermediate files...", "INFO", "Cleanup", bullet = TRUE)

      # List of intermediate directories to remove (same pattern as DRS)
      intermediate_dirs <- c(
        file.path(save_dir, "bam_parts_dir"),
        file.path(save_dir, "sequence_files_dir"),
        file.path(save_dir, "polya_signal_dir"),
        file.path(save_dir, "polya_processing_dir"),
        file.path(save_dir, "polyt_processing_dir"),
        file.path(save_dir, "polya_temp_dir"),
        file.path(save_dir, "polyt_temp_dir"),
        file.path(save_dir, "polya_chunks_dir"),
        file.path(save_dir, "polyt_chunks_dir")
      )

      cleanup_summary <- data.frame(
        directory = character(0),
        status = character(0),
        files_removed = integer(0),
        stringsAsFactors = FALSE
      )

      for (dir_path in intermediate_dirs) {
        if (dir.exists(dir_path)) {
          tryCatch({
            files_count <- length(list.files(dir_path, recursive = TRUE))
            unlink(dir_path, recursive = TRUE)

            if (!dir.exists(dir_path)) {
              cli_log(sprintf("Removed: %s (%d files)", basename(dir_path), files_count), "SUCCESS")
              cleanup_summary <- rbind(cleanup_summary,
                                       data.frame(directory = basename(dir_path),
                                                  status = "removed",
                                                  files_removed = files_count,
                                                  stringsAsFactors = FALSE))
            }
          }, error = function(e) {
            cli_log(sprintf("Failed to remove %s: %s", basename(dir_path), e$message), "WARNING")
            cleanup_summary <- rbind(cleanup_summary,
                                     data.frame(directory = basename(dir_path),
                                                status = "failed",
                                                files_removed = 0,
                                                stringsAsFactors = FALSE))
          })
        }
      }

      if (nrow(cleanup_summary) > 0) {
        total_removed <- sum(cleanup_summary$files_removed)
        successful_removals <- sum(cleanup_summary$status == "removed")
        cli_log(sprintf("Cleanup completed: %d directories processed, %d files removed",
                        successful_removals, total_removed), "SUCCESS")
      }

    } else {
      cli_log("Cleanup disabled - intermediate files preserved for inspection", "INFO", bullet = TRUE)
    }

    # Pipeline completion
    log_end_time <- Sys.time()
    runtime <- difftime(log_end_time, log_start_time, units = "mins")

    cli_log("Pipeline Statistics", "INFO")
    cli_log(sprintf("Total runtime: %.2f minutes", runtime), bullet = TRUE)

    if (!is.null(outputs$read_classes) && nrow(outputs$read_classes) > 0) {
      total_reads <- nrow(outputs$read_classes)
      cli_log(sprintf("Total reads processed: %d", total_reads), bullet = TRUE)

      if ("tail_type" %in% colnames(outputs$read_classes)) {
        polya_count <- sum(outputs$read_classes$tail_type == "polyA", na.rm = TRUE)
        polyt_count <- sum(outputs$read_classes$tail_type == "polyT", na.rm = TRUE)
        cli_log(sprintf("PolyA reads: %d, PolyT reads: %d", polya_count, polyt_count), bullet = TRUE)
      }
    }

    if (!is.null(outputs$nonadenosine_residues) && nrow(outputs$nonadenosine_residues) > 0) {
      total_modifications <- nrow(outputs$nonadenosine_residues)
      cli_log(sprintf("Non-adenosine modifications detected: %d", total_modifications), bullet = TRUE)
    }

    cli_log("Pipeline completed", "SUCCESS")
    cli_log(sprintf("Log file saved at: %s", log_filepath), "INFO", bullet = TRUE)

    # Final status messages (same as DRS pipeline)
    cli::cli_rule()
    if (warn_message) {
      cli::cli_alert_warning(paste(
        "Ninetails cDNA pipeline completed with warnings.",
        "Please check the log file for details."
      ))
    } else {
      cli_log("cDNA Pipeline completed successfully", "SUCCESS")
    }

    cli::cli_text()
    cli::cli_text(cli::col_grey("Thank you for using Ninetails."))

    return(outputs)

  }, error = function(e) {
    cli::cli_alert_danger("Error: {e$message}")
    cli_log(sprintf("Pipeline error: %s", e$message), "ERROR")
    cli::cli_alert_danger("Ninetails cDNA pipeline aborted")
    return(invisible(NULL))
  })

}


