################################################################################
# WRAPPER FUNCTION FOR DORADO OUTPUTS
################################################################################
#' Complete Oxford Nanopore poly(A) tail analysis pipeline for Dorado DRS data.
#'
#' This comprehensive wrapper function orchestrates the complete analysis of Oxford Nanopore
#' direct RNA sequencing (DRS) data processed with Dorado basecaller (>=1.0.0) using POD5
#' file format. The pipeline identifies and characterizes non-adenosine nucleotides within
#' poly(A) tails through advanced signal processing, machine learning-based classification,
#' and statistical analysis.
#'
#' @section Pipeline Overview:
#' The analysis pipeline consists of several integrated stages:
#' \enumerate{
#'   \item \strong{Input Validation}: Validates Dorado summary files and POD5 directories
#'   \item \strong{Data Preprocessing}: Splits large datasets and applies quality filters
#'   \item \strong{Signal Extraction}: Extracts poly(A) tail signals from POD5 files
#'   \item \strong{Feature Engineering}: Computes pseudomoves and signal characteristics
#'   \item \strong{Segmentation}: Identifies candidate modification regions
#'   \item \strong{GAF Creation}: Creates Gramian Angular Fields for CNN input
#'   \item \strong{Classification}: Applies trained neural networks for prediction
#'   \item \strong{Output Creation}: Produces comprehensive results and statistics
#'   \item \strong{Cleanup}: Optionally removes intermediate files to save disk space
#' }
#'
#' @section Input Requirements:
#' This pipeline requires specific input formats:
#' \itemize{
#'   \item \strong{Dorado Summary}: Must contain poly(A) information columns:
#'     \code{read_id}, \code{filename}, \code{poly_tail_length},
#'     \code{poly_tail_start}, \code{poly_tail_end}
#'   \item \strong{POD5 Files}: Raw signal files corresponding to reads in summary
#' }
#'
#' @section Output Structure:
#' The function creates several output subdirectories:
#' \describe{
#'   \item{\code{dorado_summary_dir/}}{Split summary files for parallel processing}
#'   \item{\code{polya_signal_dir/}}{Extracted poly(A) signals in RDS format}
#'   \item{\code{nonA_temp_dir/}}{Intermediate non-A prediction files}
#'   \item{\code{polya_chunks_dir/}}{Signal segmentation data}
#' }
#'
#' @section Performance Characteristics:
#' \itemize{
#'   \item \strong{Scalability}: Efficient parallel processing across multiple cores
#'   \item \strong{Memory Management}: Smart data partitioning prevents memory overflow
#'   \item \strong{Progress Tracking}: Comprehensive logging and progress indicators
#'   \item \strong{Error Handling}: Robust error recovery and informative diagnostics
#' }
#'
#' @section Quality Control:
#' When \code{qc = TRUE}, the pipeline applies several quality filters:
#' \itemize{
#'   \item Poly(A) tail length filtering (minimum 10 nucleotides)
#'   \item Coordinate validation (start < end positions)
#'   \item Terminal position masking to reduce false positives
#'   \item Duplicate read detection and handling
#' }
#'
#' @section Algorithm Details:
#' The pipeline employs several sophisticated algorithms:
#' \itemize{
#'   \item \strong{Pseudomove Computation}: Z-score based outlier detection with adaptive windowing
#'   \item \strong{Signal Segmentation}: Run-length encoding for modification region identification
#'   \item \strong{GAF Transformation}: Gramian Angular Field generation for CNN compatibility
#'   \item \strong{Neural Classification}: Pre-trained CNNs for nucleotide type prediction
#' }
#'
#' @param dorado_summary Character string or data frame. Either the full path to a
#'   Dorado summary file (.txt, .tsv, .csv) or an in-memory data frame containing
#'   summary information. \strong{Required columns}: \code{read_id} (unique identifier),
#'   \code{filename} (POD5 file name), \code{poly_tail_length} (tail length in nucleotides),
#'   \code{poly_tail_start} (start coordinate), \code{poly_tail_end} (end coordinate).
#'   Additional columns such as \code{alignment_genome}, \code{alignment_mapq} are
#'   automatically included if present.
#'
#' @param pod5_dir Character string. Full path to the directory containing POD5 files.
#'   The directory should contain POD5 files referenced in the \code{filename} column
#'   of the Dorado summary. Files can be organized in subdirectories (recursive search
#'   is performed). Typical ONT directory structures are automatically handled.
#'
#' @param num_cores Integer [1]. Number of physical CPU cores to use for parallel
#'   processing. Should not exceed \code{parallel::detectCores() - 1} to maintain
#'   system responsiveness. Higher core counts significantly reduce processing time
#'   for large datasets. Memory usage scales approximately linearly with core count.
#'
#' @param qc Logical [TRUE]. Enable comprehensive quality control filtering.
#'   When \code{TRUE}, applies stringent filters including:
#'   \itemize{
#'     \item Minimum poly(A) tail length (10 nucleotides)
#'     \item Coordinate validation and sanitization
#'     \item Terminal position masking (reduces false positives)
#'     \item Statistical outlier detection and removal
#'   }
#'   Set to \code{FALSE} only for preliminary analysis or when using pre-filtered data.
#'
#' @param save_dir Character string. Full path to the output directory where all
#'   results, intermediate files, and logs will be saved. The directory will be
#'   created if it doesn't exist. If the directory contains existing files, the
#'   user will be prompted to choose whether to abort or overwrite. Requires
#'   sufficient disk space (depending on dataset size).
#'
#' @param prefix Character string [""]. Optional prefix for output filenames.
#'   When provided, this string is prepended to all output files between the
#'   timestamp and file type identifier. Useful for organizing multiple analyses
#'   or experimental conditions. Should contain only filesystem-safe characters
#'   (alphanumeric, underscore, hyphen).
#'
#' @param part_size Integer [40000]. Maximum number of reads to process in each
#'   file partition. Must be >= 1. Larger values increase memory usage but may
#'   improve processing efficiency. Smaller values reduce memory footprint and
#'   enable processing of very large datasets on memory-constrained systems.
#'   Optimal values typically range from 10,000 to 100,000 depending on available RAM.
#'
#' @param cleanup Logical [FALSE]. Controls removal of intermediate files after
#'   successful analysis completion. When \code{TRUE}, removes all temporary
#'   subdirectories (\code{dorado_summary_dir}, \code{polya_signal_dir},
#'   \code{nonA_temp_dir}, \code{polya_chunks_dir}) keeping only final results
#'   and log files. When \code{FALSE}, preserves all intermediate files for
#'   debugging or detailed inspection. Recommended to keep \code{TRUE} to save disk space.
#'
#' @return A named list containing comprehensive analysis results:
#' \describe{
#'   \item{\code{read_classes}}{Data frame with per-read classification results including:
#'     \itemize{
#'       \item \code{readname}: Unique read identifier
#'       \item \code{contig}: Reference genome/transcript name (if aligned)
#'       \item \code{polya_length}: Estimated poly(A) tail length
#'       \item \code{qc_tag}: Quality control status
#'       \item \code{class}: Classification result ("decorated", "blank", "unclassified")
#'       \item \code{comments}: Detailed classification rationale using standard codes
#'     }
#'   }
#'   \item{\code{nonadenosine_residues}}{Data frame with detailed information about
#'     detected non-A nucleotides including:
#'     \itemize{
#'       \item \code{readname}: Read identifier
#'       \item \code{contig}: Reference information
#'       \item \code{prediction}: Predicted nucleotide type (C, G, U)
#'       \item \code{est_nonA_pos}: Estimated position within poly(A) tail
#'       \item \code{polya_length}: Total tail length
#'       \item \code{qc_tag}: Quality metrics
#'     }
#'   }
#' }
#'
#' @section Classification Codes:
#' The \code{comments} column in \code{read_classes} uses standardized codes:
#' \describe{
#'   \item{YAY}{Non-A residue detected (successful classification)}
#'   \item{MAU}{Move transition absent, no non-A residue detected}
#'   \item{MPU}{Move transition present, but no non-A residue detected}
#'   \item{IRL}{Insufficient read length for reliable analysis}
#'   \item{QCF}{Quality control filtering failed}
#' }
#'
#' @section System Requirements:
#' \itemize{
#'   \item \strong{R Version}: >= 3.5.0
#'   \item \strong{Python}: >= 3.6 with pod5 module (\code{pip install pod5})
#'   \item \strong{Memory}: >= 8GB RAM (16GB+ recommended for large datasets)
#'   \item \strong{Storage}: Temporary space ~2-5x input file size
#'   \item \strong{Dependencies}: Tensorflow/Keras for neural network inference
#' }
#'
#' @section Error Handling:
#' The function implements comprehensive error handling:
#' \itemize{
#'   \item Input validation with informative error messages
#'   \item Graceful handling of corrupted or missing files
#'   \item Automatic retry mechanisms for transient failures
#'   \item Detailed logging of all errors and warnings
#'   \item Safe cleanup of temporary files on failure
#' }
#'
#' @section Performance Tips:
#' \itemize{
#'   \item Use SSDs for \code{save_dir} to improve I/O performance
#'   \item Adjust \code{part_size} based on available RAM (larger = faster, more memory)
#'   \item Use \code{num_cores = parallel::detectCores() - 1} for maximum speed
#'   \item Ensure POD5 files and output directory are on fast storage
#'   \item Consider preprocessing very large datasets (>1M reads) in batches
#' }
#'
#' @family pipeline_functions
#' @family dorado_functions
#'
#' @seealso
#' \code{\link{check_tails_guppy}} for legacy Guppy-based analysis,
#' \code{\link{preprocess_inputs}} for input preprocessing,
#' \code{\link{process_dorado_signal_files}} for signal processing,
#' \code{\link{create_outputs_dorado}} for output generation,
#' \code{\link{filter_signal_by_threshold}} for signal analysis
#'
#' @export
#'
#' @note
#' \strong{Important considerations}:
#' \itemize{
#'   \item This function may take several hours for large datasets (>100K reads)
#'   \item Ensure sufficient disk space (typically 2-5x input size) in \code{save_dir}
#'   \item The function generates detailed log files for troubleshooting
#'   \item Results should always be assigned to a variable to prevent console overflow
#'   \item POD5 files must correspond exactly to reads in the Dorado summary
#'   \item For datasets >1M reads, consider batch processing or increased \code{part_size}
#' }
#'
#' @examples
#' \dontrun{
#' results <- check_tails_dorado_DRS(
#'   dorado_summary = "path/to/alignment_summary.txt",
#'   pod5_dir = "path/to/pod5/",
#'   num_cores = 2,
#'   qc = TRUE,
#'   save_dir = "path/to/output/",
#'   prefix = "experiment1",
#'   part_size = 40000,
#'   cleanup = TRUE
#' )
#' }
check_tails_dorado_DRS <- function(dorado_summary,
                                   pod5_dir,
                                   num_cores = 1,
                                   qc = TRUE,
                                   save_dir,
                                   prefix = "",
                                   part_size = 40000,
                                   cleanup = FALSE) {

  # Initialize warning flag
  warn_message <- FALSE

  # Initialize log collector (needed for directory checking)
  log_collector <- vector("character")
  log_start_time <- Sys.time()
  log_filename <- format(log_start_time,
                         paste0("%Y-%m-%d_%H-%M-%S",
                                if(nchar(prefix) > 0) paste0("_", prefix) else "",
                                "_ninetails.log"))

  # Create a temporary log message function for directory checking
  temp_log_message <- function(message, type = "INFO", section = NULL) {
    timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")

    if (!is.null(section)) {
      header <- paste0("\n-- ", section, " --\n")
      log_collector <<- c(log_collector, header)
    }

    log_entry <- sprintf("%s %s: %s", timestamp, type, message)
    log_collector <<- c(log_collector, log_entry)
    # For now, also print to console since file doesn't exist yet
    cat(sprintf("%s\n", log_entry))
  }

  # Check and handle output directory BEFORE creating it
  #####################################################
  should_proceed <- ninetails::check_output_directory(save_dir, temp_log_message)

  if (!should_proceed) {
    cat("Analysis aborted by user.\n")
    return(invisible(NULL))
  }

  # Now set up the proper log file path
  log_filepath <- file.path(save_dir, log_filename)

  # Helper functions for logging (same as original)
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

  #####################################################
  # PIPELINE START
  #####################################################
  tryCatch({
    # Initialize pipeline and log file with header
    ###################################################
    cli::cli_h1("Ninetails {utils::packageVersion('ninetails')}")
    cli::cli_alert_info(paste0("Welcome to Ninetails! ", {utils::packageVersion('ninetails')}))
    cli::cli_alert_info("Analysis toolkit for poly(A) tail composition in ONT data")

    # Write log file header with welcome message (only once)
    cat(sprintf("Ninetails Analysis Log - %s\n", format(log_start_time, "%Y-%m-%d %H:%M:%S")),
        file = log_filepath, append = FALSE)
    cat(sprintf("Ninetails v%s - Analysis toolkit for poly(A) tail composition in ONT data\n",
                utils::packageVersion('ninetails')),
        file = log_filepath, append = TRUE)
    cat("=====================================\n\n", file = log_filepath, append = TRUE)

    # Write any messages from directory checking to the log file
    if (length(log_collector) > 0) {
      for (entry in log_collector) {
        cat(paste0(entry, "\n"), file = log_filepath, append = TRUE)
      }
    }

    cli_log("Launching *check_tails_dorado_DRS* pipeline", "INFO","Launching *check_tails_dorado_DRS* pipeline")

    # Configuration logging display
    #####################################################
    cli_log("Configuration", "INFO", "Configuration")
    cli_log(sprintf("Dorado summary:              %s",
                    if(!is.object(dorado_summary)) dorado_summary else deparse(substitute(dorado_summary))), bullet = TRUE)
    cli_log(sprintf("Pod5 files directory:        %s", pod5_dir), bullet = TRUE)
    cli_log(sprintf("Number of cores:             %d", num_cores), bullet = TRUE)
    cli_log(sprintf("Output quality control:      %s", qc), bullet = TRUE)
    cli_log(sprintf("Output directory:            %s", save_dir), bullet = TRUE)
    cli_log(sprintf("Poly(A) processed at once:   %s", part_size), bullet = TRUE)
    cli_log(sprintf("Cleanup intermediate files:  %s", cleanup), bullet = TRUE)

    # Preprocess input files
    #####################################################
    processed_files <- ninetails::preprocess_inputs(
      dorado_summary = dorado_summary,
      pod5_dir = pod5_dir,
      num_cores = num_cores,
      qc = qc,
      save_dir = save_dir,
      prefix = prefix,
      part_size = part_size,
      cli_log = cli_log
    )

    # Process Dorado data
    #####################################################
    polya_signal_files <- processed_files$polya_signal_files
    nonA_temp_dir <- file.path(save_dir, "nonA_temp_dir")
    polya_chunks_dir <- file.path(save_dir, "polya_chunks_dir")
    nonA_results <- ninetails::process_dorado_signal_files(
      polya_signal_files,
      nonA_temp_dir,
      polya_chunks_dir,
      num_cores,
      cli_log)

    # Run create_outputs to generate tabular outputs
    #####################################################
    cli_log("Creating final outputs...", "INFO", "Creating Output", bullet = TRUE)
    dorado_summary_dir <- file.path(save_dir, "dorado_summary_dir")

    # Ensure create_outputs_dorado processes all files and returns proper format
    outputs <- tryCatch({
      result <- ninetails::create_outputs_dorado(
        dorado_summary_dir = dorado_summary_dir,
        nonA_temp_dir = nonA_temp_dir,
        polya_chunks_dir = polya_chunks_dir,
        num_cores = num_cores,
        qc = qc,
        original_summary = dorado_summary
      )

      # Validate output format - should be a list with read_classes and nonadenosine_residues
      if (!is.list(result)) {
        stop("create_outputs_dorado must return a list", call. = FALSE)
      }

      # Check for expected components
      expected_names <- c("read_classes", "nonadenosine_residues")
      if (!all(expected_names %in% names(result))) {
        cli_log("Warning: Output does not contain expected components. Found components:", "WARNING")
        cli_log(paste(names(result), collapse = ", "), "WARNING")
      }

      result
    }, error = function(e) {
      cli_log(sprintf("Error in output creation: %s", e$message), "ERROR")
      stop(e$message, call. = FALSE)
    })
    cli_log("Output creation completed", "SUCCESS")

    # Save final outputs
    #####################################################
    cli_log("Saving final outputs...", "INFO", "Saving Results", bullet = TRUE)

    # Use the same save_outputs function as in guppy pipeline
    ninetails::save_outputs(outputs, save_dir, prefix)

    # Cleanup intermediate files if requested
    #####################################################
    if (cleanup) {
      cli_log("Cleaning up intermediate files...", "INFO", "Cleanup", bullet = TRUE)

      # List of intermediate directories to remove
      intermediate_dirs <- c(
        file.path(save_dir, "dorado_summary_dir"),
        file.path(save_dir, "polya_signal_dir"),
        file.path(save_dir, "nonA_temp_dir"),
        file.path(save_dir, "polya_chunks_dir")
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
            # Count files before removal
            files_count <- length(list.files(dir_path, recursive = TRUE))

            # Remove directory and all contents
            unlink(dir_path, recursive = TRUE)

            # Verify removal
            if (!dir.exists(dir_path)) {
              cli_log(sprintf("Removed: %s (%d files)", basename(dir_path), files_count), "SUCCESS")
              cleanup_summary <- rbind(cleanup_summary,
                                       data.frame(directory = basename(dir_path),
                                                  status = "removed",
                                                  files_removed = files_count,
                                                  stringsAsFactors = FALSE))
            } else {
              cli_log(sprintf("Failed to remove: %s", basename(dir_path)), "WARNING")
              cleanup_summary <- rbind(cleanup_summary,
                                       data.frame(directory = basename(dir_path),
                                                  status = "failed",
                                                  files_removed = 0,
                                                  stringsAsFactors = FALSE))
            }
          }, error = function(e) {
            cli_log(sprintf("Error removing %s: %s", basename(dir_path), e$message), "WARNING")
            cleanup_summary <<- rbind(cleanup_summary,
                                      data.frame(directory = basename(dir_path),
                                                 status = "error",
                                                 files_removed = 0,
                                                 stringsAsFactors = FALSE))
          })
        } else {
          cli_log(sprintf("Directory not found: %s", basename(dir_path)), "INFO")
        }
      }

      # Summary of cleanup
      total_removed <- sum(cleanup_summary$files_removed)
      successful_removals <- sum(cleanup_summary$status == "removed")

      if (successful_removals > 0) {
        cli_log(sprintf("Cleanup completed: %d directories removed, %d files deleted",
                        successful_removals, total_removed), "SUCCESS")
      } else {
        cli_log("No intermediate directories were removed", "INFO")
      }

    } else {
      cli_log("Cleanup disabled - all intermediate files preserved", "INFO", "Cleanup", bullet = TRUE)
      cli_log("Intermediate directories available for inspection:", "INFO")

      # List preserved directories
      intermediate_dirs <- c("dorado_summary_dir", "polya_signal_dir", "nonA_temp_dir", "polya_chunks_dir")
      for (dir_name in intermediate_dirs) {
        dir_path <- file.path(save_dir, dir_name)
        if (dir.exists(dir_path)) {
          file_count <- length(list.files(dir_path, recursive = TRUE))
          cli_log(sprintf("- %s (%d files)", dir_name, file_count), "INFO")
        }
      }
    }

    # Pipeline statistics
    #####################################################
    log_end_time <- Sys.time()
    runtime <- difftime(log_end_time, log_start_time, units = "mins")

    cli_log("Pipeline Statistics", "INFO")
    cli_log(sprintf("Total runtime: %.2f minutes", runtime), bullet = TRUE)
    cli_log("Pipeline completed", "SUCCESS")
    cli_log(sprintf("Log file saved at: %s", log_filepath), "INFO", bullet = TRUE)

    # Final status messages
    cli::cli_rule()
    cli_log("Pipeline completed successfully", "SUCCESS")
    cli_log(sprintf("Output files saved in: %s", save_dir), "INFO", bullet = TRUE)

    if (warn_message) {
      cli::cli_alert_warning(paste(
        "Ninetails exited with WARNING.\n\n",
        "Check your Dorado outputs and POD5 files for consistency. Some reads may have been omitted from analysis."
      ))
    }

    cli::cli_text()
    cli::cli_text(cli::col_grey("Thank you for using Ninetails."))

    return(outputs)

  }, error = function(e) {
    cli::cli_alert_danger("Error: {e$message}")
    cli_log(sprintf("Pipeline error: %s", e$message), "ERROR")
    cli::cli_alert_danger("Ninetails aborted")
    return(invisible(NULL))
  })
}


#' Process Dorado poly(A) signal files for non-A prediction and tail chunk extraction
#'
#' This function takes a set of Dorado poly(A) signal files (RDS format),
#' extracts features, segments tails, generates Gramian Angular Fields (GAFs),
#' applies classification models to predict non-A residues, and saves both
#' prediction results and tail chunk lists as RDS files in specified output directories.
#' Progress is logged using a cli logging function. This function is not intended
#' to be used outside the pipeline wrapper.
#'
#' @param polya_signal_files Character vector of paths to poly(A) signal files
#' (RDS format) generated from Dorado outputs.
#' @param nonA_temp_dir Character path to the directory where non-A prediction
#' results will be written. Directory will be created if it does not exist.
#' @param polya_chunks_dir Character path to the directory where tail chunk list
#' RDS files will be written. Directory will be created if it does not exist.
#' @param num_cores Integer number of CPU cores to use for parallel feature
#' extraction and downstream computations.
#' @param cli_log Function used for logging messages.
#'
#' @returns An (invisible) list summarizing processing results for each input file.
#' Each element contains:
#'   \itemize{
#'     \item \code{success} (logical): whether the file was processed successfully
#'     \item \code{file} (character): path to the non-A prediction RDS file if successful
#'     \item \code{chunks_file} (character): path to the tail chunk list RDS file if successful
#'     \item \code{error} (character): error message if processing failed
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage
#' signal_files <- list.files("signals/", pattern = "polya_signal.*\\.rds$", full.names = TRUE)
#' results <- process_dorado_signal_files(
#'   polya_signal_files = signal_files,
#'   nonA_temp_dir = "nonA_predictions/",
#'   polya_chunks_dir = "tail_chunks/",
#'   num_cores = 4,
#'   cli_log = function(msg, level, ...) message(sprintf("[%s] %s", level, msg))
#' )
#' }
process_dorado_signal_files <- function(polya_signal_files,
                                        nonA_temp_dir,
                                        polya_chunks_dir,
                                        num_cores,
                                        cli_log) {

  # Assertions
  if (!dir.exists(nonA_temp_dir)) {
    dir.create(nonA_temp_dir, recursive = TRUE)
  }
  if (!dir.exists(polya_chunks_dir)) {
    dir.create(polya_chunks_dir, recursive = TRUE)
  }
  results_summary <- list()
  for (signal_file in polya_signal_files) {
    input_basename <- basename(signal_file)
    output_basename <- sub("polya_signal", "nonA_pred", input_basename)
    output_file <- file.path(nonA_temp_dir, output_basename)
    chunks_basename <- sub("polya_signal", "tail_chunks", input_basename)
    chunks_file <- file.path(polya_chunks_dir, chunks_basename)

    # Section header for each file
    #####################################################
    cli_log(sprintf("-- SIGNAL DATA ANALYSIS: %s --", input_basename), "INFO", sprintf("SIGNAL DATA ANALYSIS: %s", input_basename))
    cli_log(sprintf("Processing %s for nonA prediction", input_basename), "INFO")

    result <- tryCatch({
      # Computing Pseudomoves
      #####################################################
      cli_log("Computing pseudomoves...", "INFO", "Computing Pseudomoves", bullet = TRUE)
      signal_list <- readRDS(signal_file)
      features <- ninetails::create_tail_features_list_dorado(signal_list, num_cores = num_cores)
      if (!is.null(features$zeromoved_readnames) && length(features$zeromoved_readnames) > 0) {
        cli_log(sprintf("%d reads had zero moves", length(features$zeromoved_readnames)), "WARNING")
      }
      if (!is.null(features$nonpseudomoved_readnames) && length(features$nonpseudomoved_readnames) > 0) {
        cli_log(sprintf("%d reads did not meet pseudomove criteria", length(features$nonpseudomoved_readnames)), "WARNING")
      }
      cli_log("Feature extraction completed", "SUCCESS")

      # Creating Chunks
      #####################################################
      cli_log("Creating tail segmentation data...", "INFO", "Creating Chunks", bullet = TRUE)
      chunks <- ninetails::create_tail_chunk_list_dorado(features, num_cores = num_cores)
      cli_log("Chunk creation completed", "SUCCESS")
      saveRDS(chunks, chunks_file)
      cli_log(sprintf("Saved tail chunk list to %s", chunks_file), "SUCCESS")

      # Generating GAFs
      #####################################################
      cli_log("Computing gramian angular fields...", "INFO", "Generating GAFs", bullet = TRUE)
      gafs <- ninetails::create_gaf_list(chunks, num_cores = num_cores)
      cli_log("GAF computation completed", "SUCCESS")

      # Predicting Classes
      #####################################################
      cli_log("Running predictions...", "INFO", "Predicting Classes", bullet = TRUE)
      preds <- ninetails::predict_gaf_classes(gafs)
      cli_log("Predictions completed", "SUCCESS")

      saveRDS(preds, output_file)
      cli_log(sprintf("Saved predictions to %s", output_file), "SUCCESS")
      list(success = TRUE, file = output_file, chunks_file = chunks_file)
    }, error = function(e) {
      cli_log(sprintf("Error processing %s: %s", input_basename, e$message), "ERROR")
      list(success = FALSE, error = e$message)
    })
    results_summary[[input_basename]] <- result
  }
  cli_log("nonA prediction completed for all signal files", "SUCCESS")
  invisible(results_summary)
}
