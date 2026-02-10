################################################################################
# WRAPPER FUNCTION FOR GUPPY/NANOPOLISH(TAILFINDR) OUTPUTS
################################################################################
#' Wrapper function for complete DRS processing by ninetails package
#' (legacy mode).
#'
#' Orchestrates the full ninetails analysis pipeline for direct RNA
#' sequencing (DRS) data basecalled with Guppy and characterised by
#' nanopolish (or tailfindr). The function validates inputs, extracts
#' poly(A) tail features from multi-Fast5 files, segments tails into
#' chunks, computes Gramian Angular Fields (GAFs), classifies reads
#' with a pretrained CNN, and writes output files.
#'
#' @details
#' \strong{Legacy mode.} Due to Oxford Nanopore Technologies' transition
#' from Guppy to Dorado, from Fast5 to POD5, and from R9 to R10
#' chemistry, this pipeline is maintained for backward compatibility
#' only and will not be further optimised. For current data, use
#' \code{\link{check_tails_dorado_DRS}} instead.
#'
#' The function accepts either nanopolish or tailfindr output (DRS only).
#' Tailfindr input is automatically converted to nanopolish-like format
#' via \code{\link{check_polya_length_filetype}}.
#'
#' For large input files (exceeding \code{part_size} rows), the poly(A)
#' table is split into parts and processed sequentially via
#' \code{\link{process_polya_parts}} to avoid memory overflow.
#'
#' The pipeline steps are:
#' \enumerate{
#'   \item Input validation and Fast5 format checking.
#'   \item Feature extraction
#'     (\code{\link{create_tail_feature_list}}).
#'   \item Tail segmentation
#'     (\code{\link{create_tail_chunk_list}}).
#'   \item GAF computation (\code{\link{create_gaf_list}}).
#'   \item CNN classification
#'     (\code{\link{predict_gaf_classes}}).
#'   \item Output creation (\code{\link{create_outputs}}).
#' }
#'
#' The \code{comments} column in the output encodes detailed read
#' classification:
#' \describe{
#'   \item{IRL}{Insufficient read length.}
#'   \item{QCF}{Nanopolish QC failed.}
#'   \item{MAU}{Move transition absent, non-A residue undetected.}
#'   \item{MPU}{Move transition present, non-A residue undetected.}
#'   \item{NIN}{Not included in analysis (\code{pass_only = TRUE}).}
#'   \item{YAY}{Move transition present, non-A residue detected.}
#' }
#'
#' Filtering criteria: basecaller move value = 1,
#' \code{qc_tag \%in\% c("PASS")} (or \code{c("PASS", "SUFFCLIP")}),
#' and nanopolish-estimated tail length >= 10 nt.
#'
#' @param polya_data Character string or data frame. Full path of the
#'   \code{.tsv} file produced by nanopolish polya or tailfindr (DRS
#'   only), or an in-memory data frame. Tailfindr input is converted
#'   automatically.
#'
#' @param sequencing_summary Character string or data frame. Full path
#'   of the \code{.txt} sequencing summary file, or an in-memory data
#'   frame.
#'
#' @param workspace Character string. Full path of the directory
#'   containing basecalled multi-Fast5 files.
#'
#' @param num_cores Numeric \code{[1]}. Number of physical cores to use.
#'   Do not exceed 1 less than the number of cores at your disposal.
#'
#' @param basecall_group Character string \code{["Basecall_1D_000"]}.
#'   Name of the level in the Fast5 file hierarchy from which data
#'   should be extracted.
#'
#' @param pass_only Logical \code{[TRUE]}. If \code{TRUE}, only reads
#'   tagged by nanopolish as \code{"PASS"} are included. If
#'   \code{FALSE}, reads tagged \code{"PASS"} or \code{"SUFFCLIP"} are
#'   included.
#'
#' @param qc Logical \code{[TRUE]}. If \code{TRUE}, terminal non-A
#'   residue positions likely arising from nanopolish segmentation
#'   errors are labelled with a \code{"-WARN"} suffix. It is up to the
#'   user whether to include or discard such reads.
#'
#' @param save_dir Character string. Full path of the directory where
#'   output files should be stored.
#'
#' @param prefix Character string (optional, default \code{""}). If
#'   provided, inserted into output file names between the timestamp and
#'   the file-type suffix.
#'
#' @param part_size Numeric \code{[1000000]}. Maximum number of rows
#'   processed at once. Must be >= 1000. If the input exceeds this
#'   value, it is split and processed sequentially.
#'
#' @return A named list with two data frames:
#' \describe{
#'   \item{read_classes}{Data frame. Per-read classification including
#'     columns \code{class} and \code{comments}.}
#'   \item{nonadenosine_residues}{Data frame. Detailed positional
#'     information for all detected non-A residues.}
#' }
#' A log file and TSV output files (\code{read_classes},
#' \code{nonadenosine_residues}) are also written to \code{save_dir}.
#'
#' @seealso \code{\link{check_tails_dorado_DRS}} for the current
#'   (Dorado-based) pipeline,
#'   \code{\link{process_polya_complete}} and
#'   \code{\link{process_polya_parts}} for the internal processing
#'   functions,
#'   \code{\link{check_polya_length_filetype}} for input format
#'   detection,
#'   \code{\link{create_outputs}} for the output assembly step.
#'
#' @importFrom foreach %dopar%
#' @importFrom utils head
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' results <- ninetails::check_tails_guppy(
#'   polya_data = system.file('extdata', 'test_data',
#'                            'nanopolish_output.tsv',
#'                            package = 'ninetails'),
#'   sequencing_summary = system.file('extdata', 'test_data',
#'                                    'sequencing_summary.txt',
#'                                    package = 'ninetails'),
#'   workspace = system.file('extdata', 'test_data',
#'                           'basecalled_fast5',
#'                           package = 'ninetails'),
#'   num_cores = 2,
#'   basecall_group = 'Basecall_1D_000',
#'   pass_only = TRUE,
#'   qc = TRUE,
#'   save_dir = '~/Downloads',
#'   prefix = "prefix",
#'   part_size = 2000)
#'
#' }
check_tails_guppy <- function(polya_data,
                              sequencing_summary,
                              workspace,
                              num_cores=1,
                              basecall_group="Basecall_1D_000",
                              pass_only=TRUE,
                              qc=TRUE,
                              save_dir,
                              prefix="",
                              part_size=1000000) {

  # Initialize warning flag
  warn_message <- FALSE

  # Create output directory if needed
  tryCatch({
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }
  }, error = function(e) {
    cli::cli_alert_danger("Failed to create output directory: {e$message}")
    stop(e$message, call. = FALSE)
  })

  # Initialize log collector
  log_collector <- vector("character")
  log_start_time <- Sys.time()
  log_filename <- format(log_start_time,
                         paste0("%Y-%m-%d_%H-%M-%S",
                                if(nchar(prefix) > 0) paste0("_", prefix) else "",
                                "_ninetails.log"))
  log_filepath <- file.path(save_dir, log_filename)

  # Helper functions for logging
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
    # Initialize pipeline
    ###################################################
    cli::cli_h1("Ninetails {utils::packageVersion('ninetails')}")
    cli_log(paste0("Welcome to Ninetails! ", {utils::packageVersion('ninetails')}), bullet = TRUE)
    cli_log("Analysis toolkit for poly(A) tail composition in ONT data", bullet = TRUE)

    cli_log("Launching *check_tails_guppy* pipeline", "INFO","Launching *check_tails_guppy* pipeline")

    # Configuration logging display
    #####################################################
    cli_log("Configuration", "INFO", "Configuration")
    cli_log(sprintf("Poly(A) length file:         %s",
                    if(!is.object(polya_data)) polya_data else deparse(substitute(polya_data))), bullet = TRUE)
    cli_log(sprintf("Sequencing summary:          %s",
                    if(!is.object(sequencing_summary)) sequencing_summary else deparse(substitute(sequencing_summary))), bullet = TRUE)
    cli_log(sprintf("Fast5 files directory:       %s", workspace), bullet = TRUE)
    cli_log(sprintf("Number of cores:             %d", num_cores), bullet = TRUE)
    cli_log(sprintf("Basecall group:              %s", basecall_group), bullet = TRUE)
    cli_log(sprintf("Only PASS reads:             %s", pass_only), bullet = TRUE)
    cli_log(sprintf("Output quality control:      %s", qc), bullet = TRUE)
    cli_log(sprintf("Output directory:            %s", save_dir), bullet = TRUE)
    cli_log(sprintf("Poly(A) processed at once:   %s", part_size), bullet = TRUE)


    # Input validation and configuration
    cli_log("Validating input parameters", "INFO", "Validating Inputs")
    # Assertions
    ###################################################
    assert_condition(is.numeric(num_cores) && num_cores > 0,
                     "Number of cores must be a positive numeric value")

    assert_condition(is.character(workspace),
                     "Fast5 files directory path must be a character string")
    assert_dir_exists(workspace, "Fast5 files")

    assert_condition(is_string(basecall_group, min.chars = 1),
                     "Basecall group must be a non-empty character string")

    assert_condition(is.character(save_dir),
                     "Output directory path must be a character string")

    assert_condition(is.logical(pass_only),
                     "pass_only must be logical [TRUE/FALSE]")

    assert_condition(is.logical(qc),
                     "qc must be logical [TRUE/FALSE]")

    # Validate prefix (if provided)
    if (nchar(prefix) > 0) {
      assert_condition(is.character(prefix),
                       "File name prefix must be a character string")
    }

    # Validate nanopolish input
    if (is_string(polya_data)) {
      assert_file_exists(polya_data, "Poly(A) data")
    } else {
      assert_condition(is.data.frame(polya_data) && nrow(polya_data) > 0,
                       "Poly(A) data input must be a non-empty data frame or valid file path")
    }

    # Validate sequencing summary input
    if (is_string(sequencing_summary)) {
      assert_file_exists(sequencing_summary, "Sequencing summary")
    } else {
      assert_condition(is.data.frame(sequencing_summary) && nrow(sequencing_summary) > 0,
                       "Sequencing summary must be a non-empty data frame or valid file path")
    }
    cli_log("Provided input files/paths are in correct format", "SUCCESS")


    # Validate FAST5 format whether multifast5 & DRS & basecalled
    cli_log("Checking FAST5 format...", "INFO", "Validating FAST5 Data", bullet = TRUE)
    fast5_output <- utils::capture.output(
      ninetails::check_fast5_filetype(workspace, basecall_group)
    )
    log_message(paste(fast5_output, collapse = "\n"), "INFO")
    cli_log("FAST5 format validated; all data blocks present", "SUCCESS")


    # Check poly(A) file input size
    cli_log("Checking poly(A) data input size...", "INFO", "Poly(A) Input Size Check", bullet = TRUE)
    data_size <- if(is.character(polya_data)) {
      nrow(vroom::vroom(polya_data, show_col_types = FALSE))
    } else {
      nrow(polya_data)
    }

    cli_log("Input validation passed", "SUCCESS")


    # Process based on size using part_size parameter
    if (data_size > part_size) {
      cli_log(sprintf("Large input detected (%d rows). Processing in parts of %d rows...",
                      data_size, part_size), "INFO")

      # Split input
      part_files <- ninetails::split_polya_data(polya_data, part_size, save_dir)
      cli_log(sprintf("Input split into %d parts", length(part_files)), "INFO")

      # Process parts
      outputs <- ninetails::process_polya_parts(
        part_files = part_files,
        sequencing_summary = sequencing_summary,
        workspace = workspace,
        num_cores = num_cores,
        basecall_group = basecall_group,
        pass_only = pass_only,
        qc = qc,
        save_dir = save_dir,
        prefix = prefix,
        cli_log = cli_log
      )
    } else {
      cli_log(sprintf("Input size (%d rows) within limits. Processing as single file.",
                      data_size), "INFO")

      # Process complete dataset
      outputs <- ninetails::process_polya_complete(
        polya_data = polya_data,
        sequencing_summary = sequencing_summary,
        workspace = workspace,
        num_cores = num_cores,
        basecall_group = basecall_group,
        pass_only = pass_only,
        qc = qc,
        save_dir = save_dir,
        prefix = prefix,
        cli_log = cli_log
      )
    }

    # Save final outputs
    cli_log("Saving final outputs...", "INFO", "Saving Results", bullet = TRUE)
    ninetails::save_outputs(outputs, save_dir, prefix)

    # Pipeline statistics
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
        "Check your Nanopolish and Guppy outputs for consistency. It seems like some reads passing quality criteria, ",
        "which were provided in nanopolish output file (*.tsv) are absent from provided fast5 files. ",
        "They were omitted in this analysis. However they might constitute a significant fraction of your data. ",
        "If you are sure that everything was provided correctly, just ignore that info."
      ))
    }

    cli::cli_text()
    cli::cli_text(cli::col_grey("Thank you for using Ninetails."))

    return(outputs)

  }, error = function(e) {
    cli::cli_alert_danger("Error: {e$message}")
    cli::cli_alert_danger("Ninetails aborted")
  })
}



################################################################################
# INTERNAL FUNCTIONS REQUIRED FOR PIPELINE
################################################################################

#' Split large poly(A) data file into smaller parts.
#'
#' Divides a poly(A) length table (nanopolish or tailfindr output) into
#' smaller files of at most \code{part_size} rows each, saving them to a
#' \code{polya_data_parts} subdirectory within \code{save_dir}. This
#' prevents memory overflow when processing large datasets.
#'
#' @details
#' This function is called internally by
#' \code{\link{check_tails_guppy}} when the input exceeds the
#' \code{part_size} threshold. It is not intended to be called directly
#' by the user.
#'
#' @param polya_data Character string or data frame. Full path of the
#'   poly(A) length file (\code{.tsv}), or an in-memory data frame.
#'
#' @param part_size Numeric \code{[100000]}. Maximum number of rows per
#'   output part file.
#'
#' @param save_dir Character string. Full path of the directory where
#'   the \code{polya_data_parts} subdirectory will be created.
#'
#' @return Character vector of file paths to the created part files,
#'   named \code{polya_data_part_<i>_of_<n>.tsv}.
#'
#' @seealso \code{\link{check_tails_guppy}} which calls this function,
#'   \code{\link{process_polya_parts}} which processes the resulting
#'   parts.
#'
#' @keywords internal
#'
#' @export
split_polya_data <- function(polya_data, part_size = 100000, save_dir) {
  # Create subfolder for polya_parts
  parts_dir <- file.path(save_dir, "polya_data_parts")
  if (!dir.exists(parts_dir)) {
    dir.create(parts_dir, recursive = TRUE)
  }

  # Read data if file path provided
  if (is.character(polya_data)) {
    data <- vroom::vroom(polya_data, show_col_types = FALSE)
  } else {
    data <- polya_data
  }

  # Calculate number of polya_parts needed
  total_rows <- nrow(data)
  num_parts <- ceiling(total_rows / part_size)
  splitted_files <- character(num_parts)

  # Split and save polya_parts
  for (i in 1:num_parts) {
    start_idx <- ((i-1) * part_size) + 1
    end_idx <- min(i * part_size, total_rows)
    polya_part <- data[start_idx:end_idx, ]

    part_file <- file.path(parts_dir,
                           sprintf("polya_data_part_%d_of_%d.tsv", i, num_parts))
    utils::write.table(polya_part,
                       file = part_file,
                       sep = "\t",
                       row.names = FALSE,
                       quote = FALSE)
    splitted_files[i] <- part_file
  }

  return(splitted_files)
}



#' Save pipeline outputs to files.
#'
#' Writes the ninetails pipeline output data frames (\code{read_classes}
#' and \code{nonadenosine_residues}) to tab-separated files in
#' \code{save_dir}. File names include a timestamp and an optional
#' user-defined prefix.
#'
#' @details
#' This function is called internally by
#' \code{\link{check_tails_guppy}} at the end of the pipeline. It is
#' not intended to be called directly by the user.
#'
#' @param outputs Named list containing data frames to save. Expected
#'   names are \code{read_classes} and \code{nonadenosine_residues}.
#'
#' @param save_dir Character string. Full path of the directory where
#'   output files should be stored.
#'
#' @param prefix Character string (optional, default \code{""}). If
#'   provided, inserted into file names between the timestamp and the
#'   file-type suffix.
#'
#' @return Invisible \code{NULL}. Called for the side effect of writing
#'   files.
#'
#' @seealso \code{\link{check_tails_guppy}} which calls this function.
#'
#' @keywords internal
#'
#' @export
save_outputs <- function(outputs, save_dir, prefix = "") {
  tryCatch({
    mapply(function(x, y) {
      # Construct filename with optional prefix
      filename <- if (nchar(prefix) > 0) {
        paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), '_', prefix, '_', y, '.tsv')
      } else {
        paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), '_', y, '.tsv')
      }
      filepath <- file.path(save_dir, filename)
      utils::write.table(x,
                         file = filepath,
                         row.names = FALSE,
                         sep = "\t",
                         quote = FALSE)
      cli::cli_alert_success(sprintf("Saved %s", filename))
    }, outputs, names(outputs))
    invisible(NULL)
  }, error = function(e) {
    cli::cli_alert_danger(sprintf("Error saving outputs: %s", e$message))
    stop(e$message, call. = FALSE)
  })
}


#' Process a single (unsplit) poly(A) data file through the Guppy
#' pipeline.
#'
#' Runs the full ninetails analysis on a small-to-moderate poly(A) data
#' set (up to \code{part_size} rows). Performs format validation,
#' feature extraction, tail segmentation, GAF computation, CNN
#' classification, and output assembly.
#'
#' @details
#' This function is called internally by
#' \code{\link{check_tails_guppy}} (for inputs within the
#' \code{part_size} limit) and by
#' \code{\link{process_polya_parts}} (for each part of a split input).
#' It requires the \code{cli_log} closure defined inside
#' \code{\link{check_tails_guppy}} for formatted logging, and therefore
#' should not be called directly by the user.
#'
#' @param polya_data Character string or data frame. Full path of the
#'   poly(A) length file (\code{.tsv}), or an in-memory data frame.
#'
#' @param sequencing_summary Character string or data frame. Full path
#'   of the sequencing summary file, or an in-memory data frame.
#'
#' @param workspace Character string. Full path of the directory
#'   containing basecalled multi-Fast5 files.
#'
#' @param num_cores Numeric \code{[1]}. Number of physical cores.
#'
#' @param basecall_group Character string \code{["Basecall_1D_000"]}.
#'   Fast5 hierarchy level for data extraction.
#'
#' @param pass_only Logical \code{[TRUE]}. If \code{TRUE}, only
#'   \code{"PASS"} reads are included.
#'
#' @param qc Logical \code{[TRUE]}. If \code{TRUE}, terminal artefact
#'   positions are labelled with \code{"-WARN"}.
#'
#' @param save_dir Character string. Output directory path.
#'
#' @param prefix Character string (optional). Output file name prefix.
#'
#' @param cli_log Function. Logging closure defined in
#'   \code{\link{check_tails_guppy}} for formatted console and log file
#'   output.
#'
#' @param ... Additional arguments (currently unused).
#'
#' @return A named list with two data frames:
#' \describe{
#'   \item{read_classes}{Per-read classification data.}
#'   \item{nonadenosine_residues}{Detailed non-A residue positional
#'     data.}
#' }
#'
#' @seealso \code{\link{check_tails_guppy}} which calls this function,
#'   \code{\link{process_polya_parts}} for the split-input variant,
#'   \code{\link{check_polya_length_filetype}} for format detection,
#'   \code{\link{create_tail_feature_list}},
#'   \code{\link{create_tail_chunk_list}},
#'   \code{\link{create_gaf_list}},
#'   \code{\link{predict_gaf_classes}},
#'   \code{\link{create_outputs}} for the individual pipeline steps.
#'
#' @keywords internal
#'
#' @export
process_polya_complete <- function(polya_data,
                                   sequencing_summary,
                                   workspace,
                                   num_cores,
                                   basecall_group,
                                   pass_only,
                                   qc,
                                   save_dir,
                                   prefix,
                                   cli_log,
                                   ...) {

  # Initialize warning flag
  warn_message <- FALSE

  # Validate and convert poly(A) data format
  cli_log("Checking poly(A) length file format...", "INFO", "Processing Poly(A) length data", bullet = TRUE)

  polya_data <- tryCatch({
    result <- check_polya_length_filetype(polya_data)
    if (result$file_type == "nanopolish") {
      cli_log("Nanopolish input format detected", "INFO")
    } else if (result$file_type == "tailfindr_drs") {
      cli_log("Tailfindr DRS output detected. Converting to nanopolish-compatible format.", "WARNING")
      cli_log("NOTE: Tailfindr does not provide quality metrics for poly(A) tail calling.", "WARNING")
      cli_log("Results should be interpreted with caution.", "WARNING")
      cli_log("Quality tags are assigned based on tail length only (PASS for >= 10nt, NOREGION for < 10nt).", "WARNING")
    }
    cli_log("Input poly(A) format validated", "SUCCESS")
    result$data
  }, error = function(e) {
    cli_log(sprintf("Error in poly(A) input validation: %s", e$message), "ERROR")
    stop(e$message, call. = FALSE)
  })

  #####################################################
  # CREATING TAIL FEATURE LIST
  #####################################################
  cli_log("Creating tail feature list...", "INFO", "Creating Features", bullet = TRUE)
  tail_feature_list <- tryCatch({
    # Suppress text output while keeping progress bar visible
    invisible(utils::capture.output(
      features <- ninetails::create_tail_feature_list(
        nanopolish = polya_data,
        sequencing_summary = sequencing_summary,
        workspace = workspace,
        num_cores = num_cores,
        basecall_group = basecall_group,
        pass_only = pass_only
      )
    ))

    # Sanity checks
    #####################################################
    # Basic sanity check
    if (length(features$tail_feature_list) == 0) {
      stop("Produced feature list is of length 0. Check your input (nanopolish, fast5 files) for integrity.",
           call. = FALSE)
    }

    # Check for missing reads
    squiggle_names <- names(features$tail_feature_list)
    if (!is.null(squiggle_names) &&
        length(squiggle_names) != length(unique(squiggle_names))) {
      warn_message <<- TRUE  # Using <<- to modify the parent scope variable
    }

    # Report filtering results
    if (length(features$zeromoved_readnames) > 0) {
      cli_log(sprintf("%d reads had zero moves",
                      length(features$zeromoved_readnames)), "WARNING")
    }
    if (length(features$nonpseudomoved_readnames) > 0) {
      cli_log(sprintf("%d reads did not meet pseudomove criteria",
                      length(features$nonpseudomoved_readnames)), "WARNING")
    }

    features
  }, error = function(e) {
    cli_log(sprintf("Error in feature extraction: %s", e$message), "ERROR")
    stop(e$message, call. = FALSE)
  })
  cli_log("Feature extraction completed", "SUCCESS")

  #####################################################
  # CREATING TAIL CHUNK LIST
  #####################################################
  cli_log("Creating tail segmentation data...", "INFO", "Creating Chunks", bullet = TRUE)
  tail_chunk_list <- tryCatch({
    invisible(utils::capture.output(
      result <- ninetails::create_tail_chunk_list(
        tail_feature_list = tail_feature_list,
        num_cores = num_cores
      )
    ))
    result
  }, error = function(e) {
    cli_log(sprintf("Error in chunk creation: %s", e$message), "ERROR")
    stop(e$message, call. = FALSE)
  })
  cli_log("Chunk creation completed", "SUCCESS")

  #####################################################
  # CREATING GAF LIST
  #####################################################
  cli_log("Computing gramian angular fields...", "INFO", "Generating GAFs", bullet = TRUE)
  gaf_list <- tryCatch({
    invisible(utils::capture.output(
      result <- ninetails::create_gaf_list(
        tail_chunk_list = tail_chunk_list,
        num_cores = num_cores
      )
    ))
    result
  }, error = function(e) {
    cli_log(sprintf("Error in GAF computation: %s", e$message), "ERROR")
    stop(e$message, call. = FALSE)
  })
  cli_log("GAF computation completed", "SUCCESS")

  #####################################################
  # PREDICT CLASSES
  #####################################################
  cli_log("Running predictions...", "INFO", "Predicting Classes", bullet = TRUE)
  predicted_list <- tryCatch({
    invisible(utils::capture.output(
      result <- ninetails::predict_gaf_classes(gaf_list),
      type = "message"  # Captures message output from TensorFlow
    ))
    result
  }, error = function(e) {
    cli_log(sprintf("Error in class prediction: %s", e$message), "ERROR")
    stop(e$message, call. = FALSE)
  })
  cli_log("Predictions completed", "SUCCESS")

  #####################################################
  # CREATE OUTPUT
  #####################################################
  cli_log("Creating final outputs...", "INFO", "Creating Output", bullet = TRUE)
  outputs <- tryCatch({
    invisible(utils::capture.output(
      result <- ninetails::create_outputs(
        tail_feature_list = tail_feature_list,
        tail_chunk_list = tail_chunk_list,
        nanopolish = polya_data,
        predicted_list = predicted_list,
        num_cores = num_cores,
        pass_only = pass_only,
        qc = qc
      )
    ))
    result
  }, error = function(e) {
    cli_log(sprintf("Error in output creation: %s", e$message), "ERROR")
    stop(e$message, call. = FALSE)
  })
  cli_log("Output creation completed", "SUCCESS")

  return(outputs)
}

#' Process poly(A) data split into multiple parts through the Guppy
#' pipeline.
#'
#' Iterates over a set of poly(A) data part files (produced by
#' \code{\link{split_polya_data}}), processing each sequentially via
#' \code{\link{process_polya_complete}}, saving per-part outputs, and
#' merging all results into a single output list.
#'
#' @details
#' This function is called internally by
#' \code{\link{check_tails_guppy}} when the input exceeds the
#' \code{part_size} threshold. It requires the \code{cli_log} closure
#' defined inside \code{\link{check_tails_guppy}} for formatted
#' logging, and therefore should not be called directly by the user.
#'
#' For each part, intermediate results are saved as TSV files in a
#' dedicated subdirectory (\code{part_<i>_of_<n>}) within
#' \code{save_dir}. After all parts are processed, \code{read_classes}
#' and \code{nonadenosine_residues} tables are merged with
#' \code{rbind}.
#'
#' @param part_files Character vector. File paths to poly(A) data part
#'   files (as returned by \code{\link{split_polya_data}}).
#'
#' @param sequencing_summary Character string or data frame. Full path
#'   of the sequencing summary file, or an in-memory data frame.
#'
#' @param workspace Character string. Full path of the directory
#'   containing basecalled multi-Fast5 files.
#'
#' @param num_cores Numeric \code{[1]}. Number of physical cores.
#'
#' @param basecall_group Character string \code{["Basecall_1D_000"]}.
#'   Fast5 hierarchy level for data extraction.
#'
#' @param pass_only Logical \code{[TRUE]}. If \code{TRUE}, only
#'   \code{"PASS"} reads are included.
#'
#' @param qc Logical \code{[TRUE]}. If \code{TRUE}, terminal artefact
#'   positions are labelled with \code{"-WARN"}.
#'
#' @param save_dir Character string. Output directory path.
#'
#' @param prefix Character string (optional). Output file name prefix.
#'
#' @param cli_log Function. Logging closure defined in
#'   \code{\link{check_tails_guppy}}.
#'
#' @param ... Additional arguments passed to
#'   \code{\link{process_polya_complete}}.
#'
#' @return A named list with two data frames (merged across all parts):
#' \describe{
#'   \item{read_classes}{Per-read classification data.}
#'   \item{nonadenosine_residues}{Detailed non-A residue positional
#'     data.}
#' }
#'
#' @seealso \code{\link{check_tails_guppy}} which calls this function,
#'   \code{\link{split_polya_data}} for the input splitting step,
#'   \code{\link{process_polya_complete}} for the per-part processing.
#'
#' @keywords internal
#'
#' @export
process_polya_parts <- function(part_files,
                                sequencing_summary,
                                workspace,
                                num_cores,
                                basecall_group,
                                pass_only,
                                qc,
                                save_dir,
                                prefix,
                                cli_log,
                                ...) {

  # Initialize results containers
  all_read_classes <- list()
  all_nonadenosine_residues <- list()
  num_parts <- length(part_files)

  cli_log(sprintf("Starting sequential processing of %d poly(A) data parts", num_parts),
          "INFO", "Poly(A) Parts Processing", bullet = TRUE)

  #####################################################
  # LOOP FOR EACH FILE PART PROCESSING
  #####################################################
  for (i in seq_along(part_files)) {
    part_name <- sprintf("part_%d_of_%d", i, num_parts)
    part_dir <- file.path(save_dir, part_name)

    if (!dir.exists(part_dir)) {
      dir.create(part_dir, recursive = TRUE)
    }

    cli_log(sprintf("Processing part %d of %d", i, num_parts),
            "INFO", sprintf("Part %d Processing", i), bullet = TRUE)

    # Process each part (small file)
    part_output <- ninetails::process_polya_complete(
      polya_data = part_files[i],
      sequencing_summary = sequencing_summary,
      workspace = workspace,
      num_cores = num_cores,
      basecall_group = basecall_group,
      pass_only = pass_only,
      qc = qc,
      save_dir = part_dir,
      prefix = part_name,
      cli_log = cli_log)

    # Save part-specific outputs
    tryCatch({
      for (output_type in names(part_output)) {
        filename <- paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"),
                           '_', part_name, '_', output_type, '.tsv')
        filepath <- file.path(part_dir, filename)
        utils::write.table(part_output[[output_type]],
                           file = filepath,
                           row.names = FALSE,
                           sep = "\t",
                           quote = FALSE)
        cli_log(sprintf("Saved part output: %s", filename), "SUCCESS")
      }
    }, error = function(e) {
      cli_log(sprintf("Error saving part outputs: %s", e$message), "ERROR")
      stop(e$message, call. = FALSE)
    })

    # Collect results for final merge
    all_read_classes[[i]] <- part_output$read_classes
    all_nonadenosine_residues[[i]] <- part_output$nonadenosine_residues

    cli_log(sprintf("Poly(A) input file part %d processed successfully", i), "SUCCESS")
    cli_log(sprintf("Poly(A) input file part %d results saved in: %s", i, part_dir), "INFO")
  }

  # Merge results for final output
  cli_log("Merging results from all Poly(A) input file parts...", "INFO", "Final Processing")

  merged_output <- list(
    read_classes = do.call(rbind, all_read_classes),
    nonadenosine_residues = do.call(rbind, all_nonadenosine_residues)
  )

  # Log merge statistics
  cli_log(sprintf("Total reads processed: %d", nrow(merged_output$read_classes)), "INFO")
  cli_log(sprintf("Total residues identified: %d", nrow(merged_output$nonadenosine_residues)), "INFO")
  cli_log("All Poly(A) input file parts merged successfully", "SUCCESS")

  return(merged_output)
}
