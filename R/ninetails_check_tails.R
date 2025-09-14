
################################################################################
# WRAPPER FUNCTION FOR GUPPY/NANOPOLISH(TAILFINDR) OUTPUTS
################################################################################
#' Wrapper function for complete DRS processing by ninetails package (legacy mode).
#'
#' Important note: Due to updates from Oxford Nanopore Technologies
#' including the transition from the Guppy basecaller to Dorado,
#' a shift in the preferred native data format from FAST5 to POD5,
#' and a change in sequencing chemistry from R9 to R10,
#' the analysis pipeline based on Guppy outputs and Nanopolish/Tailfindr
#' poly(A) tail estimates has been moved to legacy mode.
#' To ensure backward compatibility with direct RNA sequencing (DRS) data
#' generated using earlier chemistries and algorithms,
#' this processing pipeline (previously implemented as
#' `check_tails()` in earlier versions of *ninetails*)
#' is now maintained as `check_tails_guppy()`.
#'
#' This function accepts either nanopolish or tailfindr outputs.
#'
#' This function allows to perform all of the steps required to discover
#' nonadenosine nucleotides within the given dataset using ninetails.
#' Please keep in mind, that during computations the function creates large
#' segmentation data. Therefore it may be wise to split the polya data table
#' beforehand and then run this function on table chunks.
#'
#' The output of this function is a list of 2 dataframes, containing:\itemize{
#' \item read_classes - classification of reads based on applied criteria
#' \item nonadenosine_residues - detailed positional info regarding all
#' potential nonadenosine residues detected.
#' }
#'
#' The more detailed info regarding read classification is stored within
#' 'comments' column. To make the output more compact, it contains codes
#' as follows:\itemize{
#' \item IRL - insufficient read length
#' \item QCF - nanopolish qc failed
#' \item MAU - move transition absent, nonA residue undetected
#' \item MPU - move transition present, nonA residue undetected
#' \item NIN - not included in the analysis (pass only = T)
#' \item YAY - move transition present, nonA residue detected
#' }
#'
#' The filtering criteria applied by ninetails are as follows: move of value =1
#' present, qc_tag = "PASS" or "PASS" & "SUFFCLIP", length estimated
#' by nanopolish >=10 nt).
#'
#' @param polya_data character string. Full path of the .tsv file produced
#' by either nanopolish polya function or tailfindr. Only DRS-derived tables are
#' accepted. In case of tailfindr output, it would be immediately converted to
#' nanopolish-like format.
#'
#' @param sequencing_summary character string. Full path of the .txt file
#' with sequencing summary.
#'
#' @param workspace character string. Full path of the directory to search the
#' basecalled fast5 files in. The Fast5 files have to be multi-fast5 file.
#'
#' @param num_cores numeric [1]. Number of physical cores to use in processing
#' the data. Do not exceed 1 less than the number of cores at your disposal.
#' This parameter is set to 1 by default.
#'
#' @param basecall_group character string ["Basecall_1D_000"]. Name of the
#' level in the Fast5 file hierarchy from which the data should be extracted.
#' This parameter is set to 'Basecall_1D_000' by default.
#'
#' @param pass_only logical [TRUE/FALSE]. If TRUE, only reads tagged by
#' nanopolish as "PASS" would be taken into consideration. Otherwise, reads
#' tagged as "PASS" & "SUFFCLIP" will be taken into account in analysis.
#' As a default, "TRUE" value is set.
#'
#' @param qc logical [TRUE/FALSE]. If TRUE, the quality control of the output
#' predictions would be performed. This means that the reads/non-A residue
#' positions in terminal nucleotides, which are most likely artifacts, are
#' labeled accordingly as "-WARN" (residues recognized as non-A due to
#' nanopolish segmentation error which is inherited from nanopolish,
#' as ninetails uses nanopolish segmentation). It is then up to user, whether
#' they would like to include or discard such reads from their pipeline. However,
#' it is advised to treat them with caution. By default, the qc option is enabled
#' (this parameter is set to TRUE).
#'
#' @param save_dir character string. Full path of the directory where the output
#' files containing the tail composition information should be stored.
#'
#' @param prefix character string (optional). If provided, it will be inserted
#' into the names of the output files, positioned between the timestamp
#' and the suffix indicating the file type (either read classes or
#' nonadenosine residues).
#'
#' @param part_size numeric [100000] (optional). If provided, defines maximum
#' number of rows in the poly(A) length containing input processed at once.
#' This value must be >=1000 (defaults to 1000000).If the size of
#' nanopolish/tailfindr output exceeds default 100000 rows, the input dataframe
#' would be splitted into parts and processed sequentially to avoid memory
#' overflow and pipeline crush.
#'
#' @export
#'
#' @return A list containing tail information organized by the read ID
#' is returned. Always assign this returned list to a variable, otherwise
#' the long list will be printed to the console, which may crash your R session.
#' Also a log file as well as read_classes & nonadenosine_residues files
#' are created in the user-specified directory.
#'
#' @importFrom foreach %dopar%
#' @importFrom utils head
#'
#' @examples
#' \dontrun{
#'
#'results <- ninetails::check_tails(
#'  polya_data = system.file('extdata',
#'                           'test_data',
#'                           'nanopolish_output.tsv',
#'                           package = 'ninetails'),
#'  sequencing_summary = system.file('extdata',
#'                                   'test_data',
#'                                   'sequencing_summary.txt',
#'                                   package = 'ninetails'),
#'  workspace = system.file('extdata',
#'                          'test_data',
#'                          'basecalled_fast5',
#'                          package = 'ninetails'),
#'  num_cores = 2,
#'  basecall_group = 'Basecall_1D_000',
#'  pass_only=TRUE,
#'  qc=TRUE,
#'  save_dir = '~/Downloads',
#'  prefix = "prefix",
#'  part_size=2000)
#'
#' }
check_tails <- function(polya_data,
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
    assertthat::assert_that(is.numeric(num_cores), num_cores > 0,
                            msg = "Number of cores must be a positive numeric value")

    assertthat::assert_that(is.character(workspace), dir.exists(workspace),
                            msg = "Fast5 files directory does not exist or path is invalid")

    assertthat::assert_that(checkmate::test_string(basecall_group, min.chars = 1),
                            msg = "Basecall group must be a non-empty character string")

    assertthat::assert_that(is.character(save_dir),
                            msg = "Output directory path must be a character string")

    assertthat::assert_that(is.logical(pass_only),
                            msg = "pass_only must be logical [TRUE/FALSE]")

    assertthat::assert_that(is.logical(qc),
                            msg = "qc must be logical [TRUE/FALSE]")

    # Validate prefix (if provided)
    if (nchar(prefix) > 0) {
      assertthat::assert_that(is.character(prefix),
                              msg = "File name prefix must be a character string")
    }

    # Validate nanopolish input
    if (checkmate::test_string(polya_data)) {
      assertthat::assert_that(file.exists(polya_data),
                              msg = "Poly(A) data file does not exist")
    } else {
      assertthat::assert_that(is.data.frame(polya_data) && nrow(polya_data) > 0,
                              msg = "Poly(A) data input must be a non-empty data frame or valid file path")
    }

    # Validate sequencing summary input
    if (checkmate::test_string(sequencing_summary)) {
      assertthat::assert_that(file.exists(sequencing_summary),
                              msg = "Sequencing summary file does not exist")
    } else {
      assertthat::assert_that(is.data.frame(sequencing_summary) && nrow(sequencing_summary) > 0,
                              msg = "Sequencing summary must be a non-empty data frame or valid file path")
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

#' Split large poly(A) data file into smaller polya_parts
#'
#' This function is not intended to be used as a standalone function, i.e., outside
#' the \code{\link{check_tails_guppy}} function. It does not affect data
#' processing; it only provides code to store results in files. It is enclosed
#' as a separate function, called by the wrapper, to keep the code clean and
#' modular, facilitating further development. Therefore, it should not be called
#' manually by the user.
#'
#' @param polya_data character string. Full path of the .tsv file produced
#' by either nanopolish polya function or tailfindr. Only DRS-derived tables are
#' accepted. In case of tailfindr output, it would be immediately converted to
#' nanopolish-like format.
#'
#' @param part_size numeric [100000] (optional). If provided, defines maximum
#' number of rows in the poly(A) length containing input processed at once.
#' This value must be >=1000 (defaults to 1000000).If the size of
#' nanopolish/tailfindr output exceeds default 100000 rows, the input dataframe
#' would be splitted into parts and processed sequentially to avoid memory
#' overflow and pipeline crush.
#'
#' @param save_dir character string. Full path of the directory where the output
#' files containing the tail composition information should be stored.
#'
#' @return List of file paths to created polya_parts
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



#' Save pipeline outputs to files
#'
#' This function is not intended to be used as a standalone function, i.e., outside
#' the \code{\link{check_tails_guppy}} function. It does not affect data
#' processing; it only provides code to store results in files. It is enclosed
#' as a separate function, called by the wrapper, to keep the code clean and
#' modular, facilitating further development. Therefore, it should not be called
#' manually by the user.
#'
#' @param outputs List containing read_classes and nonadenosine_residues
#' data frames produced by the ninetails pipeline
#'
#' @param save_dir character string. Full path of the directory where the output
#' files containing the tail composition information should be stored.
#'
#' @param prefix character string (optional). If provided, it will be inserted
#' into the names of the output files, positioned between the timestamp
#' and the suffix indicating the file type (either read classes or
#' nonadenosine residues).
#'
#' @return Invisible NULL
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


#' Variant of Check_tails_guppy processing of unsplitted poly(A) data.
#'
#' This variant of the pipeline allows to process small to moderate size of
#' poly(A) data (up to 100000 reads by default).
#'
#' This function as well as \code{\link{process_polya_parts}}function is not
#' intended to be used outside the \code{\link{check_tails_guppy}}function,
#' as it requires cli_log which is a custom function within the pipeline wrapper
#' required for log file & console prompts formatting. Therefore, should not be
#' manually called by the user.
#'
#' @param polya_data character string. Full path of the .tsv file produced
#' by either nanopolish polya function or tailfindr. Only DRS-derived tables are
#' accepted. In case of tailfindr output, it would be immediately converted to
#' nanopolish-like format.
#'
#' @param sequencing_summary character string. Full path of the .txt file
#' with sequencing summary.
#'
#' @param workspace character string. Full path of the directory to search the
#' basecalled fast5 files in. The Fast5 files have to be multi-fast5 file.
#'
#' @param num_cores numeric [1]. Number of physical cores to use in processing
#' the data. Do not exceed 1 less than the number of cores at your disposal.
#' This parameter is set to 1 by default.
#'
#' @param basecall_group character string ["Basecall_1D_000"]. Name of the
#' level in the Fast5 file hierarchy from which the data should be extracted.
#' This parameter is set to 'Basecall_1D_000' by default.
#'
#' @param pass_only logical [TRUE/FALSE]. If TRUE, only reads tagged by
#' nanopolish as "PASS" would be taken into consideration. Otherwise, reads
#' tagged as "PASS" & "SUFFCLIP" will be taken into account in analysis.
#' As a default, "TRUE" value is set.
#'
#' @param qc logical [TRUE/FALSE]. If TRUE, the quality control of the output
#' predictions would be performed. This means that the reads/non-A residue
#' positions in terminal nucleotides, which are most likely artifacts, are
#' labeled accordingly as "-WARN" (residues recognized as non-A due to
#' nanopolish segmentation error which is inherited from nanopolish,
#' as ninetails uses nanopolish segmentation). It is then up to user, whether
#' they would like to include or discard such reads from their pipeline. However,
#' it is advised to treat them with caution. By default, the qc option is enabled
#' (this parameter is set to TRUE).
#'
#' @param save_dir character string. Full path of the directory where the output
#' files containing the tail composition information should be stored.
#'
#' @param prefix character string (optional). If provided, it will be inserted
#' into the names of the output files, positioned between the timestamp
#' and the suffix indicating the file type (either read classes or
#' nonadenosine residues).
#'
#' @param cli_log Function for logging. This function is encoded in main
#' pipeline wrapper (in \code{\link{check_tails_guppy}} function). Its purpose
#' is to provide neatly formatted & informative log file.
#'
#' @return A list containing tail information organized by the read ID
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

#' Variant of check_tails_guppy processing of poly(A) data splitted to parts.
#'
#' In this pipeline variant, the large input poly(A) measurement file is split
#' to parts and processed sequentially one by one. This is to avoid memory
#' overflow and facilitate potential debugging.
#'
#' Creates individual outputs for each part and merges them into final results.
#'
#' This function as well as \code{\link{process_polya_complete}}function is not
#' intended to be used outside the \code{\link{check_tails_guppy}}function,
#' as it requires cli_log which is a custom function within the pipeline wrapper
#' required for log file & console prompts formatting. Therefore, should not be
#' manually called by the user.
#'
#' @param part_files List of file paths to poly(A) data parts
#'
#' @param sequencing_summary character string. Full path of the .txt file
#' with sequencing summary or in-memory object with sequencing summary data.
#'
#' @param workspace character string. Full path of the directory to search the
#' basecalled fast5 files in. The Fast5 files have to be multi-fast5 file.
#'
#' @param num_cores numeric [1]. Number of physical cores to use in processing
#' the data. Do not exceed 1 less than the number of cores at your disposal.
#' This parameter is set to 1 by default.
#'
#' @param basecall_group character string ["Basecall_1D_000"]. Name of the
#' level in the Fast5 file hierarchy from which the data should be extracted.
#' This parameter is set to 'Basecall_1D_000' by default.
#'
#' @param pass_only logical [TRUE/FALSE]. If TRUE, only reads tagged by
#' nanopolish as "PASS" would be taken into consideration. Otherwise, reads
#' tagged as "PASS" & "SUFFCLIP" will be taken into account in analysis.
#' As a default, "TRUE" value is set.
#'
#' @param qc logical [TRUE/FALSE]. If TRUE, the quality control of the output
#' predictions would be performed. This means that the reads/non-A residue
#' positions in terminal nucleotides, which are most likely artifacts, are
#' labeled accordingly as "-WARN" (residues recognized as non-A due to
#' nanopolish segmentation error which is inherited from nanopolish,
#' as ninetails uses nanopolish segmentation). It is then up to user, whether
#' they would like to include or discard such reads from their pipeline. However,
#' it is advised to treat them with caution. By default, the qc option is enabled
#' (this parameter is set to TRUE).
#'
#' @param save_dir character string. Full path of the directory where the output
#' files containing the tail composition information should be stored.
#'
#' @param prefix character string (optional). If provided, it will be inserted
#' into the names of the output files, positioned between the timestamp
#' and the suffix indicating the file type (either read classes or
#' nonadenosine residues).
#'
#' @param cli_log Function for logging
#'
#' @return A list containing tail information organized by the read ID.
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

