################################################################################
# WRAPPER FUNCTION FOR DORADO OUTPUTS
################################################################################
#' Wrapper function for complete DRS processing by ninetails package (Dorado mode).
#'
#' This function processes Oxford Nanopore direct RNA sequencing (DRS) data
#' generated using Dorado basecaller (>=1.0.0) and POD5 file format. It extracts
#' and analyzes poly(A) tail information from BAM files containing 'pt' and 'pa'
#' tags, which store poly(A) tail length and coordinates respectively.
#'
#' The function performs several key steps:\itemize{
#' \item Validates and preprocesses input files
#' \item Splits large BAM and summary files into manageable parts
#' \item Extracts poly(A) information from BAM files
#' \item Performs quality control if enabled
#' \item Generates comprehensive output files and logs
#' }
#'
#' The processing pipeline creates several subdirectories in the specified
#' output directory:\itemize{
#' \item dorado_summary_dir - Contains split summary files
#' \item bam_dir - Contains split BAM files
#' \item polya_data_dir - Contains extracted poly(A) information
#' }
#'
#' @param bam_file Character string. Full path to the BAM file produced by Dorado
#' (>=1.0.0) containing poly(A) information in 'pt' and 'pa' tags.
#'
#' @param dorado_summary Character string or data frame. Either the full path to
#' the Dorado summary file (.txt) or a data frame containing the summary
#' information with required columns (filename, read_id).
#'
#' @param pod5_dir Character string. Full path to the directory containing POD5
#' files. These files must correspond to the reads in the BAM file.
#'
#' @param num_cores Numeric [1]. Number of physical cores to use for processing.
#' Should not exceed one less than the total available cores. Defaults to 1.
#'
#' @param qc Logical [TRUE/FALSE]. If TRUE, performs quality control on the
#' output data. This includes checking for read duplicates and validating
#' poly(A) coordinates. Enabled by default.
#'
#' @param save_dir Character string. Full path to the directory where output
#' files and subdirectories will be created.
#'
#' @param prefix Character string (optional). When provided, this string is
#' prepended to output filenames, positioned between the timestamp and the
#' file type suffix.
#'
#' @param part_size Numeric [40000]. Defines the maximum number of reads to
#' process in each file part. Must be >=1000. Used to split large inputs
#' into manageable parts to prevent memory overflow.
#'
#' @return A list containing processed file paths and analysis results. Always
#' assign the returned list to a variable to prevent console overflow. The
#' function also generates a log file and output files in the specified
#' save_dir.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' results <- check_tails_dorado_DRS(
#'   bam_file = "path/to/basecalled.bam",
#'   dorado_summary = "path/to/alignment_summary.txt",
#'   pod5_dir = "path/to/pod5/",
#'   num_cores = 2,
#'   qc = TRUE,
#'   save_dir = "path/to/output/",
#'   prefix = "experiment1",
#'   part_size = 40000
#' )
#' }
check_tails_dorado_DRS <- function(bam_file,
                                   dorado_summary,
                                   pod5_dir,
                                   num_cores = 1,
                                   qc = TRUE,
                                   save_dir,
                                   prefix = "",
                                   part_size = 1000) {

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

    cli_log("Launching *check_tails_dorado_DRS* pipeline", "INFO","Launching *check_tails_dorado_DRS* pipeline")

    # Configuration logging display
    #####################################################
    cli_log("Configuration", "INFO", "Configuration")
    cli_log(sprintf("BAM file:                    %s", bam_file), bullet = TRUE)
    cli_log(sprintf("Dorado summary:              %s",
                    if(!is.object(dorado_summary)) dorado_summary else deparse(substitute(dorado_summary))), bullet = TRUE)
    cli_log(sprintf("Pod5 files directory:        %s", pod5_dir), bullet = TRUE)
    cli_log(sprintf("Number of cores:             %d", num_cores), bullet = TRUE)
    cli_log(sprintf("Output quality control:      %s", qc), bullet = TRUE)
    cli_log(sprintf("Output directory:            %s", save_dir), bullet = TRUE)
    cli_log(sprintf("Poly(A) processed at once:   %s", part_size), bullet = TRUE)

    # Preprocess input files
    processed_files <- preprocess_dorado_input(
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

    # Save final outputs
    cli_log("Saving final outputs...", "INFO", "Saving Results", bullet = TRUE)
    # ninetails::save_outputs(outputs, save_dir, prefix)

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
        "HERE SHOULD BE WARNING"
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
