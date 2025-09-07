################################################################################
# WRAPPER FUNCTION FOR DRS DORADO OUTPUTS (POD5/BAM)
################################################################################
#' Wrapper function for complete DRS processing by ninetails package (POD5/Dorado mode).
#'
#' This function processes direct RNA sequencing (DRS) data generated using
#' R10 chemistry SQK-RNA004 flowcells and analyzed with Dorado basecaller.
#' It accepts BAM files produced by Dorado's align function and POD5 files
#' containing raw signal data.
#'
#' This function allows to perform all of the steps required to discover
#' nonadenosine nucleotides within poly(A) tails using raw signal analysis.
#' Please keep in mind, that during computations the function creates large
#' segmentation data. Therefore it may be wise to process data in parts,
#' which is handled automatically by the function.
#'
#' The output of this function is a list of 2 dataframes, containing:\itemize{
#' \item read_classes - classification of reads based on applied criteria
#' \item nonadenosine_residues - detailed positional info regarding all
#' potential nonadenosine residues detected.
#' }
#'
#' @param bam_file character string. Full path to the BAM file produced by
#' Dorado's align function. The BAM file must be sorted and indexed.
#'
#' @param pod5_dir character string. Full path of the directory containing
#' POD5 files. The POD5 files must correspond to the reads in the BAM file.
#'
#' @param dorado_summary character string or data.frame. Either full path of the .txt file
#' with dorado summary or an in-memory data.frame containing dorado summary data.
#'
#' @param num_cores numeric [1]. Number of physical cores to use in processing
#' the data. Do not exceed 1 less than the number of cores at your disposal.
#' This parameter is set to 1 by default.
#'
#' @param save_dir character string. Full path of the directory where the output
#' files containing the tail composition information should be stored.
#'
#' @param prefix character string (optional). If provided, it will be inserted
#' into the names of the output files, positioned between the timestamp
#' and the suffix indicating the file type.
#'
#' @param part_size numeric [100000] (optional). If provided, defines maximum
#' number of reads to process in a single batch. Default is 100000 reads.
#'
#' @param qc logical [TRUE/FALSE]. If TRUE, quality control of nonadenosine
#' residue positions will be performed. As a default, "TRUE" value is set.
#'
#' @return A list containing tail information organized by the read ID.
#' Additionally, creates timestamped read_classes & nonadenosine_residues
#' files in the specified output directory.
#'
#' @export
#'
#' @examples
#'\dontrun{
#'
#' results <- ninetails::check_tails_dorado_DRS(
#'   bam_file = '/path/to/aligned.bam',
#'   pod5_dir = '/path/to/pod5/files',
#'   dorado_summary = '/path/to/summary.txt',
#'   num_cores = 2,
#'   save_dir = '~/Downloads',
#'   prefix = "prefix",
#'   part_size = 2000
#' )
#'
#'}
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
    #####################################################
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


