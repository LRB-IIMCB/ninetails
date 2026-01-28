
################################################################################
# BAM SPLITTING FUNCTION FOR cDNA PIPELINE
################################################################################

#' Split BAM file into parts based on read IDs from summary file
#'
#' This function splits a large BAM file into smaller parts based on read IDs
#' from corresponding Dorado summary files. This is essential for memory management
#' when processing large cDNA datasets. The function filters the BAM file to include
#' only reads present in the summary file and creates appropriately sized output files.
#'
#' @param bam_file Character string. Path to input BAM file to be split.
#' @param dorado_summary Character string. Path to corresponding Dorado summary file
#' containing read IDs to include in this part.
#' @param part_size Integer. Target number of reads per output file part.
#' @param save_dir Character string. Directory where split BAM files will be saved.
#' @param part_number Integer. Part number for naming output files.
#' @param cli_log Function for logging messages and progress.
#'
#' @return Character vector of output BAM file paths created.
#' @export
#'
#' @examples
#' \dontrun{
#' bam_files <- split_bam_file_cdna(
#'   bam_file = "large_dataset.bam",
#'   dorado_summary = "summary_part1.txt",
#'   part_size = 40000,
#'   save_dir = "bam_parts/",
#'   part_number = 1,
#'   cli_log = message
#' )
#' }
split_bam_file_cdna <- function(bam_file,
                                dorado_summary,
                                part_size = 100000,
                                save_dir,
                                part_number,
                                cli_log = message) {

  # Check for required Bioconductor packages
  if (!requireNamespace("Rsamtools", quietly = TRUE)) {
    stop("Package 'Rsamtools' is required for split_bam_file_cdna(). Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("S4Vectors", quietly = TRUE)) {
    stop("Package 'S4Vectors' is required for split_bam_file_cdna(). Please install it.",
         call. = FALSE)
  }

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

  if (missing(part_number)) {
    stop("Part number is missing. Please provide a valid part_number argument.",
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
  cli_log(sprintf("Reading summary file: %s", basename(dorado_summary)), "INFO")

  t1 <- Sys.time()
  # Get read IDs efficiently with vroom
  read_ids <- tryCatch({
    vroom::vroom(dorado_summary,
                 col_select = "read_id",
                 show_col_types = FALSE,
                 altrep = TRUE)
  }, error = function(e) {
    stop(sprintf("Error reading summary file %s: %s", basename(dorado_summary), e$message))
  })

  total_alignments <- nrow(read_ids)
  unique_reads <- length(unique(read_ids$read_id))

  cli_log(sprintf("Summary processed in %.2f seconds",
                  as.numeric(difftime(Sys.time(), t1, units = "secs"))), "INFO")

  cli_log(sprintf("Found %d alignments for %d unique reads",
                  total_alignments, unique_reads), "INFO")

  if(total_alignments > unique_reads) {
    cli_log(sprintf("Average alignments per read: %.2f",
                    total_alignments/unique_reads), "INFO")
  }

  # Get unique read IDs and clean up
  read_ids <- unique(read_ids$read_id)
  total_reads <- length(read_ids)
  gc()

  # Calculate how many output files we'll create
  num_parts <- ceiling(total_reads / part_size)

  # Calculate number of digits needed for file numbering
  num_digits <- nchar(as.character(num_parts))

  cli_log(sprintf("Splitting into %d files with ~%d reads each",
                  num_parts, part_size), "INFO")

  # Prepare output file paths with dynamic padding
  output_files <- sprintf(paste0("%s_%0", num_digits, "d.bam"),
                          output_prefix, seq_len(num_parts))

  # Split read IDs into parts of desired size
  read_parts <- split(read_ids, ceiling(seq_along(read_ids) / part_size))

  # Set up BAM reading parameters - we need read names for filtering
  param <- Rsamtools::ScanBamParam(what = "qname")

  # Process each part
  for (i in seq_along(read_parts)) {
    # Log progress
    cli_log(sprintf("Processing part %d/%d (expecting %d reads)",
                    i, num_parts, length(read_parts[[i]])), "INFO")

    # Prepare filter function for the current part of read IDs
    current_part <- read_parts[[i]]
    filter <- function(x) {
      matches <- x$qname %in% current_part
      return(matches)
    }

    # Filter BAM file to create part
    tryCatch({
      Rsamtools::filterBam(bam_file, output_files[i],
                           filter = S4Vectors::FilterRules(filter),
                           param = param)
    }, error = function(e) {
      stop(sprintf("Error filtering BAM file for part %d: %s", i, e$message))
    })

    # Count actual number of reads in output file
    actual_reads <- tryCatch({
      Rsamtools::countBam(output_files[i])$records
    }, error = function(e) {
      cli_log(sprintf("Warning: Could not count reads in %s", basename(output_files[i])), "WARNING")
      0
    })

    cli_log(sprintf("Created: %s (contains %d alignments)",
                    basename(output_files[i]), actual_reads), "INFO")

    if(actual_reads > length(current_part)) {
      cli_log(sprintf("Note: %d additional alignments included for these reads",
                      actual_reads - length(current_part)), "INFO")
    }
  }

  cli_log("BAM splitting completed successfully", "SUCCESS")

  return(output_files)
}

################################################################################
# DATA (SEQUENCE) EXTRACTION FUNCTION FROM BAM FILES
################################################################################

#' Extract data from BAM file for cDNA analysis
#'
#' This function extracts various informations from BAM files for reads that
#' are present in the corresponding Dorado summary file. It returns a data frame
#' with read IDs and their corresponding additional info, such as basecalled
#' sequences, poly(A) lengths, poly(A) coordinates etc., which can then be used
#' in downstream analyses.
#'
#' @param bam_file Character string. Path to BAM file containing basecalled sequences.
#' @param summary_file Character string. Path to corresponding Dorado summary file
#' containing read IDs to extract.
#' @param seq_only logical [TRUE]. When \code{TRUE}, only minimal information is
#' extracted, including read id, pod5 file name, and basecalled sequence.
#' If \code{FALSE}, a more comprehensive information is extracted, including poly(A)
#' tail length, coordinates etc.
#' @param cli_log Function for logging messages and progress.
#'
#' @return Data frame with columns:
#'   \describe{
#'     \item{read_id}{Character. Read identifier}
#'     \item{sequence}{Character. Basecalled sequence}
#'     \item{sequence_length}{Integer. Length of the sequence}
#'     \item{mapping_quality}{Integer. Mapping quality score}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' sequences <- extract_data_from_bam(
#'   bam_file = "aligned_reads.bam",
#'   summary_file = "dorado_summary.txt",
#'   seq_only = TRUE,
#'   cli_log = message
#' )
#' }
extract_data_from_bam <- function(bam_file, summary_file, seq_only = TRUE, cli_log = message) {

  # Check for required Bioconductor packages
  if (!requireNamespace("Rsamtools", quietly = TRUE)) {
    stop("Package 'Rsamtools' is required for extract_sequences_from_bam(). Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("S4Vectors", quietly = TRUE)) {
    stop("Package 'S4Vectors' is required for extract_sequences_from_bam(). Please install it.",
         call. = FALSE)
  }

  # Assertions
  if (!file.exists(bam_file)) {
    stop("BAM file does not exist: ", bam_file)
  }
  if (!file.exists(summary_file)) {
    stop("Summary file does not exist: ", summary_file)
  }

  # added parameter to filter only crucial info
  if (!is.logical(seq_only)) {
    stop("`seq_only` must be logical")
  }

  # Read summary file efficiently (only needed columns)
  cli_log(sprintf("Reading summary file: %s", base::basename(summary_file)), "INFO")
  summary_data <- vroom::vroom(summary_file,
                               col_select = c("filename", "read_id"),
                               show_col_types = FALSE)
  # Create lookup for pod5 files
  pod5_lookup <- summary_data$filename
  names(pod5_lookup) <- summary_data$read_id
  # Clear summary_data to free memory
  rm(summary_data)
  gc()

  # Create parameter settings for scanning BAM file
  param <- Rsamtools::ScanBamParam(
    tag = c("pt", "pa"),  # polya tags to extract
    what = c("qname", "rname", "pos", "mapq", "flag", "cigar", "seq"))

  # Process BAM file
  cli_log(sprintf("Reading BAM file: %s", base::basename(bam_file)), "INFO")
  bam_data <- Rsamtools::scanBam(bam_file, param = param)[[1]]
  initial_count <- length(bam_data$qname)
  cli_log(sprintf("Total read count: %d", initial_count), "INFO")

  # Check if pa tag exists in the BAM file
  if (is.null(bam_data$tag$pa) || all(sapply(bam_data$tag$pa, is.null))) {
    cli_log("poly(A) coordinate information (pa tag) not found in BAM file.", "ERROR")
    cli_log("Please ensure data were basecalled with dorado version >= 1.0.0", "ERROR")
    stop("poly(A) coordinate information (pa tag) not found in BAM file")
  }

  # Pre-allocate vectors with maximum possible size
  read_ids <- character(initial_count)
  references <- character(initial_count)
  mapqs <- integer(initial_count)
  poly_tail_lengths <- numeric(initial_count)
  pod5_files <- character(initial_count)  # New vector for pod5 files
  seqs <- character(initial_count) # basecalled sequences

  # conditionally for variant with more detailed info extraction
  if (seq_only == FALSE){
    ref_starts <- integer(initial_count)
    ref_ends <- integer(initial_count)
    anchor_positions <- integer(initial_count)
    poly_tail_starts <- integer(initial_count)
    polya_ends <- integer(initial_count)
    poly_tail2_starts <- integer(initial_count)
    poly_tail2_ends <- integer(initial_count)
  }

  # Counter for valid entries
  valid_count <- 0

  # Counters for logging
  mapped_count <- 0
  primary_count <- 0
  pt_tag_count <- 0

  # Process reads directly from BAM data
  cli_log("Processing reads...", "INFO")
  for (i in seq_along(bam_data$qname)) {
    # Skip if no pt tag or if it's NA
    if (is.null(bam_data$tag$pt[[i]]) || is.na(bam_data$tag$pt[[i]])) {
      next
    }

    # Filter unmapped reads
    if (is.na(bam_data$rname[i]) || bitwAnd(bam_data$flag[i], 0x4) || is.na(bam_data$pos[i])) {
      next
    }
    mapped_count <- mapped_count + 1

    # Filter non-primary alignments
    if (bitwAnd(bam_data$flag[i], 0x100) || bitwAnd(bam_data$flag[i], 0x800)) {
      next
    }
    primary_count <- primary_count + 1

    # Get polyA length from pt tag
    poly_tail_length <- bam_data$tag$pt[[i]]
    pt_tag_count <- pt_tag_count + 1

    # check whether extracting detailed info or not
    if (seq_only == FALSE){
      # Initialize pa tag values
      #anchor_pos <- -1L
      poly_tail_start <- -1L
      polya_end <- -1L
      poly_tail2_start <- -1L
      poly_tail2_end <- -1L

      # Add pa tag information if it exists and is valid
      if (!is.null(bam_data$tag$pa[[i]]) && length(bam_data$tag$pa[[i]]) == 5) {
        pa_values <- bam_data$tag$pa[[i]]
        if (!any(is.na(pa_values))) {
          #anchor_pos <- pa_values[1]
          poly_tail_start <- pa_values[2]
          polya_end <- pa_values[3]
          poly_tail2_start <- pa_values[4]
          poly_tail2_end <- pa_values[5]
        }
      }
    }

    # Get pod5 file from lookup
    current_read_id <- bam_data$qname[i]
    pod5_file <- pod5_lookup[current_read_id]
    if (is.na(pod5_file)) {
      next  # Skip if no matching pod5 file found
    }


    # Increment counter and store values
    valid_count <- valid_count + 1
    read_ids[valid_count] <- current_read_id
    references[valid_count] <- as.character(bam_data$rname[i])
    mapqs[valid_count] <- bam_data$mapq[i]
    pod5_files[valid_count] <- pod5_file
    seqs[valid_count] <- bam_data$seq[i] # basecalled sequences

    # optional values
    if (seq_only == FALSE){
      ref_starts[valid_count] <- bam_data$pos[i]
      ref_ends[valid_count] <- bam_data$pos[i] + nchar(bam_data$cigar[i])
      poly_tail_lengths[valid_count] <- poly_tail_length
      #anchor_positions[valid_count] <- anchor_pos
      poly_tail_starts[valid_count] <- poly_tail_start
      polya_ends[valid_count] <- polya_end
      poly_tail2_starts[valid_count] <- poly_tail2_start
      poly_tail2_ends[valid_count] <- poly_tail2_end
    }

  }

  # Output filtering statistics for logging
  cli_log(sprintf("Mapped reads filtered: %d", mapped_count), "INFO")
  cli_log(sprintf("Primary alignment filtered: %d", primary_count), "INFO")
  cli_log(sprintf("Reads with pt tag filtered: %d", pt_tag_count), "INFO")

  # Create data frame if we have any results
  if (valid_count > 0) {
    # Trim vectors to actual size
    read_ids <- read_ids[1:valid_count]
    references <- references[1:valid_count]
    mapqs <- mapqs[1:valid_count]
    pod5_files <- pod5_files[1:valid_count]
    seqs <- seqs[1:valid_count] # basecalled sequences

    if (seq_only == FALSE){
      ref_starts <- ref_starts[1:valid_count]
      ref_ends <- ref_ends[1:valid_count]
      poly_tail_lengths <- poly_tail_lengths[1:valid_count]
      #anchor_positions <- anchor_positions[1:valid_count]
      poly_tail_starts <- poly_tail_starts[1:valid_count]
      polya_ends <- polya_ends[1:valid_count]
      poly_tail2_starts <- poly_tail2_starts[1:valid_count]
      poly_tail2_ends <- poly_tail2_ends[1:valid_count]
    }


    # Check for duplicates
    duplicate_reads <- duplicated(read_ids)
    if (any(duplicate_reads)) {
      cli_log(sprintf("WARNING: Found %d duplicate read names, keeping first occurrence only",
                      sum(duplicate_reads)), "WARNING")
      # Keep only unique entries
      keep_idx <- !duplicate_reads
      read_ids <- read_ids[keep_idx]
      references <- references[keep_idx]
      pod5_files <- pod5_files[keep_idx]
      seqs <- seqs[keep_idx]
      # add more data categories if extended output is expected
      if (seq_only == FALSE){
        ref_starts <- ref_starts[keep_idx]
        ref_ends <- ref_ends[keep_idx]
        mapqs <- mapqs[keep_idx]
        poly_tail_lengths <- poly_tail_lengths[keep_idx]
        anchor_positions <- anchor_positions[keep_idx]
        poly_tail_starts <- poly_tail_starts[keep_idx]
        polya_ends <- polya_ends[keep_idx]
        poly_tail2_starts <- poly_tail2_starts[keep_idx]
        poly_tail2_ends <- poly_tail2_ends[keep_idx]
      }
    }

    # Create the final data frame
    if (seq_only == TRUE){
      polya_df <- tibble::tibble(
        read_id = read_ids,
        pod5_file = pod5_files,  # Added pod5 file information
        reference = references,
        sequence = seqs # bug fix: column name adjustment
      )
    } else if (seq_only == FALSE) {
      polya_df <- tibble::tibble(
        read_id = read_ids,
        pod5_file = pod5_files,  # Added pod5 file information
        reference = references,
        ref_start = ref_starts,
        ref_end = ref_ends,
        mapq = mapqs,
        poly_tail_length = poly_tail_lengths,
        #anchor_pos = anchor_positions,
        poly_tail_start = poly_tail_starts,
        polya_end = polya_ends,
        poly_tail2_start = poly_tail2_starts,
        poly_tail2_end = poly_tail2_ends,
        sequence = seqs
      )
    } else {
      cli_log("Wrong type of seq_only variable introduced (logical - TRUE/FALSE - required)", "ERROR")
      stop("seq_only is not logical")
    }



    cli_log("Processing complete", "SUCCESS")
    return(polya_df)
  } else {
    cli_log("No reads with extracted info found after filtering", "WARNING")
    return(tibble::tibble())
  }
}


################################################################################
# PREPROCESSING FUNCTION FOR cDNA INPUTS
################################################################################

#' Preprocess Dorado inputs for ninetails cDNA analysis
#'
#' This function prepares inputs for the cDNA pipeline by processing BAM files,
#' Dorado summary files, and extracting both basecalled sequences and poly(A)
#' signals from POD5 files. Unlike the DRS pipeline, this function additionally
#' extracts basecalled sequences from BAM files for read orientation classification.
#'
#' @param bam_file Character string. Path to BAM file containing aligned cDNA reads
#' with basecalled sequences.
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
#'     \item{bam_files}{Paths to split BAM files}
#'     \item{sequence_files}{Paths to extracted sequence files}
#'     \item{polya_signal_files}{Paths to extracted poly(A) signal files}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' processed_files <- preprocess_inputs_cdna(
#'   bam_file = "path/to/aligned.bam",
#'   dorado_summary = "path/to/summary.txt",
#'   pod5_dir = "path/to/pod5/",
#'   num_cores = 4,
#'   qc = TRUE,
#'   save_dir = "path/to/output/",
#'   prefix = "experiment1",
#'   part_size = 40000,
#'   cli_log = message
#' )
#' }
preprocess_inputs_cdna <- function(bam_file,
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
  } else {
    stop("BAM file must be provided as a character string path")
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


  # DORADO SUMMARY PROCESSING
  ################################################################################

  # Process dorado summary
  cli_log("Processing dorado summary...", "INFO", "Dorado Summary Processing", bullet = TRUE)

  # Create dorado summary directory
  summary_dir <- file.path(save_dir, "dorado_summary_dir")
  if (!dir.exists(summary_dir)) {
    dir.create(summary_dir, recursive = TRUE)
  }

  # Process summary file - handles size checking internally
  part_files <- ninetails::process_dorado_summary(
    dorado_summary = dorado_summary,
    save_dir = summary_dir,
    part_size = part_size,
    cli_log = cli_log
  )

  # Verify summary parts were created
  if (length(part_files) == 0) {
    stop("No summary parts were created")
  }


  # BAM FILE PROCESSING AND SPLITTING
  ################################################################################

  cli_log("Processing BAM file...", "INFO", "BAM Processing", bullet = TRUE)

  # Create BAM directory
  bam_dir <- file.path(save_dir, "bam_parts_dir")
  if (!dir.exists(bam_dir)) {
    dir.create(bam_dir, recursive = TRUE)
  }

  # Process BAM file based on whether summary was split
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

      bam_files[i] <- ninetails::split_bam_file_cdna(
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
    # Copy the original BAM file for consistency
    new_bam_path <- file.path(bam_dir, basename(bam_file))
    file.copy(from = bam_file, to = new_bam_path)
    bam_files <- new_bam_path

    cli_log("BAM file copied as single file", "SUCCESS")
  }


  # SEQUENCE EXTRACTION FROM BAM FILES
  ################################################################################

  cli_log("Extracting basecalled sequences from BAM files...", "INFO", "Sequence Extraction", bullet = TRUE)

  # Create sequence files directory
  sequence_dir <- file.path(save_dir, "sequence_files_dir")
  if (!dir.exists(sequence_dir)) {
    dir.create(sequence_dir, recursive = TRUE)
  }

  # Extract sequences from BAM files
  sequence_files <- character(length(bam_files))

  for (i in seq_along(bam_files)) {
    current_bam <- bam_files[i]
    current_summary <- part_files[i]

    cli_log(sprintf("Extracting sequences from BAM file %d/%d", i, length(bam_files)), "INFO")

    # Extract sequences from current BAM file
    sequence_data <- ninetails::extract_data_from_bam(
      bam_file = current_bam,
      summary_file = current_summary,
      seq_only = TRUE,
      cli_log = cli_log
    )

    # Save results if we have any
    if (nrow(sequence_data) > 0) {
      output_file <- file.path(sequence_dir,
                               sprintf("%ssequences_part%d.tsv",
                                       ifelse(nchar(prefix) > 0, paste0(prefix, "_"), ""),
                                       i))

      vroom::vroom_write(sequence_data, output_file, delim = "\t")
      sequence_files[i] <- output_file
      cli_log(sprintf("Saved sequences for BAM %d to %s (%d reads)",
                      i, basename(output_file), nrow(sequence_data)), "INFO")
    } else {
      cli_log(sprintf("No sequences extracted from BAM %d", i), "WARNING")
      sequence_files[i] <- ""  # Keep placeholder for consistency
    }

    # Explicitly trigger garbage collection after processing each BAM file
    gc()
  }

  # Filter out any empty results
  sequence_files <- sequence_files[file.exists(sequence_files)]

  if (length(sequence_files) > 0) {
    cli_log(sprintf("Successfully extracted sequences from %d BAM files", length(sequence_files)), "SUCCESS")
  } else {
    stop("No sequences were extracted from any BAM file")
  }


  # SIGNAL EXTRACTION FROM POD5 FILES
  ################################################################################

  cli_log("Extracting poly(A) signals from POD5 files...", "INFO", "Signal Extraction", bullet = TRUE)

  # Create polya_signal_dir for extracted signals
  polya_signal_dir <- file.path(save_dir, "polya_signal_dir")
  if (!dir.exists(polya_signal_dir)) {
    dir.create(polya_signal_dir, recursive = TRUE)
  }

  # Extract signals from POD5 files for each summary part
  polya_signal_files <- character(length(part_files))

  for (i in seq_along(part_files)) {
    cli_log(sprintf("Extracting signals from POD5 for summary part %d/%d", i, length(part_files)), "INFO")

    # Read summary data for this part
    summary_data <- vroom::vroom(part_files[i], show_col_types = FALSE)

    # Extract signals using adapted function
    signal_list <- ninetails::extract_tails_from_pod5(
      polya_data = summary_data,
      pod5_dir = pod5_dir
    )

    # Save extracted signals
    if (length(signal_list) > 0) {
      output_file <- file.path(polya_signal_dir,
                               sprintf("%spolya_signal_part%d.rds",
                                       ifelse(nchar(prefix) > 0, paste0(prefix, "_"), ""),
                                       i))

      saveRDS(signal_list, output_file)
      polya_signal_files[i] <- output_file
      cli_log(sprintf("Saved extracted signals for part %d to %s (%d signals)",
                      i, basename(output_file), length(signal_list)), "INFO")
    } else {
      cli_log(sprintf("No signals extracted for summary part %d", i), "WARNING")
      polya_signal_files[i] <- ""  # Keep placeholder
    }

    # Trigger garbage collection
    gc()
  }

  # Filter out empty results
  polya_signal_files <- polya_signal_files[file.exists(polya_signal_files)]

  if (length(polya_signal_files) > 0) {
    cli_log(sprintf("Successfully extracted signals for %d parts", length(polya_signal_files)), "SUCCESS")
  } else {
    cli_log("Warning: No signals were extracted from POD5 files", "WARNING")
  }

  # RETURN PROCESSED FILES
  ################################################################################

  cli_log("Preprocessing completed successfully", "SUCCESS")

  # Return organized file paths for downstream processing
  return(list(
    summary_files = part_files,
    bam_files = bam_files,
    sequence_files = sequence_files,
    polya_signal_files = polya_signal_files
  ))
}


################################################################################
# Detect orientation single
################################################################################
#' Detect poly tail type for a single sequence using Dorado-style algorithm
#'
#' This function implements the Dorado-style poly tail detection algorithm
#' for a single cDNA sequence. It uses edit distance matching of SSP and VNP
#' primer sequences at both ends of the read to determine orientation,
#' testing both forward (polyA) and reverse (polyT) strand configurations.
#'
#' @param sequence Character string containing the DNA sequence to classify
#'
#' @return Character string: "A" for polyA orientation, "T" for polyT orientation,
#'   or "unknown" for unclassified sequences
#' @export
#'
#' @examples
#' \dontrun{
#' detect_orientation_single("TTTCTGTTGGTGCTGATATTGCTTT...")  # Returns "A"
#' detect_orientation_single("ACTTGCCTGTCGCTCTATCTTCAG...")   # Returns "T"
#' }
detect_orientation_single <- function(sequence) {
  # Dorado primer sequences and parameters
  front_primer <- "TTTCTGTTGGTGCTGATATTGCTTT"  # SSP from dorado
  rear_primer <- "ACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAGTTTT"  # VNP from dorado
  primer_window <- 150
  flank_threshold <- 0.6
  min_primer_separation <- 10

  # Handle empty or very short sequences
  if (is.na(sequence) || nchar(sequence) < 50) {
    return("unknown")
  }

  tryCatch({
    # Trim trailing Ts from rear_primer
    trailing_Ts <- count_trailing_chars(rear_primer, 'T')
    rear_primer_trimmed <- if (trailing_Ts > 0) {
      substr(rear_primer, 1, nchar(rear_primer) - trailing_Ts)
    } else {
      rear_primer
    }

    # Generate reverse complements
    front_primer_rc <- reverse_complement(front_primer)
    rear_primer_rc <- reverse_complement(rear_primer_trimmed)

    # Extract windows
    seq_len <- nchar(sequence)
    read_top <- substr(sequence, 1, min(primer_window, seq_len))
    bottom_start <- max(1, seq_len - primer_window + 1)
    read_bottom <- substr(sequence, bottom_start, seq_len)

    # Test forward strand (poly-A)
    dist_v1 <- edit_distance_hw(front_primer, read_top) +
      edit_distance_hw(rear_primer_rc, read_bottom)

    # Test reverse strand (poly-T)
    dist_v2 <- edit_distance_hw(rear_primer_trimmed, read_top) +
      edit_distance_hw(front_primer_rc, read_bottom)

    # Determine winner
    fwd <- dist_v1 < dist_v2

    # Calculate score
    total_primer_length <- nchar(front_primer) + nchar(rear_primer_trimmed)
    flank_score <- 1.0 - (min(dist_v1, dist_v2) / total_primer_length)

    # Validate
    score_passes <- flank_score >= flank_threshold
    separation_passes <- abs(dist_v1 - dist_v2) > min_primer_separation
    proceed <- score_passes && separation_passes

    # Return tail type
    if (!proceed) {
      return("unknown")
    } else if (fwd) {
      return("A")  # Forward strand -> poly-A
    } else {
      return("T")  # Reverse strand -> poly-T
    }

  }, error = function(e) {
    return("unknown")
  })
}



################################################################################
# Detect orientation multiple
################################################################################
#' Classify multiple cDNA read orientations using Dorado-style poly tail detection
#'
#' This function analyzes Nanopore cDNA sequences using the Dorado-style approach
#' for poly tail type detection. It uses edit distance-based matching of SSP and VNP
#' primer sequences at both ends of reads, testing forward and reverse orientations
#' to determine if reads contain poly(A) or poly(T) tails. The method uses sliding
#' window matching and requires both score and separation thresholds to be met.
#' The tail_type column is added to original sequence files and all sequences
#' are returned as a single tibble.
#'
#' @param sequence_files Character vector. Paths to sequence files (TSV format)
#' containing read_id and sequence columns.
#' @param num_cores Integer. Number of CPU cores to use for parallel processing.
#' @param cli_log Function for logging messages and progress.
#'
#' @return Data frame/tibble containing all sequences with added tail_type column
#' ("polyA", "polyT", or "unidentified") based on Dorado-style edit distance matching
#' of SSP and VNP primers with score and separation validation.
#' The tail_type column is also written back to the original sequence files for persistence.
#' @export
#'
#' @examples
#' \dontrun{
#' classified_sequences <- detect_orientation_multiple(
#'   sequence_files = c("sequences_part1.tsv", "sequences_part2.tsv"),
#'   num_cores = 4,
#'   cli_log = message
#' )
#' }
detect_orientation_multiple <- function(sequence_files,
                                        num_cores = 1,
                                        cli_log = message) {

  # Input validation
  if (missing(sequence_files)) {
    stop("Sequence files are missing. Please provide valid sequence_files argument.",
         call. = FALSE)
  }

  assertthat::assert_that(is.character(sequence_files),
                          msg = "sequence_files must be a character vector of file paths")
  assertthat::assert_that(all(file.exists(sequence_files)),
                          msg = "All sequence files must exist")
  assertthat::assert_that(is.numeric(num_cores), num_cores > 0,
                          msg = "Number of cores must be a positive numeric value")

  cli_log("Starting read orientation classification", "INFO", "Sequence Classification")

  # Initialize collection for all sequences
  all_classified_sequences <- tibble::tibble()

  # Initialize statistics
  total_polya <- 0
  total_polyt <- 0
  total_unidentified <- 0

  # Process each sequence file
  ################################################################################
  for (i in seq_along(sequence_files)) {
    current_file <- sequence_files[i]

    cli_log(sprintf("Processing sequence file %d/%d: %s",
                    i, length(sequence_files), basename(current_file)), "INFO")

    tryCatch({
      # Read sequence file
      sequence_data <- vroom::vroom(current_file, show_col_types = FALSE)

      if (!"sequence" %in% colnames(sequence_data)) {
        cli_log(sprintf("Warning: No 'sequence' column found in %s", basename(current_file)), "WARNING")
        next
      }

      cli_log(sprintf("Loaded %d sequences for classification", nrow(sequence_data)), "INFO")

      # Apply Dorado-style classification to all sequences
      cli_log("Applying Dorado-style poly tail classification...", "INFO")

      # Apply classification and map results to ninetails naming convention
      raw_classifications <- sapply(sequence_data$sequence, ninetails::detect_orientation_single)

      # Map Dorado results to ninetails naming
      classifications <- sapply(raw_classifications, function(result) {
        switch(result,
               "A" = "polyA",
               "T" = "polyT",
               "unknown" = "unidentified",
               "unidentified")  # Default fallback
      })

      # Add tail_type column
      sequence_data$tail_type <- classifications

      # Update statistics
      file_polya <- sum(classifications == "polyA")
      file_polyt <- sum(classifications == "polyT")
      file_unidentified <- sum(classifications == "unidentified")

      total_polya <- total_polya + file_polya
      total_polyt <- total_polyt + file_polyt
      total_unidentified <- total_unidentified + file_unidentified

      cli_log(sprintf("Classification complete: %d polyA, %d polyT, %d unidentified",
                      file_polya, file_polyt, file_unidentified), "INFO")

      # Write back to original file with tail_type column
      vroom::vroom_write(sequence_data, current_file, delim = "\t")
      cli_log(sprintf("Updated %s with tail_type column", basename(current_file)), "INFO")

      # Add to collection
      all_classified_sequences <- dplyr::bind_rows(all_classified_sequences, sequence_data)

    }, error = function(e) {
      cli_log(sprintf("Error processing file %s: %s", basename(current_file), e$message), "ERROR")
    })

    # Trigger garbage collection
    gc()
  }

  # Create classification statistics from collected data
  classification_stats <- list(
    total_reads = nrow(all_classified_sequences),
    polya_reads = total_polya,
    polyt_reads = total_polyt,
    unidentified_reads = total_unidentified
  )

  # Add percentages
  if (classification_stats$total_reads > 0) {
    classification_stats$polya_percentage <- round(classification_stats$polya_reads / classification_stats$total_reads * 100, 2)
    classification_stats$polyt_percentage <- round(classification_stats$polyt_reads / classification_stats$total_reads * 100, 2)
    classification_stats$unidentified_percentage <- round(classification_stats$unidentified_reads / classification_stats$total_reads * 100, 2)
  } else {
    classification_stats$polya_percentage <- 0
    classification_stats$polyt_percentage <- 0
    classification_stats$unidentified_percentage <- 0
  }

  # Final logging
  ################################################################################
  cli_log("Classification Summary", "INFO", "Classification Summary")
  cli_log(sprintf("Total reads processed: %d", classification_stats$total_reads), bullet = TRUE)
  cli_log(sprintf("PolyA reads: %d (%.2f%%)",
                  classification_stats$polya_reads, classification_stats$polya_percentage), bullet = TRUE)
  cli_log(sprintf("PolyT reads: %d (%.2f%%)",
                  classification_stats$polyt_reads, classification_stats$polyt_percentage), bullet = TRUE)
  cli_log(sprintf("Unidentified reads: %d (%.2f%%)",
                  classification_stats$unidentified_reads, classification_stats$unidentified_percentage), bullet = TRUE)

  cli_log("Read orientation classification completed", "SUCCESS")

  # Return the classified sequences tibble
  return(all_classified_sequences)
}



################################################################################
# POLYA READS PROCESSING FUNCTION FOR cDNA (UPDATED)
################################################################################

#' Process polyA reads using standard ninetails pipeline
#'
#' This function processes reads that have been classified as polyA-containing
#' through the cDNA classification pipeline. Updated to work with tibble input
#' instead of separate files. It applies the standard ninetails analysis
#' pipeline to identify non-adenosine residues within polyA tails.
#'
#' @param polya_sequences Data frame/tibble. PolyA-classified sequences with
#' read_id and tail_type columns.
#' @param signal_files Character vector. Paths to signal files containing
#' poly(A) tail signals extracted from POD5 files.
#' @param num_cores Integer. Number of CPU cores to use for parallel processing.
#' @param qc Logical. Whether to apply quality control filtering.
#' @param save_dir Character. Directory where processing results will be saved.
#' @param prefix Character. Optional prefix for output file names.
#' @param cli_log Function for logging messages and progress.
#'
#' @return List containing polyA processing results:
#'   \describe{
#'     \item{read_classes}{Data frame with read classification results (with tail_type preserved)}
#'     \item{nonadenosine_residues}{Data frame with predicted modifications (with tail_type preserved)}
#'     \item{processing_stats}{Summary statistics for polyA processing}
#'     \item{tail_type}{Character indicating this is "polyA" data}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' polya_results <- process_polya_reads_cdna(
#'   polya_sequences = polya_tibble,
#'   signal_files = c("polya_signal_part1.rds"),
#'   num_cores = 4,
#'   qc = TRUE,
#'   save_dir = "path/to/output/",
#'   prefix = "experiment1",
#'   cli_log = message
#' )
#' }
process_polya_reads_cdna <- function(polya_sequences,
                                     signal_files,
                                     num_cores = 1,
                                     qc = TRUE,
                                     save_dir,
                                     prefix = "",
                                     cli_log = message) {

  # Input validation
  if (missing(polya_sequences)) {
    stop("PolyA sequences are missing. Please provide valid polya_sequences argument.",
         call. = FALSE)
  }

  if (missing(signal_files)) {
    stop("Signal files are missing. Please provide valid signal_files argument.",
         call. = FALSE)
  }

  if (missing(save_dir)) {
    stop("Output directory is missing. Please provide a valid save_dir argument.",
         call. = FALSE)
  }

  assertthat::assert_that(is.data.frame(polya_sequences),
                          msg = "polya_sequences must be a data frame")
  assertthat::assert_that("read_id" %in% colnames(polya_sequences),
                          msg = "polya_sequences must contain a 'read_id' column")
  assertthat::assert_that(is.character(signal_files),
                          msg = "signal_files must be a character vector of file paths")
  assertthat::assert_that(all(file.exists(signal_files)),
                          msg = "All signal files must exist")

  cli_log("Starting polyA reads processing", "INFO", "PolyA Processing")

  # Create output directories for polyA processing
  polya_processing_dir <- file.path(save_dir, "polya_processing_dir")
  polya_temp_dir <- file.path(polya_processing_dir, "polya_temp_dir")

  for (dir in c(polya_processing_dir, polya_temp_dir)) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
  }

  ################################################################################
  # FILTER SIGNALS FOR POLYA READS
  ################################################################################

  cli_log("Filtering signals for polyA reads...", "INFO", bullet = TRUE)

  # Get polyA read IDs
  polya_read_ids <- polya_sequences$read_id
  cli_log(sprintf("Found %d polyA read IDs for processing", length(polya_read_ids)), "INFO")

  # Load and filter signal files to keep only polyA reads
  all_polya_signals <- list()

  for (i in seq_along(signal_files)) {
    signal_file <- signal_files[i]
    cli_log(sprintf("Filtering signal file %d/%d: %s",
                    i, length(signal_files), basename(signal_file)), "INFO")

    tryCatch({
      # Load signal data
      signal_list <- readRDS(signal_file)

      # Filter to keep only polyA reads
      polya_signals <- signal_list[names(signal_list) %in% polya_read_ids]

      if (length(polya_signals) > 0) {
        all_polya_signals <- c(all_polya_signals, polya_signals)
        cli_log(sprintf("Added %d polyA signals from %s",
                        length(polya_signals), basename(signal_file)), "INFO")
      }

    }, error = function(e) {
      cli_log(sprintf("Error filtering signal file %s: %s", basename(signal_file), e$message), "ERROR")
    })
  }

  if (length(all_polya_signals) == 0) {
    cli_log("No polyA signals were found for processing", "WARNING")
    return(list(
      read_classes = data.frame(),
      nonadenosine_residues = data.frame(),
      processing_stats = list(total_reads_processed = 0),
      tail_type = "polyA"
    ))
  }

  cli_log(sprintf("Total polyA signals loaded: %d", length(all_polya_signals)), "SUCCESS")

  ################################################################################
  # NINETAILS PROCESSING PIPELINE FOR POLYA (DRS-STYLE)
  ################################################################################

  # Create intermediate directories following DRS pattern
  nonA_temp_dir <- file.path(polya_processing_dir, "nonA_temp_dir")
  polya_chunks_dir <- file.path(polya_processing_dir, "polya_chunks_dir")
  dorado_summary_dir <- file.path(polya_processing_dir, "dorado_summary_dir")

  for (dir in c(nonA_temp_dir, polya_chunks_dir, dorado_summary_dir)) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
  }

  # Step 1: Create and save summary file (following DRS pattern)
  cli_log("Creating summary file for polyA reads...", "INFO", bullet = TRUE)

  temp_summary_data <- data.frame(
    read_id = names(all_polya_signals),
    readname = names(all_polya_signals),
    contig = "polya_tail",
    poly_tail_length = sapply(all_polya_signals, length),
    poly_tail_start = 1,  # Dummy values for compatibility
    poly_tail_end = sapply(all_polya_signals, length),
    alignment_genome = "polya_tail",
    alignment_mapq = 60,  # Good quality score
    stringsAsFactors = FALSE
  )

  # Save summary file to directory
  summary_file <- file.path(dorado_summary_dir, "polya_summary_part1.txt")
  vroom::vroom_write(temp_summary_data, summary_file, delim = "\t")

  cli_log(sprintf("Created summary file with %d polyA reads", nrow(temp_summary_data)), "INFO")

  # Step 2: Process signals using DRS-style approach
  cli_log("Processing polyA signals through ninetails pipeline...", "INFO", bullet = TRUE)

  # Save signals to RDS file for processing
  signal_file <- file.path(polya_temp_dir, "polya_signal_part1.rds")
  saveRDS(all_polya_signals, signal_file)

  # Use the DRS signal processing function
  tryCatch({
    ninetails::process_dorado_signal_files(
      polya_signal_files = signal_file,  # Note: parameter name is polya_signal_files in the function
      nonA_temp_dir = nonA_temp_dir,
      polya_chunks_dir = polya_chunks_dir,
      num_cores = num_cores,
      cli_log = cli_log
    )
  }, error = function(e) {
    cli_log(sprintf("Error in polyA signal processing: %s", e$message), "ERROR")
    stop(e$message, call. = FALSE)
  })

  # Step 3: Create outputs using cDNA-compatible function
  cli_log("Creating final polyA outputs...", "INFO", bullet = TRUE)

  outputs <- tryCatch({
    result <- ninetails::create_outputs_dorado_cdna(
      dorado_summary_dir = dorado_summary_dir,
      nonA_temp_dir = nonA_temp_dir,
      polya_chunks_dir = polya_chunks_dir,
      num_cores = num_cores,
      qc = qc
    )
    result
  }, error = function(e) {
    cli_log(sprintf("Error in polyA output creation: %s", e$message), "ERROR")
    stop(e$message, call. = FALSE)
  })

  # Cleanup signal file
  if (file.exists(signal_file)) {
    file.remove(signal_file)
  }

  cli_log("PolyA processing completed using cDNA-compatible pipeline", "SUCCESS")

  cli_log("PolyA output creation completed", "SUCCESS")

  ################################################################################
  # COMPILE FINAL RESULTS
  ################################################################################

  # Calculate processing statistics
  processing_stats <- list(
    total_reads_processed = if (!is.null(outputs$read_classes)) nrow(outputs$read_classes) else 0,
    decorated_reads = if (!is.null(outputs$read_classes)) sum(outputs$read_classes$class == "decorated", na.rm = TRUE) else 0,
    blank_reads = if (!is.null(outputs$read_classes)) sum(outputs$read_classes$class == "blank", na.rm = TRUE) else 0,
    unclassified_reads = if (!is.null(outputs$read_classes)) sum(outputs$read_classes$class == "unclassified", na.rm = TRUE) else 0,
    total_modifications = if (!is.null(outputs$nonadenosine_residues)) nrow(outputs$nonadenosine_residues) else 0,
    tail_type = "polyA"
  )

  cli_log("PolyA Processing Summary", "INFO", "PolyA Summary")
  cli_log(sprintf("Total polyA reads processed: %d", processing_stats$total_reads_processed), bullet = TRUE)
  cli_log(sprintf("Decorated reads: %d", processing_stats$decorated_reads), bullet = TRUE)
  cli_log(sprintf("Blank reads: %d", processing_stats$blank_reads), bullet = TRUE)
  cli_log(sprintf("Unclassified reads: %d", processing_stats$unclassified_reads), bullet = TRUE)
  cli_log(sprintf("Total modifications detected: %d", processing_stats$total_modifications), bullet = TRUE)

  cli_log("PolyA reads processing completed successfully", "SUCCESS")

  # Return combined results with tail_type
  return(list(
    read_classes = outputs$read_classes,
    nonadenosine_residues = outputs$nonadenosine_residues,
    processing_stats = processing_stats,
    tail_type = "polyA"
  ))
}



################################################################################
# POLYT READS PROCESSING FUNCTION FOR cDNA (UPDATED)
################################################################################

#' Process polyT reads using ninetails pipeline
#'
#' This function processes reads that have been classified as polyT-containing
#' through the cDNA classification pipeline. Updated to work with tibble input
#' instead of separate files. It applies the ninetails analysis pipeline to
#' identify non-thymidine residues within polyT tails. Currently uses the same
#' model as polyA processing until a polyT-specific model is trained.
#'
#' @param polyt_sequences Data frame/tibble. PolyT-classified sequences with
#' read_id and tail_type columns.
#' @param signal_files Character vector. Paths to signal files containing
#' poly(A) tail signals extracted from POD5 files (will be filtered for polyT reads).
#' @param num_cores Integer. Number of CPU cores to use for parallel processing.
#' @param qc Logical. Whether to apply quality control filtering.
#' @param save_dir Character. Directory where processing results will be saved.
#' @param prefix Character. Optional prefix for output file names.
#' @param cli_log Function for logging messages and progress.
#'
#' @return List containing polyT processing results:
#'   \describe{
#'     \item{read_classes}{Data frame with read classification results (with tail_type preserved)}
#'     \item{nonadenosine_residues}{Data frame with predicted modifications (with tail_type preserved)}
#'     \item{processing_stats}{Summary statistics for polyT processing}
#'     \item{tail_type}{Character indicating this is "polyT" data}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' polyt_results <- process_polyt_reads_cdna(
#'   polyt_sequences = polyt_tibble,
#'   signal_files = c("polya_signal_part1.rds"),
#'   num_cores = 4,
#'   qc = TRUE,
#'   save_dir = "path/to/output/",
#'   prefix = "experiment1",
#'   cli_log = message
#' )
#' }
process_polyt_reads_cdna <- function(polyt_sequences,
                                     signal_files,
                                     num_cores = 1,
                                     qc = TRUE,
                                     save_dir,
                                     prefix = "",
                                     cli_log = message) {

  # Input validation
  if (missing(polyt_sequences)) {
    stop("PolyT sequences are missing. Please provide valid polyt_sequences argument.",
         call. = FALSE)
  }

  if (missing(signal_files)) {
    stop("Signal files are missing. Please provide valid signal_files argument.",
         call. = FALSE)
  }

  if (missing(save_dir)) {
    stop("Output directory is missing. Please provide a valid save_dir argument.",
         call. = FALSE)
  }

  assertthat::assert_that(is.data.frame(polyt_sequences),
                          msg = "polyt_sequences must be a data frame")
  assertthat::assert_that("read_id" %in% colnames(polyt_sequences),
                          msg = "polyt_sequences must contain a 'read_id' column")
  assertthat::assert_that(is.character(signal_files),
                          msg = "signal_files must be a character vector of file paths")
  assertthat::assert_that(all(file.exists(signal_files)),
                          msg = "All signal files must exist")

  cli_log("Starting polyT reads processing", "INFO", "PolyT Processing")

  # NOTE: Currently using same model as polyA until polyT-specific model is trained
  cli_log("NOTE: Currently using polyA model for polyT reads (temporary)", "INFO", bullet = TRUE)

  # Create output directories for polyT processing
  polyt_processing_dir <- file.path(save_dir, "polyt_processing_dir")
  polyt_temp_dir <- file.path(polyt_processing_dir, "polyt_temp_dir")

  for (dir in c(polyt_processing_dir, polyt_temp_dir)) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
  }

  ################################################################################
  # FILTER SIGNALS FOR POLYT READS
  ################################################################################

  cli_log("Filtering signals for polyT reads...", "INFO", bullet = TRUE)

  # Get polyT read IDs
  polyt_read_ids <- polyt_sequences$read_id
  cli_log(sprintf("Found %d polyT read IDs for processing", length(polyt_read_ids)), "INFO")

  # Load and filter signal files to keep only polyT reads
  all_polyt_signals <- list()

  for (i in seq_along(signal_files)) {
    signal_file <- signal_files[i]
    cli_log(sprintf("Filtering signal file %d/%d: %s",
                    i, length(signal_files), basename(signal_file)), "INFO")

    tryCatch({
      # Load signal data
      signal_list <- readRDS(signal_file)

      # Filter to keep only polyT reads
      polyt_signals <- signal_list[names(signal_list) %in% polyt_read_ids]

      if (length(polyt_signals) > 0) {
        all_polyt_signals <- c(all_polyt_signals, polyt_signals)
        cli_log(sprintf("Added %d polyT signals from %s",
                        length(polyt_signals), basename(signal_file)), "INFO")
      }

    }, error = function(e) {
      cli_log(sprintf("Error filtering signal file %s: %s", basename(signal_file), e$message), "ERROR")
    })
  }

  if (length(all_polyt_signals) == 0) {
    cli_log("No polyT signals were found for processing", "WARNING")
    return(list(
      read_classes = data.frame(),
      nonadenosine_residues = data.frame(),
      processing_stats = list(total_reads_processed = 0),
      tail_type = "polyT"
    ))
  }

  cli_log(sprintf("Total polyT signals loaded: %d", length(all_polyt_signals)), "SUCCESS")

  ################################################################################
  # NINETAILS PROCESSING PIPELINE FOR POLYT (DRS-STYLE)
  ################################################################################

  # Create intermediate directories following DRS pattern
  nonA_temp_dir <- file.path(polyt_processing_dir, "nonA_temp_dir")
  polyt_chunks_dir <- file.path(polyt_processing_dir, "polyt_chunks_dir")
  dorado_summary_dir <- file.path(polyt_processing_dir, "dorado_summary_dir")

  for (dir in c(nonA_temp_dir, polyt_chunks_dir, dorado_summary_dir)) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
  }

  # Step 1: Create and save summary file (following DRS pattern)
  cli_log("Creating summary file for polyT reads...", "INFO", bullet = TRUE)

  temp_summary_data <- data.frame(
    read_id = names(all_polyt_signals),
    readname = names(all_polyt_signals),
    contig = "polyt_tail",
    poly_tail_length = sapply(all_polyt_signals, length),
    poly_tail_start = 1,  # Dummy values for compatibility
    poly_tail_end = sapply(all_polyt_signals, length),
    alignment_genome = "polyt_tail",
    alignment_mapq = 60,  # Good quality score
    stringsAsFactors = FALSE
  )

  # Save summary file to directory
  summary_file <- file.path(dorado_summary_dir, "polyt_summary_part1.txt")
  vroom::vroom_write(temp_summary_data, summary_file, delim = "\t")

  cli_log(sprintf("Created summary file with %d polyT reads", nrow(temp_summary_data)), "INFO")

  # Step 2: Process signals using DRS-style approach
  cli_log("Processing polyT signals through ninetails pipeline...", "INFO", bullet = TRUE)
  cli_log("NOTE: Using polyA model for polyT reads (temporary)", "INFO", bullet = TRUE)

  # Save signals to RDS file for processing
  signal_file <- file.path(polyt_temp_dir, "polyt_signal_part1.rds")
  saveRDS(all_polyt_signals, signal_file)

  # Use the DRS signal processing function (temporarily with polyA model)
  tryCatch({
    ninetails::process_dorado_signal_files(
      polya_signal_files = signal_file,  # Note: parameter name is polya_signal_files in the function
      nonA_temp_dir = nonA_temp_dir,
      polya_chunks_dir = polyt_chunks_dir,
      num_cores = num_cores,
      cli_log = cli_log
    )
  }, error = function(e) {
    cli_log(sprintf("Error in polyT signal processing: %s", e$message), "ERROR")
    stop(e$message, call. = FALSE)
  })

  # Step 3: Create outputs using cDNA-compatible function
  cli_log("Creating final polyT outputs...", "INFO", bullet = TRUE)

  outputs <- tryCatch({
    result <- ninetails::create_outputs_dorado_cdna(
      dorado_summary_dir = dorado_summary_dir,
      nonA_temp_dir = nonA_temp_dir,
      polya_chunks_dir = polyt_chunks_dir,
      num_cores = num_cores,
      qc = qc
    )
    result
  }, error = function(e) {
    cli_log(sprintf("Error in polyT output creation: %s", e$message), "ERROR")
    stop(e$message, call. = FALSE)
  })

  # Cleanup signal file
  if (file.exists(signal_file)) {
    file.remove(signal_file)
  }

  cli_log("PolyT processing completed using cDNA-compatible pipeline", "SUCCESS")

  cli_log("PolyT output creation completed", "SUCCESS")

  ################################################################################
  # COMPILE FINAL RESULTS
  ################################################################################

  # Calculate processing statistics
  processing_stats <- list(
    total_reads_processed = if (!is.null(outputs$read_classes)) nrow(outputs$read_classes) else 0,
    decorated_reads = if (!is.null(outputs$read_classes)) sum(outputs$read_classes$class == "decorated", na.rm = TRUE) else 0,
    blank_reads = if (!is.null(outputs$read_classes)) sum(outputs$read_classes$class == "blank", na.rm = TRUE) else 0,
    unclassified_reads = if (!is.null(outputs$read_classes)) sum(outputs$read_classes$class == "unclassified", na.rm = TRUE) else 0,
    total_modifications = if (!is.null(outputs$nonadenosine_residues)) nrow(outputs$nonadenosine_residues) else 0,
    tail_type = "polyT"
  )

  cli_log("PolyT Processing Summary", "INFO", "PolyT Summary")
  cli_log(sprintf("Total polyT reads processed: %d", processing_stats$total_reads_processed), bullet = TRUE)
  cli_log(sprintf("Decorated reads: %d", processing_stats$decorated_reads), bullet = TRUE)
  cli_log(sprintf("Blank reads: %d", processing_stats$blank_reads), bullet = TRUE)
  cli_log(sprintf("Unclassified reads: %d", processing_stats$unclassified_reads), bullet = TRUE)
  cli_log(sprintf("Total modifications detected: %d", processing_stats$total_modifications), bullet = TRUE)
  cli_log("NOTE: Results generated using polyA model (temporary)", "INFO", bullet = TRUE)

  cli_log("PolyT reads processing completed successfully", "SUCCESS")

  # Return combined results with tail_type
  # Note: Using 'nonadenosine_residues' name for consistency, but these represent non-T residues in polyT context
  return(list(
    read_classes = outputs$read_classes,
    nonadenosine_residues = outputs$nonadenosine_residues,  # Non-thymidine residues in polyT tails
    processing_stats = processing_stats,
    tail_type = "polyT"
  ))
}


################################################################################
# cDNA-SPECIFIC OUTPUT CREATION FUNCTION
################################################################################

#' Create Ninetails output tables for Dorado cDNA pipeline
#'
#' This function extends create_outputs_dorado to handle cDNA data with additional
#' tail_type information. It produces the same output structure as the DRS pipeline
#' but preserves the tail_type column for read orientation information.
#'
#' @param dorado_summary_dir Character string. Path to a directory containing Dorado
#' summary files (.txt, .tsv, or .csv) with per-read poly(A) tail information and tail_type.
#'
#' @param nonA_temp_dir Character string. Path to a directory containing non-adenosine
#' prediction RDS files, generated from temporary models.
#'
#' @param polya_chunks_dir Character string. Path to a directory containing poly(A) chunk
#' RDS files used for position inference of predictions.
#'
#' @param num_cores Integer. Number of cores to use for parallelized file loading
#' and processing. Must be a positive integer. Default is 1.
#'
#' @param qc Logical. Whether to apply quality control filtering of terminal
#' predictions (removing predictions near the ends of poly(A) tails). Default is TRUE.
#'
#' @return A named list with two data frames (identical to DRS output + tail_type):
#' \describe{
#' \item{read_classes}{Data frame with per-read classification results, including
#' columns for read name, contig, poly(A) length, QC tag, class, comments, and tail_type.}
#' \item{nonadenosine_residues}{Data frame with per-chunk predictions of
#' non-adenosine residues, including read name, contig, predicted base,
#' estimated position within the poly(A) tail, poly(A) length, QC tag, and tail_type.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' results <- create_outputs_dorado_cdna(
#' dorado_summary_dir = "data/dorado_summaries",
#' nonA_temp_dir = "data/nonA_predictions",
#' polya_chunks_dir = "data/polya_chunks",
#' num_cores = 4,
#' qc = TRUE
#' )
#'
#' # Access read classifications with tail_type
#' head(results$read_classes)
#'
#' # Access non-adenosine residues with tail_type
#' head(results$nonadenosine_residues)
#' }
create_outputs_dorado_cdna <- function(dorado_summary_dir,
                                       nonA_temp_dir,
                                       polya_chunks_dir,
                                       num_cores = 1,
                                       qc = TRUE) {

  # Variable binding for R CMD check
  read_id <- alignment_genome <- alignment_mapq <- poly_tail_length <- NULL
  poly_tail_start <- poly_tail_end <- chunkname <- prediction <- NULL
  centr_signal_pos <- signal_length <- est_nonA_pos <- class <- comments <- tail_type <- NULL

  # Assertions
  if (missing(dorado_summary_dir)) stop("Dorado summary directory is missing.", call. = FALSE)
  if (missing(nonA_temp_dir)) stop("Non-A predictions directory is missing.", call. = FALSE)
  if (missing(polya_chunks_dir)) stop("Poly(A) chunks directory is missing.", call. = FALSE)

  assertthat::assert_that(is.numeric(num_cores), num_cores > 0,
                          msg = "Number of cores must be a positive integer")
  assertthat::assert_that(is.logical(qc),
                          msg = "QC must be a logical value")
  assertthat::assert_that(dir.exists(dorado_summary_dir),
                          msg = "Dorado summary directory must exist")
  assertthat::assert_that(dir.exists(nonA_temp_dir),
                          msg = "Non-A predictions directory must exist")
  assertthat::assert_that(dir.exists(polya_chunks_dir),
                          msg = "Poly(A) chunks directory must exist")

  ################################################################################
  # LOAD DATA FILES
  ################################################################################

  # Load Dorado summary data with tail_type support
  dorado_summary_files <- list.files(dorado_summary_dir, pattern = "\\.(txt|tsv|csv)$", full.names = TRUE)

  if (length(dorado_summary_files) == 0) {
    stop("No summary files found in dorado_summary_dir", call. = FALSE)
  }

  # Load and combine all summary files
  dorado_summary_list <- parallel::mclapply(dorado_summary_files, function(file) {
    vroom::vroom(file, show_col_types = FALSE)
  }, mc.cores = num_cores)

  dorado_summary <- do.call(rbind, dorado_summary_list)

  # Verify required columns exist
  required_cols <- c("read_id", "alignment_genome", "alignment_mapq", "poly_tail_length",
                     "poly_tail_start", "poly_tail_end")
  missing_cols <- setdiff(required_cols, colnames(dorado_summary))

  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns in summary data: %s", paste(missing_cols, collapse = ", ")),
         call. = FALSE)
  }

  # Check if tail_type column exists (should be present in cDNA data)
  has_tail_type <- "tail_type" %in% colnames(dorado_summary)
  if (!has_tail_type) {
    warning("tail_type column not found in summary data - adding default 'unknown' values")
    dorado_summary$tail_type <- "unknown"
  }

  # Load non-A predictions
  nonA_files <- list.files(nonA_temp_dir, pattern = "\\.rds$", full.names = TRUE)
  chunks_files <- list.files(polya_chunks_dir, pattern = "\\.rds$", full.names = TRUE)

  ################################################################################
  # PROCESS PREDICTIONS (SAME LOGIC AS DRS BUT PRESERVE tail_type)
  ################################################################################

  moved_chunks_table <- data.frame()
  moved_blank_readnames <- character(0)

  if (length(nonA_files) > 0 && length(chunks_files) > 0) {
    # Load all prediction and chunk files
    all_predictions <- parallel::mclapply(nonA_files, readRDS, mc.cores = num_cores)
    all_chunks <- parallel::mclapply(chunks_files, readRDS, mc.cores = num_cores)

    # Combine predictions and chunks
    combined_predictions <- do.call(rbind, all_predictions)
    combined_chunks <- do.call(c, all_chunks)

    if (nrow(combined_predictions) > 0 && length(combined_chunks) > 0) {
      # Convert chunks to data frame
      chunks_df <- do.call(rbind, lapply(names(combined_chunks), function(read_id) {
        chunks <- combined_chunks[[read_id]]
        if (length(chunks) > 0) {
          data.frame(
            read_id = read_id,
            chunkname = names(chunks),
            centr_signal_pos = sapply(chunks, function(x) x$centr_signal_pos),
            stringsAsFactors = FALSE
          )
        }
      }))

      # Merge with predictions
      moved_chunks_table <- merge(combined_predictions, chunks_df, by = c("read_id", "chunkname"), all.x = TRUE)
      moved_chunks_table <- moved_chunks_table[!is.na(moved_chunks_table$centr_signal_pos), ]

      if (nrow(moved_chunks_table) > 0) {
        # Clean chunkname
        moved_chunks_table$chunkname <- gsub("chunk_\\d+\\*", "", moved_chunks_table$chunkname)

        # Create position list and merge with summary
        non_a_position_list <- moved_chunks_table[, c("read_id", "chunkname")]
        non_a_position_list <- merge(non_a_position_list,
                                     dorado_summary[, c("read_id", "poly_tail_length", "poly_tail_start",
                                                        "poly_tail_end", "tail_type")],
                                     by = "read_id")
        non_a_position_list$signal_length <- 0.2 * (non_a_position_list$poly_tail_end - non_a_position_list$poly_tail_start)

        # Final merge with tail_type preservation
        moved_chunks_table <- merge(moved_chunks_table, non_a_position_list, by = c("read_id", "chunkname"))

        # Calculate position
        moved_chunks_table$est_nonA_pos <- round(
          moved_chunks_table$poly_tail_length - ((moved_chunks_table$poly_tail_length * moved_chunks_table$centr_signal_pos) / moved_chunks_table$signal_length),
          2
        )

        # Merge with full summary data to get additional columns including tail_type
        moved_chunks_table <- merge(moved_chunks_table,
                                    dorado_summary[, c("read_id", "alignment_genome", "alignment_mapq", "tail_type")],
                                    by = "read_id")

        # Select final columns (including tail_type)
        moved_chunks_table <- moved_chunks_table[, c("read_id", "alignment_genome", "prediction", "est_nonA_pos",
                                                     "poly_tail_length", "signal_length", "alignment_mapq", "tail_type")]
      }
    }
  }

  # Get all read IDs
  all_read_ids <- unique(dorado_summary$read_id)

  # Quality control & sanity check - same as DRS
  if (qc == TRUE && nrow(moved_chunks_table) > 0) {
    terminal_mask <- (moved_chunks_table$est_nonA_pos < 2) |
      (moved_chunks_table$est_nonA_pos > moved_chunks_table$poly_tail_length - 2)

    moved_chunks_table_discarded_ids <- unique(moved_chunks_table$read_id[terminal_mask])
    moved_chunks_table <- moved_chunks_table[!terminal_mask, ]

    moved_blank_readnames <- unique(c(moved_blank_readnames, moved_chunks_table_discarded_ids))
  }

  decorated_read_ids <- unique(moved_chunks_table$read_id)
  blank_read_ids <- setdiff(all_read_ids, c(decorated_read_ids, moved_blank_readnames))

  # Vectorized read classification (same logic as DRS)
  dorado_summary$class <- "blank"  # Default
  dorado_summary$comments <- "MAU"  # Default

  # Update classifications
  dorado_summary$class[dorado_summary$read_id %in% decorated_read_ids] <- "decorated"
  dorado_summary$comments[dorado_summary$read_id %in% decorated_read_ids] <- "YAY"

  dorado_summary$class[dorado_summary$read_id %in% moved_blank_readnames] <- "blank"
  dorado_summary$comments[dorado_summary$read_id %in% moved_blank_readnames] <- "MPU"

  dorado_summary$class[dorado_summary$poly_tail_length < 10] <- "unclassified"
  dorado_summary$comments[dorado_summary$poly_tail_length < 10] <- "IRL"

  read_classes <- dorado_summary
  read_classes <- read_classes[!duplicated(read_classes$read_id), ]

  # Format read_classes (same as DRS + tail_type)
  names(read_classes)[names(read_classes) == "read_id"] <- "readname"
  names(read_classes)[names(read_classes) == "alignment_genome"] <- "contig"
  names(read_classes)[names(read_classes) == "poly_tail_length"] <- "polya_length"
  names(read_classes)[names(read_classes) == "alignment_mapq"] <- "qc_tag"

  # Select columns in same order as DRS but include tail_type
  read_classes <- read_classes[, c("readname", "contig", "polya_length", "qc_tag", "class", "comments", "tail_type")]

  # Format moved_chunks_table (same as DRS + tail_type)
  if (nrow(moved_chunks_table) > 0) {
    names(moved_chunks_table)[names(moved_chunks_table) == "read_id"] <- "readname"
    names(moved_chunks_table)[names(moved_chunks_table) == "alignment_genome"] <- "contig"
    names(moved_chunks_table)[names(moved_chunks_table) == "poly_tail_length"] <- "polya_length"
    names(moved_chunks_table)[names(moved_chunks_table) == "alignment_mapq"] <- "qc_tag"

    # Select columns in same order as DRS but include tail_type
    moved_chunks_table <- moved_chunks_table[, c("readname", "contig", "prediction", "est_nonA_pos", "polya_length", "qc_tag", "tail_type")]
  } else {
    moved_chunks_table <- data.frame(
      readname = character(), contig = character(), prediction = character(),
      est_nonA_pos = numeric(), polya_length = numeric(), qc_tag = character(), tail_type = character(),
      stringsAsFactors = FALSE
    )
  }

  # Final output (identical structure to DRS + tail_type)
  ninetails_output <- list(
    read_classes = read_classes,
    nonadenosine_residues = moved_chunks_table
  )

  # Summary statistics
  cat("cDNA Processing complete!\n")
  cat("Decorated reads:", sum(read_classes$class == "decorated", na.rm = TRUE), "\n")
  cat("Blank reads:", sum(read_classes$class == "blank", na.rm = TRUE), "\n")
  cat("Unclassified reads:", sum(read_classes$class == "unclassified", na.rm = TRUE), "\n")
  cat("Total nonadenosine residues:", nrow(moved_chunks_table), "\n")

  if (has_tail_type) {
    cat("PolyA reads:", sum(read_classes$tail_type == "polyA", na.rm = TRUE), "\n")
    cat("PolyT reads:", sum(read_classes$tail_type == "polyT", na.rm = TRUE), "\n")
    cat("Unidentified reads:", sum(read_classes$tail_type == "unidentified", na.rm = TRUE), "\n")
  }

  return(ninetails_output)
}



################################################################################
# RESULTS MERGING FUNCTION FOR cDNA (UPDATED FOR SIMPLIFIED APPROACH)
################################################################################

#' Merge polyA and polyT processing results for cDNA analysis
#'
#' This function combines the results from separate polyA and polyT processing
#' paths into a unified output in standard ninetails format. Updated to work
#' with the simplified approach using tibbles instead of files.
#'
#' @param polya_results List. Results from polyA processing (from process_polya_reads_cdna).
#' Can be NULL if no polyA reads were found.
#' @param polyt_results List. Results from polyT processing (from process_polyt_reads_cdna).
#' Can be NULL if no polyT reads were found.
#' @param unidentified_reads Data frame/tibble. Reads that could not be classified
#' as polyA or polyT (with tail_type = "unidentified").
#' @param save_dir Character. Directory where merged results will be saved.
#' @param prefix Character. Optional prefix for output file names.
#' @param cli_log Function for logging messages and progress.
#'
#' @return List containing merged cDNA results in standard ninetails format:
#'   \describe{
#'     \item{read_classes}{Data frame with all read classifications including tail_type}
#'     \item{nonadenosine_residues}{Data frame with all predicted modifications including tail_type}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' merged_results <- merge_cdna_results(
#'   polya_results = polya_output,
#'   polyt_results = polyt_output,
#'   unidentified_reads = unidentified_tibble,
#'   save_dir = "path/to/output/",
#'   prefix = "experiment1",
#'   cli_log = message
#' )
#' }
merge_cdna_results <- function(polya_results = NULL,
                               polyt_results = NULL,
                               unidentified_reads = NULL,
                               save_dir,
                               prefix = "",
                               cli_log = message) {

  # Input validation
  if (missing(save_dir)) {
    stop("Output directory is missing. Please provide a valid save_dir argument.",
         call. = FALSE)
  }

  assertthat::assert_that(is.character(save_dir),
                          msg = "save_dir must be a character string")

  cli_log("Starting cDNA results merging", "INFO", "Merging cDNA Results")

  # Check if we have any results to merge
  if (is.null(polya_results) && is.null(polyt_results)) {
    cli_log("WARNING: No polyA or polyT results provided for merging", "WARNING")
  }

  ################################################################################
  # MERGE READ CLASSIFICATIONS
  ################################################################################

  cli_log("Merging read classifications...", "INFO", bullet = TRUE)

  # Initialize empty data frame with expected structure
  merged_read_classes <- data.frame()

  # Add polyA read classes
  if (!is.null(polya_results) && !is.null(polya_results$read_classes)) {
    polya_classes <- polya_results$read_classes

    # Ensure tail_type column exists
    if (!"tail_type" %in% colnames(polya_classes)) {
      polya_classes$tail_type <- "polyA"
    }

    merged_read_classes <- rbind(merged_read_classes, polya_classes)
    cli_log(sprintf("Added %d polyA read classifications", nrow(polya_classes)), "INFO")
  } else {
    cli_log("No polyA read classifications to merge", "INFO")
  }

  # Add polyT read classes
  if (!is.null(polyt_results) && !is.null(polyt_results$read_classes)) {
    polyt_classes <- polyt_results$read_classes

    # Ensure tail_type column exists
    if (!"tail_type" %in% colnames(polyt_classes)) {
      polyt_classes$tail_type <- "polyT"
    }

    merged_read_classes <- rbind(merged_read_classes, polyt_classes)
    cli_log(sprintf("Added %d polyT read classifications", nrow(polyt_classes)), "INFO")
  } else {
    cli_log("No polyT read classifications to merge", "INFO")
  }

  ################################################################################
  # MERGE MODIFICATION PREDICTIONS
  ################################################################################

  cli_log("Merging modification predictions...", "INFO", bullet = TRUE)

  # Initialize empty data frame for modifications
  merged_modifications <- data.frame()

  # Add polyA modifications
  if (!is.null(polya_results) && !is.null(polya_results$nonadenosine_residues)) {
    polya_mods <- polya_results$nonadenosine_residues

    # Ensure tail_type column exists
    if (!"tail_type" %in% colnames(polya_mods)) {
      polya_mods$tail_type <- "polyA"
    }

    merged_modifications <- rbind(merged_modifications, polya_mods)
    cli_log(sprintf("Added %d polyA modifications", nrow(polya_mods)), "INFO")
  } else {
    cli_log("No polyA modifications to merge", "INFO")
  }

  # Add polyT modifications
  if (!is.null(polyt_results) && !is.null(polyt_results$nonadenosine_residues)) {
    polyt_mods <- polyt_results$nonadenosine_residues

    # Ensure tail_type column exists
    if (!"tail_type" %in% colnames(polyt_mods)) {
      polyt_mods$tail_type <- "polyT"
    }

    merged_modifications <- rbind(merged_modifications, polyt_mods)
    cli_log(sprintf("Added %d polyT modifications", nrow(polyt_mods)), "INFO")
  } else {
    cli_log("No polyT modifications to merge", "INFO")
  }

  ################################################################################
  # FINAL SUMMARY
  ################################################################################

  # Calculate totals for summary
  polya_count <- if (!is.null(polya_results))
    ifelse(!is.null(polya_results$processing_stats), polya_results$processing_stats$total_reads_processed, 0) else 0
  polyt_count <- if (!is.null(polyt_results))
    ifelse(!is.null(polyt_results$processing_stats), polyt_results$processing_stats$total_reads_processed, 0) else 0
  unidentified_count <- if (!is.null(unidentified_reads) && is.data.frame(unidentified_reads)) nrow(unidentified_reads) else 0
  total_modifications <- nrow(merged_modifications)
  total_reads <- nrow(merged_read_classes)
  grand_total <- polya_count + polyt_count + unidentified_count

  cli_log("Merging Summary", "INFO", "Merging Summary")
  cli_log(sprintf("Total reads in dataset: %d", grand_total), bullet = TRUE)
  cli_log(sprintf("- PolyA reads processed: %d (%.2f%%)",
                  polya_count, if (grand_total > 0) polya_count / grand_total * 100 else 0), bullet = TRUE)
  cli_log(sprintf("- PolyT reads processed: %d (%.2f%%)",
                  polyt_count, if (grand_total > 0) polyt_count / grand_total * 100 else 0), bullet = TRUE)
  cli_log(sprintf("- Unidentified reads: %d (%.2f%%)",
                  unidentified_count, if (grand_total > 0) unidentified_count / grand_total * 100 else 0), bullet = TRUE)

  cli_log(sprintf("Total modifications detected: %d", total_modifications), bullet = TRUE)

  if (polyt_count > 0) {
    cli_log("NOTE: PolyT modifications detected using temporary polyA model", "INFO", bullet = TRUE)
  }

  cli_log("cDNA results merging completed successfully", "SUCCESS")

  ################################################################################
  # RETURN STANDARD NINETAILS FORMAT
  ################################################################################

  # Return in standard ninetails format (same as DRS pipeline)
  return(list(
    read_classes = merged_read_classes,
    nonadenosine_residues = merged_modifications
  ))
}


################################################################################
# OUTPUT SAVING FUNCTION FOR cDNA PIPELINE
################################################################################

#' Save cDNA pipeline outputs in standard ninetails format
#'
#' This function takes the merged cDNA results and formats them as the standard
#' ninetails output format (same as create_outputs_dorado) with two tables:
#' read_classes and nonadenosine_residues. The key difference from DRS is the
#' additional tail_type column indicating polyA or polyT.
#'
#' @param outputs List. Merged results from merge_cdna_results().
#' @param save_dir Character string. Directory where outputs will be saved.
#' @param prefix Character string. Optional prefix for output file names.
#'
#' @return List containing standard ninetails output format:
#'   \describe{
#'     \item{read_classes}{Data frame with read classifications including tail_type}
#'     \item{nonadenosine_residues}{Data frame with modifications including tail_type}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' final_results <- save_cdna_outputs(
#'   outputs = merged_results,
#'   save_dir = "path/to/output/",
#'   prefix = "experiment1"
#' )
#' }
save_cdna_outputs <- function(outputs, save_dir, prefix = "") {

  # Input validation
  if (missing(outputs)) {
    stop("Outputs are missing. Please provide valid outputs argument.", call. = FALSE)
  }

  if (missing(save_dir)) {
    stop("Output directory is missing. Please provide a valid save_dir argument.", call. = FALSE)
  }

  assertthat::assert_that(is.list(outputs),
                          msg = "outputs must be a list from merge_cdna_results()")

  # Extract the two standard ninetails tables
  read_classes <- if (!is.null(outputs$read_classes)) {
    outputs$read_classes
  } else {
    data.frame()
  }

  nonadenosine_residues <- if (!is.null(outputs$nonadenosine_residues)) {
    outputs$nonadenosine_residues
  } else {
    data.frame()
  }

  # Return in standard ninetails format (same as create_outputs_dorado)
  return(list(
    read_classes = read_classes,
    nonadenosine_residues = nonadenosine_residues
  ))
}

















