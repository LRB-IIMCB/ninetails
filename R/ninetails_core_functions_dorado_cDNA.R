
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
# SEQUENCE EXTRACTION FUNCTION FROM BAM FILES
################################################################################

#' Extract basecalled sequences from BAM file for cDNA analysis
#'
#' This function extracts basecalled sequences from BAM files for reads that
#' are present in the corresponding Dorado summary file. It returns a data frame
#' with read IDs and their corresponding basecalled sequences, which can then
#' be used for orientation classification (polyA vs polyT) in the cDNA pipeline.
#'
#' @param bam_file Character string. Path to BAM file containing basecalled sequences.
#' @param summary_file Character string. Path to corresponding Dorado summary file
#' containing read IDs to extract.
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
#' sequences <- extract_sequences_from_bam(
#'   bam_file = "aligned_reads.bam",
#'   summary_file = "dorado_summary.txt",
#'   cli_log = message
#' )
#' }
extract_sequences_from_bam <- function(bam_file, summary_file, cli_log = message) {

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
  #ref_starts <- integer(initial_count)
  #ref_ends <- integer(initial_count)
  mapqs <- integer(initial_count)
  poly_tail_lengths <- numeric(initial_count)
  #anchor_positions <- integer(initial_count)
  #poly_tail_starts <- integer(initial_count)
  #polya_ends <- integer(initial_count)
  #poly_tail2_starts <- integer(initial_count)
  #poly_tail2_ends <- integer(initial_count)
  pod5_files <- character(initial_count)  # New vector for pod5 files
  seqs <- character(initial_count) # basecalled sequences

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

    # Initialize pa tag values
    #anchor_pos <- -1L
    #poly_tail_start <- -1L
    #polya_end <- -1L
    #poly_tail2_start <- -1L
    #poly_tail2_end <- -1L

    # Add pa tag information if it exists and is valid
    # if (!is.null(bam_data$tag$pa[[i]]) && length(bam_data$tag$pa[[i]]) == 5) {
    #   pa_values <- bam_data$tag$pa[[i]]
    #   if (!any(is.na(pa_values))) {
    #     #anchor_pos <- pa_values[1]
    #     poly_tail_start <- pa_values[2]
    #     polya_end <- pa_values[3]
    #     poly_tail2_start <- pa_values[4]
    #     poly_tail2_end <- pa_values[5]
    #   }
    # }

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
    #ref_starts[valid_count] <- bam_data$pos[i]
    #ref_ends[valid_count] <- bam_data$pos[i] + nchar(bam_data$cigar[i])
    mapqs[valid_count] <- bam_data$mapq[i]
    #poly_tail_lengths[valid_count] <- poly_tail_length
    #anchor_positions[valid_count] <- anchor_pos
    #poly_tail_starts[valid_count] <- poly_tail_start
    #polya_ends[valid_count] <- polya_end
    #poly_tail2_starts[valid_count] <- poly_tail2_start
    #poly_tail2_ends[valid_count] <- poly_tail2_end
    pod5_files[valid_count] <- pod5_file
    seqs[valid_count] <- bam_data$seq[i] # basecalled sequences

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
    #ref_starts <- ref_starts[1:valid_count]
    #ref_ends <- ref_ends[1:valid_count]
    mapqs <- mapqs[1:valid_count]
    #poly_tail_lengths <- poly_tail_lengths[1:valid_count]
    #anchor_positions <- anchor_positions[1:valid_count]
    #poly_tail_starts <- poly_tail_starts[1:valid_count]
    #polya_ends <- polya_ends[1:valid_count]
    #poly_tail2_starts <- poly_tail2_starts[1:valid_count]
    #poly_tail2_ends <- poly_tail2_ends[1:valid_count]
    pod5_files <- pod5_files[1:valid_count]
    seqs <- seqs[1:valid_count] # basecalled sequences

    # Check for duplicates
    duplicate_reads <- duplicated(read_ids)
    if (any(duplicate_reads)) {
      cli_log(sprintf("WARNING: Found %d duplicate read names, keeping first occurrence only",
                      sum(duplicate_reads)), "WARNING")
      # Keep only unique entries
      keep_idx <- !duplicate_reads
      read_ids <- read_ids[keep_idx]
      references <- references[keep_idx]
      #ref_starts <- ref_starts[keep_idx]
      #ref_ends <- ref_ends[keep_idx]
      #mapqs <- mapqs[keep_idx]
      #poly_tail_lengths <- poly_tail_lengths[keep_idx]
      #anchor_positions <- anchor_positions[keep_idx]
      #poly_tail_starts <- poly_tail_starts[keep_idx]
      #polya_ends <- polya_ends[keep_idx]
      #poly_tail2_starts <- poly_tail2_starts[keep_idx]
      #poly_tail2_ends <- poly_tail2_ends[keep_idx]
      pod5_files <- pod5_files[keep_idx]
      seqs <- seqs[keep_idx]
    }

    # Create the final data frame
    polya_df <- tibble::tibble(
      read_id = read_ids,
      pod5_file = pod5_files,  # Added pod5 file information
      reference = references,
      #ref_start = ref_starts,
      #ref_end = ref_ends,
      #mapq = mapqs,
      #poly_tail_length = poly_tail_lengths,
      #anchor_pos = anchor_positions,
      #poly_tail_start = poly_tail_starts,
      #polya_end = polya_ends,
      #poly_tail2_start = poly_tail2_starts,
      #poly_tail2_end = poly_tail2_ends,
      sequence = seqs # bug fix: column name adjustment
    )

    cli_log("Processing complete", "SUCCESS")
    return(polya_df)
  } else {
    cli_log("No reads with basecalled sequences found after filtering", "WARNING")
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
#' extracts basecalled sequences from BAM files for orientation classification.
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

  ################################################################################
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

  ################################################################################
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

  ################################################################################
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
    sequence_data <- ninetails::extract_sequences_from_bam(
      bam_file = current_bam,
      summary_file = current_summary,
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



  ################################################################################
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

  ################################################################################
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
# READ ORIENTATION CLASSIFICATION FUNCTION
################################################################################

#' Classify read orientations based on primer alignment
#'
#' This function classifies cDNA reads as polyA, polyT, or unidentified based on
#' local alignment of specific primers to the first portion of each read sequence.
#' The classification determines which processing path (polyA or polyT) each read
#' should follow in the downstream analysis.
#'
#' @param sequence_files Character vector. Paths to sequence files (TSV format)
#' containing read_id and sequence columns.
#' @param save_dir Character string. Directory where classified sequence files
#' will be saved.
#' @param prefix Character string. Optional prefix for output file names.
#' @param num_cores Integer. Number of CPU cores to use for parallel processing.
#' @param cli_log Function for logging messages and progress.
#'
#' @return List containing paths to classified sequence files:
#'   \describe{
#'     \item{polya_files}{Paths to files containing polyA-classified sequences}
#'     \item{polyt_files}{Paths to files containing polyT-classified sequences}
#'     \item{unidentified_files}{Paths to files containing unidentified sequences}
#'     \item{classification_stats}{Summary statistics of classification results}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' classification_results <- classify_read_orientations(
#'   sequence_files = c("sequences_part1.tsv", "sequences_part2.tsv"),
#'   save_dir = "path/to/output/",
#'   prefix = "experiment1",
#'   num_cores = 4,
#'   cli_log = message
#' )
#' }
classify_read_orientations <- function(sequence_files,
                                       save_dir,
                                       prefix = "",
                                       num_cores = 1,
                                       cli_log = message) {

  # Check for required Bioconductor packages
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Package 'Biostrings' is required for classify_read_orientations(). Please install it from Bioconductor.",
         call. = FALSE)
  }
  if (!requireNamespace("pwalign", quietly = TRUE)) {
    stop("Package 'pwalign' is required for classify_read_orientations(). Please install it from Bioconductor.",
         call. = FALSE)
  }

  # Input validation
  if (missing(sequence_files)) {
    stop("Sequence files are missing. Please provide valid sequence_files argument.",
         call. = FALSE)
  }

  if (missing(save_dir)) {
    stop("Output directory is missing. Please provide a valid save_dir argument.",
         call. = FALSE)
  }

  assertthat::assert_that(is.character(sequence_files),
                          msg = "sequence_files must be a character vector of file paths")
  assertthat::assert_that(all(file.exists(sequence_files)),
                          msg = "All sequence files must exist")
  assertthat::assert_that(is.numeric(num_cores), num_cores > 0,
                          msg = "Number of cores must be a positive numeric value")

  cli_log("Starting read orientation classification", "INFO", "Sequence Classification")

  # Create output directories
  classification_dir <- file.path(save_dir, "classification_dir")
  if (!dir.exists(classification_dir)) {
    dir.create(classification_dir, recursive = TRUE)
  }

  polya_dir <- file.path(classification_dir, "polya_sequences")
  polyt_dir <- file.path(classification_dir, "polyt_sequences")
  unidentified_dir <- file.path(classification_dir, "unidentified_sequences")

  for (dir in c(polya_dir, polyt_dir, unidentified_dir)) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
  }

  # Define primer sequences and alignment parameters
  ################################################################################
  front_primer <- Biostrings::DNAString("TTTCTGTTGGTGCTGATATTGCTTT")
  rear_primer <- Biostrings::DNAString("ACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAGTTTT")

  # Alignment parameters
  match <- 1
  mismatch <- -1
  type <- 'local'
  gapOpening <- 0
  gapExtension <- 1

  # Create substitution matrix
  submat <- pwalign::nucleotideSubstitutionMatrix(match = match,
                                                  mismatch = mismatch,
                                                  baseOnly = TRUE)

  # Classification parameters
  threshold <- 0.6
  search_window <- 140

  cli_log(sprintf("Using alignment parameters: threshold=%.2f, search_window=%d",
                  threshold, search_window), "INFO")

  # Initialize output file vectors
  polya_files <- character(length(sequence_files))
  polyt_files <- character(length(sequence_files))
  unidentified_files <- character(length(sequence_files))

  # Initialize statistics
  total_polya <- 0
  total_polyt <- 0
  total_unidentified <- 0
  total_processed <- 0

  # Process each sequence file
  ################################################################################
  for (i in seq_along(sequence_files)) {
    current_file <- sequence_files[i]
    file_basename <- tools::file_path_sans_ext(basename(current_file))

    cli_log(sprintf("Processing sequence file %d/%d: %s",
                    i, length(sequence_files), basename(current_file)), "INFO")

    # Read sequence data
    tryCatch({
      sequence_data <- vroom::vroom(current_file, show_col_types = FALSE)

      # Validate required columns
      if (!all(c("read_id", "sequence") %in% colnames(sequence_data))) {
        cli_log(sprintf("Skipping file %s: missing required columns (read_id, sequence)",
                        basename(current_file)), "WARNING")
        next
      }

      if (nrow(sequence_data) == 0) {
        cli_log(sprintf("Skipping empty file: %s", basename(current_file)), "WARNING")
        next
      }

      cli_log(sprintf("Loaded %d sequences for classification", nrow(sequence_data)), "INFO")

      # Apply classification function to each sequence
      ##############################################################################

      # Classification function
      classify_single_read <- function(sequence_str) {
        # Handle empty or very short sequences
        if (is.na(sequence_str) || nchar(sequence_str) < 25) {
          return("unidentified")
        }

        # Extract search window from beginning of sequence
        search_seq <- substr(sequence_str, 1, min(search_window, nchar(sequence_str)))
        search_dna <- Biostrings::DNAString(search_seq)

        tryCatch({
          # Align both primers to the search window
          align_front <- pwalign::pairwiseAlignment(pattern = front_primer,
                                                    subject = search_dna,
                                                    substitutionMatrix = submat,
                                                    type = type,
                                                    scoreOnly = FALSE,
                                                    gapOpening = gapOpening,
                                                    gapExtension = gapExtension)

          align_rear <- pwalign::pairwiseAlignment(pattern = rear_primer,
                                                   subject = search_dna,
                                                   substitutionMatrix = submat,
                                                   type = type,
                                                   scoreOnly = FALSE,
                                                   gapOpening = gapOpening,
                                                   gapExtension = gapExtension)

          # Calculate normalized alignment scores
          score_front <- align_front@score / length(front_primer)
          score_rear <- align_rear@score / length(rear_primer)

          # Apply classification logic
          if (score_front > score_rear && score_front > threshold) {
            return("polyA")
          } else if (score_front < score_rear && score_rear > threshold) {
            return("polyT")
          } else {
            return("unidentified")
          }

        }, error = function(e) {
          return("unidentified")
        })
      }

      # Apply classification to all sequences
      cli_log("Applying primer alignment classification...", "INFO")

      # Use sequential processing for reliability
      classifications <- sapply(sequence_data$sequence, classify_single_read)

      # Add classification results to data frame
      sequence_data$tail_type <- classifications

      # Count results for this file
      file_polya <- sum(classifications == "polyA")
      file_polyt <- sum(classifications == "polyT")
      file_unidentified <- sum(classifications == "unidentified")

      cli_log(sprintf("Classification complete: %d polyA, %d polyT, %d unidentified",
                      file_polya, file_polyt, file_unidentified), "INFO")

      # Update totals
      total_polya <- total_polya + file_polya
      total_polyt <- total_polyt + file_polyt
      total_unidentified <- total_unidentified + file_unidentified
      total_processed <- total_processed + nrow(sequence_data)

      # Split and save classified sequences
      ##############################################################################

      # Save polyA sequences
      if (file_polya > 0) {
        polya_data <- sequence_data[sequence_data$tail_type == "polyA", ]
        polya_file <- file.path(polya_dir,
                                sprintf("%spolya_%s.tsv",
                                        ifelse(nchar(prefix) > 0, paste0(prefix, "_"), ""),
                                        file_basename))
        vroom::vroom_write(polya_data, polya_file, delim = "\t")
        polya_files[i] <- polya_file
        cli_log(sprintf("Saved %d polyA sequences to %s", file_polya, basename(polya_file)), "INFO")
      }

      # Save polyT sequences
      if (file_polyt > 0) {
        polyt_data <- sequence_data[sequence_data$tail_type == "polyT", ]
        polyt_file <- file.path(polyt_dir,
                                sprintf("%spolyt_%s.tsv",
                                        ifelse(nchar(prefix) > 0, paste0(prefix, "_"), ""),
                                        file_basename))
        vroom::vroom_write(polyt_data, polyt_file, delim = "\t")
        polyt_files[i] <- polyt_file
        cli_log(sprintf("Saved %d polyT sequences to %s", file_polyt, basename(polyt_file)), "INFO")
      }

      # Save unidentified sequences
      if (file_unidentified > 0) {
        unidentified_data <- sequence_data[sequence_data$tail_type == "unidentified", ]
        unidentified_file <- file.path(unidentified_dir,
                                       sprintf("%sunidentified_%s.tsv",
                                               ifelse(nchar(prefix) > 0, paste0(prefix, "_"), ""),
                                               file_basename))
        vroom::vroom_write(unidentified_data, unidentified_file, delim = "\t")
        unidentified_files[i] <- unidentified_file
        cli_log(sprintf("Saved %d unidentified sequences to %s",
                        file_unidentified, basename(unidentified_file)), "INFO")
      }

    }, error = function(e) {
      cli_log(sprintf("Error processing file %s: %s", basename(current_file), e$message), "ERROR")
    })

    # Trigger garbage collection
    gc()
  }

  # Filter out empty file paths and compile results
  ################################################################################
  polya_files <- polya_files[file.exists(polya_files)]
  polyt_files <- polyt_files[file.exists(polyt_files)]
  unidentified_files <- unidentified_files[file.exists(unidentified_files)]

  # Create classification statistics
  classification_stats <- list(
    total_reads = total_processed,
    polya_reads = total_polya,
    polyt_reads = total_polyt,
    unidentified_reads = total_unidentified,
    polya_percentage = round(total_polya / total_processed * 100, 2),
    polyt_percentage = round(total_polyt / total_processed * 100, 2),
    unidentified_percentage = round(total_unidentified / total_processed * 100, 2)
  )

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

  # Return results
  return(list(
    polya_files = polya_files,
    polyt_files = polyt_files,
    unidentified_files = unidentified_files,
    classification_stats = classification_stats
  ))
}




################################################################################
# POLYA READS PROCESSING FUNCTION FOR cDNA
################################################################################

#' Process polyA reads using standard ninetails pipeline
#'
#' This function processes reads that have been classified as polyA-containing
#' through the cDNA classification pipeline. It applies the standard ninetails
#' analysis pipeline including feature extraction, signal segmentation,
#' GAF creation, and neural network prediction to identify non-adenosine
#' residues within polyA tails.
#'
#' @param polya_files Character vector. Paths to polyA-classified sequence files.
#' @param polya_signal_files Character vector. Paths to signal files containing
#' poly(A) tail signals extracted from POD5 files.
#' @param num_cores Integer. Number of CPU cores to use for parallel processing.
#' @param qc Logical. Whether to apply quality control filtering.
#' @param save_dir Character. Directory where processing results will be saved.
#' @param prefix Character. Optional prefix for output file names.
#' @param cli_log Function for logging messages and progress.
#'
#' @return List containing polyA processing results:
#'   \describe{
#'     \item{read_classes}{Data frame with read classification results}
#'     \item{nonadenosine_residues}{Data frame with predicted modifications}
#'     \item{processing_stats}{Summary statistics for polyA processing}
#'     \item{tail_type}{Character indicating this is "polyA" data}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' polya_results <- process_polya_reads_cdna(
#'   polya_files = c("polya_sequences_part1.tsv"),
#'   polya_signal_files = c("polya_signal_part1.rds"),
#'   num_cores = 4,
#'   qc = TRUE,
#'   save_dir = "path/to/output/",
#'   prefix = "experiment1",
#'   cli_log = message
#' )
#' }
process_polya_reads_cdna <- function(polya_files,
                                     polya_signal_files,
                                     num_cores = 1,
                                     qc = TRUE,
                                     save_dir,
                                     prefix = "",
                                     cli_log = message) {

  # Input validation
  if (missing(polya_files)) {
    stop("PolyA files are missing. Please provide valid polya_files argument.",
         call. = FALSE)
  }

  if (missing(polya_signal_files)) {
    stop("PolyA signal files are missing. Please provide valid polya_signal_files argument.",
         call. = FALSE)
  }

  if (missing(save_dir)) {
    stop("Output directory is missing. Please provide a valid save_dir argument.",
         call. = FALSE)
  }

  assertthat::assert_that(is.character(polya_files),
                          msg = "polya_files must be a character vector of file paths")
  assertthat::assert_that(all(file.exists(polya_files)),
                          msg = "All polyA files must exist")
  assertthat::assert_that(is.character(polya_signal_files),
                          msg = "polya_signal_files must be a character vector of file paths")
  assertthat::assert_that(all(file.exists(polya_signal_files)),
                          msg = "All polyA signal files must exist")
  assertthat::assert_that(is.numeric(num_cores), num_cores > 0,
                          msg = "Number of cores must be a positive numeric value")

  cli_log("Starting polyA reads processing", "INFO", "PolyA Processing")

  # Create output directories for polyA processing
  polya_processing_dir <- file.path(save_dir, "polya_processing_dir")
  polya_temp_dir <- file.path(polya_processing_dir, "polya_temp_dir")
  polya_chunks_dir <- file.path(polya_processing_dir, "polya_chunks_dir")

  for (dir in c(polya_processing_dir, polya_temp_dir, polya_chunks_dir)) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
  }

  ################################################################################
  # FILTER SIGNALS FOR POLYA READS
  ################################################################################

  cli_log("Filtering signals for polyA reads...", "INFO", bullet = TRUE)

  # Collect all polyA read IDs from classified sequences
  polya_read_ids <- character(0)

  for (polya_file in polya_files) {
    tryCatch({
      polya_data <- vroom::vroom(polya_file, show_col_types = FALSE)
      if ("read_id" %in% colnames(polya_data)) {
        polya_read_ids <- c(polya_read_ids, polya_data$read_id)
      }
    }, error = function(e) {
      cli_log(sprintf("Error reading polyA file %s: %s", basename(polya_file), e$message), "WARNING")
    })
  }

  polya_read_ids <- unique(polya_read_ids)
  cli_log(sprintf("Found %d unique polyA read IDs for processing", length(polya_read_ids)), "INFO")

  # Filter signal files to keep only polyA reads
  filtered_signal_files <- character(length(polya_signal_files))

  for (i in seq_along(polya_signal_files)) {
    signal_file <- polya_signal_files[i]
    cli_log(sprintf("Filtering signal file %d/%d: %s",
                    i, length(polya_signal_files), basename(signal_file)), "INFO")

    tryCatch({
      # Load signal data
      signal_list <- readRDS(signal_file)

      # Filter to keep only polyA reads
      polya_signals <- signal_list[names(signal_list) %in% polya_read_ids]

      cli_log(sprintf("Filtered signals: %d polyA reads from %d total reads",
                      length(polya_signals), length(signal_list)), "INFO")

      # Save filtered signals
      if (length(polya_signals) > 0) {
        output_file <- file.path(polya_processing_dir,
                                 sprintf("%spolya_filtered_signals_part%d.rds",
                                         ifelse(nchar(prefix) > 0, paste0(prefix, "_"), ""),
                                         i))
        saveRDS(polya_signals, output_file)
        filtered_signal_files[i] <- output_file
        cli_log(sprintf("Saved filtered signals to %s", basename(output_file)), "INFO")
      } else {
        cli_log(sprintf("No polyA signals found in part %d", i), "WARNING")
      }

    }, error = function(e) {
      cli_log(sprintf("Error filtering signal file %s: %s", basename(signal_file), e$message), "ERROR")
    })
  }

  # Remove empty entries
  filtered_signal_files <- filtered_signal_files[file.exists(filtered_signal_files)]

  if (length(filtered_signal_files) == 0) {
    stop("No polyA signals were found for processing")
  }

  ################################################################################
  # STANDARD NINETAILS PROCESSING PIPELINE
  ################################################################################

  # Process each filtered signal file through the ninetails pipeline
  all_outputs <- list()

  for (i in seq_along(filtered_signal_files)) {
    signal_file <- filtered_signal_files[i]

    cli_log(sprintf("Processing polyA signals %d/%d through ninetails pipeline",
                    i, length(filtered_signal_files)), "INFO")

    tryCatch({
      # Load signal data
      signal_list <- readRDS(signal_file)

      if (length(signal_list) == 0) {
        cli_log(sprintf("Skipping empty signal file: %s", basename(signal_file)), "WARNING")
        next
      }

      #####################################################
      # CREATE TAIL FEATURES LIST
      #####################################################
      cli_log("Creating tail features list...", "INFO", bullet = TRUE)

      tail_feature_list <- tryCatch({
        invisible(utils::capture.output(
          result <- ninetails::create_tail_features_list_dorado(
            signal_list = signal_list,
            num_cores = num_cores
          )
        ))
        result
      }, error = function(e) {
        cli_log(sprintf("Error in feature extraction: %s", e$message), "ERROR")
        stop(e$message, call. = FALSE)
      })

      cli_log(sprintf("Created features for %d reads", length(tail_feature_list)), "INFO")

      #####################################################
      # CREATE TAIL CHUNK LIST
      #####################################################
      cli_log("Creating tail segmentation data...", "INFO", bullet = TRUE)

      tail_chunk_list <- tryCatch({
        invisible(utils::capture.output(
          result <- ninetails::create_tail_chunk_list_dorado(
            tail_feature_list = tail_feature_list,
            num_cores = num_cores
          )
        ))
        result
      }, error = function(e) {
        cli_log(sprintf("Error in chunk creation: %s", e$message), "ERROR")
        stop(e$message, call. = FALSE)
      })

      cli_log(sprintf("Created chunks for %d reads", length(tail_chunk_list)), "INFO")

      #####################################################
      # CREATE GAF LIST
      #####################################################
      cli_log("Computing gramian angular fields...", "INFO", bullet = TRUE)

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

      cli_log(sprintf("Created GAFs for %d chunks", length(gaf_list)), "INFO")

      #####################################################
      # PREDICT CLASSES
      #####################################################
      cli_log("Running predictions with polyA model...", "INFO", bullet = TRUE)

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

      cli_log(sprintf("Completed predictions for %d GAFs", length(predicted_list$chunkname)), "INFO")

      #####################################################
      # CREATE OUTPUTS FOR THIS PART
      #####################################################
      cli_log("Creating outputs...", "INFO", bullet = TRUE)

      # Create a temporary summary for this subset of polyA reads
      # We need to reconstruct summary information for the output creation function
      temp_summary <- data.frame(
        readname = names(signal_list),
        polya_length = sapply(signal_list, length),
        qc_tag = "PASS"  # Assume all classified polyA reads pass QC
      )

      part_outputs <- tryCatch({
        invisible(utils::capture.output(
          result <- ninetails::create_outputs(
            tail_feature_list = tail_feature_list,
            tail_chunk_list = tail_chunk_list,
            nanopolish = temp_summary,  # Use temp summary
            predicted_list = predicted_list,
            num_cores = num_cores,
            pass_only = FALSE,  # We've already filtered
            qc = qc
          )
        ))
        result
      }, error = function(e) {
        cli_log(sprintf("Error in output creation: %s", e$message), "ERROR")
        stop(e$message, call. = FALSE)
      })

      cli_log("Output creation completed", "SUCCESS")

      # Store this part's outputs
      all_outputs[[i]] <- part_outputs

      # Save intermediate results
      output_basename <- sprintf("%spolya_results_part%d.rds",
                                 ifelse(nchar(prefix) > 0, paste0(prefix, "_"), ""),
                                 i)
      output_file <- file.path(polya_temp_dir, output_basename)
      saveRDS(part_outputs, output_file)

      cli_log(sprintf("Saved intermediate results to %s", basename(output_file)), "INFO")

    }, error = function(e) {
      cli_log(sprintf("Error processing signal file %d: %s", i, e$message), "ERROR")
    })

    # Explicit garbage collection
    gc()
  }

  ################################################################################
  # MERGE ALL POLYΑ RESULTS
  ################################################################################

  cli_log("Merging polyA processing results...", "INFO", "Merging PolyA Results")

  if (length(all_outputs) == 0) {
    stop("No polyA results were generated")
  }

  # Combine read classes from all parts
  all_read_classes <- do.call(rbind, lapply(all_outputs, function(x) {
    if (!is.null(x$read_classes)) {
      return(x$read_classes)
    } else {
      return(data.frame())
    }
  }))

  # Combine nonadenosine residues from all parts
  all_nonadenosine <- do.call(rbind, lapply(all_outputs, function(x) {
    if (!is.null(x$nonadenosine_residues)) {
      return(x$nonadenosine_residues)
    } else {
      return(data.frame())
    }
  }))

  # Add tail type information to results
  if (nrow(all_read_classes) > 0) {
    all_read_classes$tail_type <- "polyA"
  }

  if (nrow(all_nonadenosine) > 0) {
    all_nonadenosine$tail_type <- "polyA"
  }

  ################################################################################
  # COMPILE FINAL RESULTS
  ################################################################################

  # Calculate processing statistics
  processing_stats <- list(
    total_reads_processed = nrow(all_read_classes),
    decorated_reads = sum(all_read_classes$class == "decorated", na.rm = TRUE),
    blank_reads = sum(all_read_classes$class == "blank", na.rm = TRUE),
    unclassified_reads = sum(all_read_classes$class == "unclassified", na.rm = TRUE),
    total_modifications = nrow(all_nonadenosine),
    tail_type = "polyA"
  )

  cli_log("PolyA Processing Summary", "INFO", "PolyA Summary")
  cli_log(sprintf("Total polyA reads processed: %d", processing_stats$total_reads_processed), bullet = TRUE)
  cli_log(sprintf("Decorated reads: %d", processing_stats$decorated_reads), bullet = TRUE)
  cli_log(sprintf("Blank reads: %d", processing_stats$blank_reads), bullet = TRUE)
  cli_log(sprintf("Unclassified reads: %d", processing_stats$unclassified_reads), bullet = TRUE)
  cli_log(sprintf("Total modifications detected: %d", processing_stats$total_modifications), bullet = TRUE)

  cli_log("PolyA reads processing completed successfully", "SUCCESS")

  # Return combined results
  return(list(
    read_classes = all_read_classes,
    nonadenosine_residues = all_nonadenosine,
    processing_stats = processing_stats,
    tail_type = "polyA"
  ))
}



################################################################################
# POLYT READS PROCESSING FUNCTION FOR cDNA
################################################################################

#' Process polyT reads using ninetails pipeline
#'
#' This function processes reads that have been classified as polyT-containing
#' through the cDNA classification pipeline. It applies the ninetails analysis
#' pipeline including feature extraction, signal segmentation, GAF creation,
#' and neural network prediction to identify non-thymidine residues within
#' polyT tails. Currently uses the same model as polyA processing until a
#' polyT-specific model is trained.
#'
#' @param polyt_files Character vector. Paths to polyT-classified sequence files.
#' @param polyt_signal_files Character vector. Paths to signal files containing
#' poly(A) tail signals extracted from POD5 files (will be filtered for polyT reads).
#' @param num_cores Integer. Number of CPU cores to use for parallel processing.
#' @param qc Logical. Whether to apply quality control filtering.
#' @param save_dir Character. Directory where processing results will be saved.
#' @param prefix Character. Optional prefix for output file names.
#' @param cli_log Function for logging messages and progress.
#'
#' @return List containing polyT processing results:
#'   \describe{
#'     \item{read_classes}{Data frame with read classification results}
#'     \item{nonadenosine_residues}{Data frame with predicted modifications (renamed for polyT context)}
#'     \item{processing_stats}{Summary statistics for polyT processing}
#'     \item{tail_type}{Character indicating this is "polyT" data}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' polyt_results <- process_polyt_reads_cdna(
#'   polyt_files = c("polyt_sequences_part1.tsv"),
#'   polyt_signal_files = c("polya_signal_part1.rds"),
#'   num_cores = 4,
#'   qc = TRUE,
#'   save_dir = "path/to/output/",
#'   prefix = "experiment1",
#'   cli_log = message
#' )
#' }
process_polyt_reads_cdna <- function(polyt_files,
                                     polyt_signal_files,
                                     num_cores = 1,
                                     qc = TRUE,
                                     save_dir,
                                     prefix = "",
                                     cli_log = message) {

  # Input validation
  if (missing(polyt_files)) {
    stop("PolyT files are missing. Please provide valid polyt_files argument.",
         call. = FALSE)
  }

  if (missing(polyt_signal_files)) {
    stop("PolyT signal files are missing. Please provide valid polyt_signal_files argument.",
         call. = FALSE)
  }

  if (missing(save_dir)) {
    stop("Output directory is missing. Please provide a valid save_dir argument.",
         call. = FALSE)
  }

  assertthat::assert_that(is.character(polyt_files),
                          msg = "polyt_files must be a character vector of file paths")
  assertthat::assert_that(all(file.exists(polyt_files)),
                          msg = "All polyT files must exist")
  assertthat::assert_that(is.character(polyt_signal_files),
                          msg = "polyt_signal_files must be a character vector of file paths")
  assertthat::assert_that(all(file.exists(polyt_signal_files)),
                          msg = "All polyT signal files must exist")
  assertthat::assert_that(is.numeric(num_cores), num_cores > 0,
                          msg = "Number of cores must be a positive numeric value")

  cli_log("Starting polyT reads processing", "INFO", "PolyT Processing")

  # NOTE: Currently using same model as polyA until polyT-specific model is trained
  cli_log("NOTE: Currently using polyA model for polyT reads (temporary)", "INFO", bullet = TRUE)

  # Create output directories for polyT processing
  polyt_processing_dir <- file.path(save_dir, "polyt_processing_dir")
  polyt_temp_dir <- file.path(polyt_processing_dir, "polyt_temp_dir")
  polyt_chunks_dir <- file.path(polyt_processing_dir, "polyt_chunks_dir")

  for (dir in c(polyt_processing_dir, polyt_temp_dir, polyt_chunks_dir)) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
  }

  ################################################################################
  # FILTER SIGNALS FOR POLYT READS
  ################################################################################

  cli_log("Filtering signals for polyT reads...", "INFO", bullet = TRUE)

  # Collect all polyT read IDs from classified sequences
  polyt_read_ids <- character(0)

  for (polyt_file in polyt_files) {
    tryCatch({
      polyt_data <- vroom::vroom(polyt_file, show_col_types = FALSE)
      if ("read_id" %in% colnames(polyt_data)) {
        polyt_read_ids <- c(polyt_read_ids, polyt_data$read_id)
      }
    }, error = function(e) {
      cli_log(sprintf("Error reading polyT file %s: %s", basename(polyt_file), e$message), "WARNING")
    })
  }

  polyt_read_ids <- unique(polyt_read_ids)
  cli_log(sprintf("Found %d unique polyT read IDs for processing", length(polyt_read_ids)), "INFO")

  # Filter signal files to keep only polyT reads
  filtered_signal_files <- character(length(polyt_signal_files))

  for (i in seq_along(polyt_signal_files)) {
    signal_file <- polyt_signal_files[i]
    cli_log(sprintf("Filtering signal file %d/%d: %s",
                    i, length(polyt_signal_files), basename(signal_file)), "INFO")

    tryCatch({
      # Load signal data
      signal_list <- readRDS(signal_file)

      # Filter to keep only polyT reads
      polyt_signals <- signal_list[names(signal_list) %in% polyt_read_ids]

      cli_log(sprintf("Filtered signals: %d polyT reads from %d total reads",
                      length(polyt_signals), length(signal_list)), "INFO")

      # Save filtered signals
      if (length(polyt_signals) > 0) {
        output_file <- file.path(polyt_processing_dir,
                                 sprintf("%spolyt_filtered_signals_part%d.rds",
                                         ifelse(nchar(prefix) > 0, paste0(prefix, "_"), ""),
                                         i))
        saveRDS(polyt_signals, output_file)
        filtered_signal_files[i] <- output_file
        cli_log(sprintf("Saved filtered signals to %s", basename(output_file)), "INFO")
      } else {
        cli_log(sprintf("No polyT signals found in part %d", i), "WARNING")
      }

    }, error = function(e) {
      cli_log(sprintf("Error filtering signal file %s: %s", basename(signal_file), e$message), "ERROR")
    })
  }

  # Remove empty entries
  filtered_signal_files <- filtered_signal_files[file.exists(filtered_signal_files)]

  if (length(filtered_signal_files) == 0) {
    stop("No polyT signals were found for processing")
  }

  ################################################################################
  # NINETAILS PROCESSING PIPELINE FOR POLYT
  ################################################################################

  # Process each filtered signal file through the ninetails pipeline
  all_outputs <- list()

  for (i in seq_along(filtered_signal_files)) {
    signal_file <- filtered_signal_files[i]

    cli_log(sprintf("Processing polyT signals %d/%d through ninetails pipeline",
                    i, length(filtered_signal_files)), "INFO")

    tryCatch({
      # Load signal data
      signal_list <- readRDS(signal_file)

      if (length(signal_list) == 0) {
        cli_log(sprintf("Skipping empty signal file: %s", basename(signal_file)), "WARNING")
        next
      }

      #####################################################
      # CREATE TAIL FEATURES LIST
      #####################################################
      cli_log("Creating tail features list for polyT...", "INFO", bullet = TRUE)

      tail_feature_list <- tryCatch({
        invisible(utils::capture.output(
          result <- ninetails::create_tail_features_list_dorado(
            signal_list = signal_list,
            num_cores = num_cores
          )
        ))
        result
      }, error = function(e) {
        cli_log(sprintf("Error in polyT feature extraction: %s", e$message), "ERROR")
        stop(e$message, call. = FALSE)
      })

      cli_log(sprintf("Created features for %d polyT reads", length(tail_feature_list)), "INFO")

      #####################################################
      # CREATE TAIL CHUNK LIST
      #####################################################
      cli_log("Creating tail segmentation data for polyT...", "INFO", bullet = TRUE)

      tail_chunk_list <- tryCatch({
        invisible(utils::capture.output(
          result <- ninetails::create_tail_chunk_list_dorado(
            tail_feature_list = tail_feature_list,
            num_cores = num_cores
          )
        ))
        result
      }, error = function(e) {
        cli_log(sprintf("Error in polyT chunk creation: %s", e$message), "ERROR")
        stop(e$message, call. = FALSE)
      })

      cli_log(sprintf("Created chunks for %d polyT reads", length(tail_chunk_list)), "INFO")

      #####################################################
      # CREATE GAF LIST
      #####################################################
      cli_log("Computing gramian angular fields for polyT...", "INFO", bullet = TRUE)

      gaf_list <- tryCatch({
        invisible(utils::capture.output(
          result <- ninetails::create_gaf_list(
            tail_chunk_list = tail_chunk_list,
            num_cores = num_cores
          )
        ))
        result
      }, error = function(e) {
        cli_log(sprintf("Error in polyT GAF computation: %s", e$message), "ERROR")
        stop(e$message, call. = FALSE)
      })

      cli_log(sprintf("Created GAFs for %d polyT chunks", length(gaf_list)), "INFO")

      #####################################################
      # PREDICT CLASSES (TEMPORARY: USING POLYA MODEL)
      #####################################################
      cli_log("Running predictions with polyA model (temporary for polyT)...", "INFO", bullet = TRUE)

      # TODO: Replace with polyT-specific model when available
      # predicted_list <- ninetails::predict_gaf_classes_polyt(gaf_list)
      predicted_list <- tryCatch({
        invisible(utils::capture.output(
          result <- ninetails::predict_gaf_classes(gaf_list),  # Using polyA model temporarily
          type = "message"  # Captures message output from TensorFlow
        ))
        result
      }, error = function(e) {
        cli_log(sprintf("Error in polyT prediction: %s", e$message), "ERROR")
        stop(e$message, call. = FALSE)
      })

      cli_log(sprintf("Completed predictions for %d polyT GAFs", length(predicted_list$chunkname)), "INFO")

      #####################################################
      # CREATE OUTPUTS FOR THIS PART
      #####################################################
      cli_log("Creating outputs for polyT...", "INFO", bullet = TRUE)

      # Create a temporary summary for this subset of polyT reads
      temp_summary <- data.frame(
        readname = names(signal_list),
        polya_length = sapply(signal_list, length),  # Note: still called polya_length in output format
        qc_tag = "PASS"  # Assume all classified polyT reads pass QC
      )

      part_outputs <- tryCatch({
        invisible(utils::capture.output(
          result <- ninetails::create_outputs(
            tail_feature_list = tail_feature_list,
            tail_chunk_list = tail_chunk_list,
            nanopolish = temp_summary,  # Use temp summary
            predicted_list = predicted_list,
            num_cores = num_cores,
            pass_only = FALSE,  # We've already filtered
            qc = qc
          )
        ))
        result
      }, error = function(e) {
        cli_log(sprintf("Error in polyT output creation: %s", e$message), "ERROR")
        stop(e$message, call. = FALSE)
      })

      cli_log("PolyT output creation completed", "SUCCESS")

      # Store this part's outputs
      all_outputs[[i]] <- part_outputs

      # Save intermediate results
      output_basename <- sprintf("%spolyt_results_part%d.rds",
                                 ifelse(nchar(prefix) > 0, paste0(prefix, "_"), ""),
                                 i)
      output_file <- file.path(polyt_temp_dir, output_basename)
      saveRDS(part_outputs, output_file)

      cli_log(sprintf("Saved intermediate polyT results to %s", basename(output_file)), "INFO")

    }, error = function(e) {
      cli_log(sprintf("Error processing polyT signal file %d: %s", i, e$message), "ERROR")
    })

    # Explicit garbage collection
    gc()
  }

  ################################################################################
  # MERGE ALL POLYT RESULTS
  ################################################################################

  cli_log("Merging polyT processing results...", "INFO", "Merging PolyT Results")

  if (length(all_outputs) == 0) {
    stop("No polyT results were generated")
  }

  # Combine read classes from all parts
  all_read_classes <- do.call(rbind, lapply(all_outputs, function(x) {
    if (!is.null(x$read_classes)) {
      return(x$read_classes)
    } else {
      return(data.frame())
    }
  }))

  # Combine nonadenosine residues from all parts
  all_nonthymidine <- do.call(rbind, lapply(all_outputs, function(x) {
    if (!is.null(x$nonadenosine_residues)) {
      return(x$nonadenosine_residues)
    } else {
      return(data.frame())
    }
  }))

  # Add tail type information to results
  if (nrow(all_read_classes) > 0) {
    all_read_classes$tail_type <- "polyT"
  }

  if (nrow(all_nonthymidine) > 0) {
    all_nonthymidine$tail_type <- "polyT"
  }

  ################################################################################
  # COMPILE FINAL RESULTS
  ################################################################################

  # Calculate processing statistics
  processing_stats <- list(
    total_reads_processed = nrow(all_read_classes),
    decorated_reads = sum(all_read_classes$class == "decorated", na.rm = TRUE),
    blank_reads = sum(all_read_classes$class == "blank", na.rm = TRUE),
    unclassified_reads = sum(all_read_classes$class == "unclassified", na.rm = TRUE),
    total_modifications = nrow(all_nonthymidine),
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

  # Return combined results
  return(list(
    read_classes = all_read_classes,
    nonadenosine_residues = all_nonthymidine,  # Non-thymidine residues in polyT tails
    processing_stats = processing_stats,
    tail_type = "polyT"
  ))
}


################################################################################
# RESULTS MERGING FUNCTION FOR cDNA
################################################################################

#' Merge polyA and polyT processing results for cDNA analysis
#'
#' This function combines the results from separate polyA and polyT processing
#' paths into a unified output. It merges read classifications, modification
#' predictions, and processing statistics while preserving tail type information.
#' It also handles unidentified reads that could not be classified.
#'
#' @param polya_results List. Results from polyA processing (from process_polya_reads_cdna).
#' Can be NULL if no polyA reads were found.
#' @param polyt_results List. Results from polyT processing (from process_polyt_reads_cdna).
#' Can be NULL if no polyT reads were found.
#' @param unidentified_files Character vector. Paths to files containing reads
#' that could not be classified as polyA or polyT.
#' @param save_dir Character. Directory where merged results will be saved.
#' @param prefix Character. Optional prefix for output file names.
#' @param cli_log Function for logging messages and progress.
#'
#' @return List containing merged cDNA results:
#'   \describe{
#'     \item{read_classes}{Data frame with all read classifications including tail_type}
#'     \item{nonadenosine_residues}{Data frame with all predicted modifications including tail_type}
#'     \item{processing_stats}{Combined summary statistics from both processing paths}
#'     \item{unidentified_stats}{Statistics for unidentified reads}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' merged_results <- merge_cdna_results(
#'   polya_results = polya_output,
#'   polyt_results = polyt_output,
#'   unidentified_files = c("unidentified_part1.tsv"),
#'   save_dir = "path/to/output/",
#'   prefix = "experiment1",
#'   cli_log = message
#' )
#' }
merge_cdna_results <- function(polya_results = NULL,
                               polyt_results = NULL,
                               unidentified_files = character(0),
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

  # Create output directory for merged results
  merged_dir <- file.path(save_dir, "merged_results_dir")
  if (!dir.exists(merged_dir)) {
    dir.create(merged_dir, recursive = TRUE)
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
  # PROCESS UNIDENTIFIED READS
  ################################################################################

  cli_log("Processing unidentified reads...", "INFO", bullet = TRUE)

  # Process unidentified reads for statistics
  total_unidentified <- 0
  unidentified_read_ids <- character(0)

  if (length(unidentified_files) > 0) {
    # Validate unidentified files exist
    existing_files <- unidentified_files[file.exists(unidentified_files)]

    if (length(existing_files) < length(unidentified_files)) {
      missing_files <- unidentified_files[!file.exists(unidentified_files)]
      cli_log(sprintf("WARNING: %d unidentified files not found: %s",
                      length(missing_files), paste(basename(missing_files), collapse=", ")), "WARNING")
    }

    for (unidentified_file in existing_files) {
      tryCatch({
        unidentified_data <- vroom::vroom(unidentified_file, show_col_types = FALSE)
        if ("read_id" %in% colnames(unidentified_data)) {
          file_count <- nrow(unidentified_data)
          total_unidentified <- total_unidentified + file_count
          unidentified_read_ids <- c(unidentified_read_ids, unidentified_data$read_id)
          cli_log(sprintf("Processed %d unidentified reads from %s",
                          file_count, basename(unidentified_file)), "INFO")
        }
      }, error = function(e) {
        cli_log(sprintf("Error reading unidentified file %s: %s",
                        basename(unidentified_file), e$message), "WARNING")
      })
    }

    # Remove duplicates
    unidentified_read_ids <- unique(unidentified_read_ids)
    cli_log(sprintf("Total unique unidentified reads: %d", length(unidentified_read_ids)), "INFO")
  } else {
    cli_log("No unidentified read files to process", "INFO")
  }

  ################################################################################
  # CREATE COMBINED STATISTICS
  ################################################################################

  cli_log("Calculating combined statistics...", "INFO", bullet = TRUE)

  # Extract individual statistics
  polya_stats <- if (!is.null(polya_results$processing_stats)) polya_results$processing_stats else list()
  polyt_stats <- if (!is.null(polyt_results$processing_stats)) polyt_results$processing_stats else list()

  # Calculate combined statistics
  combined_stats <- list(
    # Overall totals
    total_reads_processed = nrow(merged_read_classes),
    total_modifications = nrow(merged_modifications),

    # PolyA statistics
    polya_reads_processed = if (!is.null(polya_stats$total_reads_processed)) polya_stats$total_reads_processed else 0,
    polya_decorated = if (!is.null(polya_stats$decorated_reads)) polya_stats$decorated_reads else 0,
    polya_blank = if (!is.null(polya_stats$blank_reads)) polya_stats$blank_reads else 0,
    polya_unclassified = if (!is.null(polya_stats$unclassified_reads)) polya_stats$unclassified_reads else 0,
    polya_modifications = if (!is.null(polya_stats$total_modifications)) polya_stats$total_modifications else 0,

    # PolyT statistics
    polyt_reads_processed = if (!is.null(polyt_stats$total_reads_processed)) polyt_stats$total_reads_processed else 0,
    polyt_decorated = if (!is.null(polyt_stats$decorated_reads)) polyt_stats$decorated_reads else 0,
    polyt_blank = if (!is.null(polyt_stats$blank_reads)) polyt_stats$blank_reads else 0,
    polyt_unclassified = if (!is.null(polyt_stats$unclassified_reads)) polyt_stats$unclassified_reads else 0,
    polyt_modifications = if (!is.null(polyt_stats$total_modifications)) polyt_stats$total_modifications else 0,

    # Unidentified statistics
    unidentified_reads = length(unidentified_read_ids),

    # Calculate percentages
    polya_percentage = 0,
    polyt_percentage = 0,
    unidentified_percentage = 0
  )

  # Calculate percentages of total reads
  grand_total <- combined_stats$polya_reads_processed + combined_stats$polyt_reads_processed + combined_stats$unidentified_reads

  if (grand_total > 0) {
    combined_stats$polya_percentage <- round(combined_stats$polya_reads_processed / grand_total * 100, 2)
    combined_stats$polyt_percentage <- round(combined_stats$polyt_reads_processed / grand_total * 100, 2)
    combined_stats$unidentified_percentage <- round(combined_stats$unidentified_reads / grand_total * 100, 2)
  }

  # Create unidentified statistics
  unidentified_stats <- list(
    total_unidentified = combined_stats$unidentified_reads,
    unidentified_files = unidentified_files,
    unidentified_percentage = combined_stats$unidentified_percentage
  )

  ################################################################################
  # SAVE INTERMEDIATE MERGED RESULTS
  ################################################################################

  cli_log("Saving intermediate merged results...", "INFO", bullet = TRUE)

  # Save merged read classes
  if (nrow(merged_read_classes) > 0) {
    read_classes_file <- file.path(merged_dir,
                                   sprintf("%smerged_read_classes.tsv",
                                           ifelse(nchar(prefix) > 0, paste0(prefix, "_"), "")))
    vroom::vroom_write(merged_read_classes, read_classes_file, delim = "\t")
    cli_log(sprintf("Saved merged read classes to %s", basename(read_classes_file)), "INFO")
  }

  # Save merged modifications
  if (nrow(merged_modifications) > 0) {
    modifications_file <- file.path(merged_dir,
                                    sprintf("%smerged_modifications.tsv",
                                            ifelse(nchar(prefix) > 0, paste0(prefix, "_"), "")))
    vroom::vroom_write(merged_modifications, modifications_file, delim = "\t")
    cli_log(sprintf("Saved merged modifications to %s", basename(modifications_file)), "INFO")
  }


  ################################################################################
  # FINAL LOGGING AND SUMMARY
  ################################################################################

  cli_log("Merging Summary", "INFO", "Merging Summary")
  cli_log(sprintf("Total reads processed: %d", grand_total), bullet = TRUE)
  cli_log(sprintf("- PolyA reads: %d (%.2f%%)",
                  combined_stats$polya_reads_processed, combined_stats$polya_percentage), bullet = TRUE)
  cli_log(sprintf("- PolyT reads: %d (%.2f%%)",
                  combined_stats$polyt_reads_processed, combined_stats$polyt_percentage), bullet = TRUE)
  cli_log(sprintf("- Unidentified reads: %d (%.2f%%)",
                  combined_stats$unidentified_reads, combined_stats$unidentified_percentage), bullet = TRUE)

  cli_log(sprintf("Total modifications detected: %d", combined_stats$total_modifications), bullet = TRUE)
  cli_log(sprintf("- PolyA modifications: %d", combined_stats$polya_modifications), bullet = TRUE)
  cli_log(sprintf("- PolyT modifications: %d", combined_stats$polyt_modifications), bullet = TRUE)

  if (combined_stats$polyt_modifications > 0) {
    cli_log("NOTE: PolyT modifications detected using temporary polyA model", "INFO", bullet = TRUE)
  }

  cli_log("cDNA results merging completed successfully", "SUCCESS")

  ################################################################################
  # RETURN MERGED RESULTS
  ################################################################################

  return(list(
    read_classes = merged_read_classes,
    nonadenosine_residues = merged_modifications,
    processing_stats = combined_stats,
    unidentified_stats = unidentified_stats
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

















