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

  assert_file_exists(bam_file, "BAM")
  assert_file_exists(dorado_summary, "Summary")
  assert_dir_exists(save_dir, "Output")

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
    read_sequence <- as.character(bam_data$seq[i]) # basecalled sequences
    # BAM stores reads aligned to the reverse strand as their reverse complement.
    # Restore the sequencing orientation, otherwise the primer matching in
    # detect_orientation_single() sees every aligned read in the same orientation
    if (bitwAnd(bam_data$flag[i], 0x10)) {
      read_sequence <- chartr("ACGTacgt", "TGCAtgca",
                              paste(rev(strsplit(read_sequence, "")[[1]]), collapse = ""))
    }
    seqs[valid_count] <- read_sequence

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
  assert_condition(is.numeric(num_cores) && num_cores > 0,
                   "Number of cores must be a positive numeric value")

  assert_condition(is.character(pod5_dir),
                   "Pod5 directory path must be a character string")
  assert_dir_exists(pod5_dir, "Pod5 files")

  assert_condition(is.character(save_dir),
                   "Output directory path must be a character string")

  assert_condition(is.logical(qc),
                   "qc must be logical [TRUE/FALSE]")

  # Validate prefix (if provided)
  if (nchar(prefix) > 0) {
    assert_condition(is.character(prefix),
                     "File name prefix must be a character string")
  }

  # Validate BAM input
  if (is_string(bam_file)) {
    assert_file_exists(bam_file, "BAM")
  } else {
    stop("BAM file must be provided as a character string path", call. = FALSE)
  }

  # Validate sequencing summary input
  if (is_string(dorado_summary)) {
    assert_file_exists(dorado_summary, "Dorado summary")
  } else {
    assert_condition(is.data.frame(dorado_summary) && nrow(dorado_summary) > 0,
                     "Dorado summary must be a non-empty data frame or valid file path")
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

    # Extract raw flanked tail signals (no winsorization/downsampling: the
    # cDNA models work on the raw signal)
    signal_list <- ninetails::extract_tail_signals_cdna(
      tail_coordinates = ninetails::read_tail_coordinates_cdna(summary_data, source = "dorado"),
      pod5_dir = pod5_dir,
      num_cores = num_cores
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
    warning(
      "detect_orientation_single failed on a sequence - returning 'unknown': ",
      conditionMessage(e)
    )
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
#' The augmented table (with tail_type) is also written to a sibling
#' "<name>_classified.tsv" file for persistence; the original input files are
#' left untouched.
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

  assert_condition(is.character(sequence_files),
                   "sequence_files must be a character vector of file paths")
  assert_condition(all(file.exists(sequence_files)),
                   "All sequence files must exist")
  assert_condition(is.numeric(num_cores) && num_cores > 0,
                   "Number of cores must be a positive numeric value")

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

      # Persist results WITHOUT mutating the input: write the augmented table
      # (with tail_type) to a sibling "<name>_classified.tsv". Overwriting the
      # original sequence file risks corrupting user inputs on partial failure
      # and makes re-runs non-idempotent (re-classifying an already-tagged file)
      classified_file <- file.path(
        dirname(current_file),
        paste0(tools::file_path_sans_ext(basename(current_file)), "_classified.tsv")
      )
      vroom::vroom_write(sequence_data, classified_file, delim = "\t")
      cli_log(sprintf("Wrote classified sequences to %s", basename(classified_file)), "INFO")

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
# POLYA/POLYT READS PROCESSING FUNCTION FOR cDNA
################################################################################

# Internal: Lachesis QC tags other than PASS mark unreliable tail coordinates.
# Dorado input carries the mapping quality (a number) in qc_tag and is never
# rejected here.
.is_qc_failed_cdna <- function(qc_tag) {
  !is.na(qc_tag) & grepl("[A-Za-z]", qc_tag) & qc_tag != "PASS"
} # .is_qc_failed_cdna


# Internal: shared body of process_polya_reads_cdna() and
# process_polyt_reads_cdna(); the two orientations differ only in the model
# and calling rule, which are selected by tail_type.
.process_reads_cdna <- function(sequences,
                                signal_files,
                                tail_type,
                                sequences_arg,
                                num_cores,
                                qc,
                                save_dir,
                                prefix,
                                cli_log) {

  assert_condition(is.data.frame(sequences),
                   paste0(sequences_arg, " must be a data frame"))
  assert_condition("read_id" %in% colnames(sequences),
                   paste0(sequences_arg, " must contain a 'read_id' column"))
  assert_condition(is.character(signal_files),
                   "signal_files must be a character vector of file paths")
  assert_condition(all(file.exists(signal_files)),
                   "All signal files must exist")

  cli_log(sprintf("Starting %s reads processing", tail_type), "INFO",
          sprintf("%s Processing", tail_type))

  #coerce tibble to df
  coordinates <- as.data.frame(sequences)
  coordinates <- coordinates[!duplicated(coordinates$read_id), , drop = FALSE]

  # columns absent from the input are reported as NA
  for (column_name in c("contig", "qc_tag")) {
    if (!column_name %in% colnames(coordinates)) {
      coordinates[[column_name]] <- rep(NA_character_, nrow(coordinates))
    }
  }
  coordinates$tail_type <- rep(tail_type, nrow(coordinates))

  has_coordinates <- all(c("tail_start", "tail_end", "polya_length") %in% colnames(coordinates))

  # reads eligible for classification: valid coordinates, tail of at least
  # 10 nt and (with qc) Lachesis QC tag PASS
  if (has_coordinates) {
    eligible <- !is.na(coordinates$tail_start) & !is.na(coordinates$tail_end) &
      !is.na(coordinates$polya_length) &
      coordinates$tail_start > 0 & coordinates$tail_end > coordinates$tail_start &
      coordinates$polya_length >= 10
    if (qc) {
      eligible <- eligible & !.is_qc_failed_cdna(coordinates$qc_tag)
    }
  } else {
    eligible <- rep(TRUE, nrow(coordinates))
  }
  eligible_read_ids <- coordinates$read_id[eligible]

  ################################################################################
  # INTERMEDIATE DIRECTORIES (DRS-STYLE: RESULTS OF EACH PART ARE KEPT ON DISK)
  ################################################################################

  file_tag <- ifelse(tail_type == "polyA", "polya", "polyt")
  file_prefix <- ifelse(nchar(prefix) > 0, paste0(prefix, "_"), "")

  processing_dir <- file.path(save_dir, paste0(file_tag, "_processing_dir"))
  nonA_temp_dir <- file.path(processing_dir, "nonA_temp_dir")
  polya_chunks_dir <- file.path(processing_dir, "polya_chunks_dir")
  dorado_summary_dir <- file.path(processing_dir, "dorado_summary_dir")

  for (dir in c(nonA_temp_dir, polya_chunks_dir, dorado_summary_dir)) {
    # leftovers of a previous run would be merged into the outputs
    unlink(dir, recursive = TRUE)
    dir.create(dir, recursive = TRUE)
  }

  # tail coordinates of all reads of this orientation
  vroom::vroom_write(coordinates[, intersect(c("read_id", "contig", "tail_start", "tail_end",
                                               "polya_length", "qc_tag", "tail_type"),
                                             colnames(coordinates))],
                     file.path(dorado_summary_dir,
                               sprintf("%s%s_summary_part1.txt", file_prefix, file_tag)),
                     delim = "\t")

  ################################################################################
  # CLASSIFY SIGNALS (ONE SIGNAL FILE AT A TIME TO LIMIT MEMORY USAGE)
  ################################################################################

  n_classified <- 0

  for (i in seq_along(signal_files)) {
    signal_file <- signal_files[i]
    cli_log(sprintf("Processing signal file %d/%d: %s",
                    i, length(signal_files), basename(signal_file)), "INFO")

    signal_list <- readRDS(signal_file)
    signal_list <- signal_list[names(signal_list) %in% eligible_read_ids]

    if (length(signal_list) == 0) {
      next
    }

    if (!has_coordinates) {
      stop(sprintf("Required columns missing: %s",
                   paste(setdiff(c("tail_start", "tail_end", "polya_length"), colnames(coordinates)),
                         collapse = ", ")),
           call. = FALSE)
    }

    nonA_calls <- ninetails::predict_nonA_cdna(
      signal_list = signal_list,
      tail_coordinates = coordinates,
      tail_type = tail_type
    )

    # calls of this part and IDs of the reads which were classified go to disk
    saveRDS(list(nonA_calls = nonA_calls,
                 classified_reads = names(signal_list)),
            file.path(nonA_temp_dir,
                      sprintf("%s%s_nonA_part%d.rds", file_prefix, file_tag, i)))
    n_classified <- n_classified + length(signal_list)

    rm(signal_list, nonA_calls)
    gc()
  } # for signal_file

  if (n_classified == 0) {
    cli_log(sprintf("No %s signals were found for processing", tail_type), "WARNING")
    return(list(
      read_classes = data.frame(),
      nonadenosine_residues = data.frame(),
      processing_stats = list(total_reads_processed = 0),
      tail_type = tail_type
    )) # No signals, return empty result
  }

  ################################################################################
  # CREATE OUTPUTS
  ################################################################################

  outputs <- ninetails::create_outputs_dorado_cdna(
    dorado_summary_dir = dorado_summary_dir,
    nonA_temp_dir = nonA_temp_dir,
    polya_chunks_dir = polya_chunks_dir,
    num_cores = num_cores,
    qc = qc
  )

  processing_stats <- list(
    total_reads_processed = nrow(outputs$read_classes),
    decorated_reads = sum(outputs$read_classes$class == "decorated", na.rm = TRUE),
    blank_reads = sum(outputs$read_classes$class == "blank", na.rm = TRUE),
    unclassified_reads = sum(outputs$read_classes$class == "unclassified", na.rm = TRUE),
    total_modifications = nrow(outputs$nonadenosine_residues),
    tail_type = tail_type
  )

  cli_log(sprintf("%s Processing Summary", tail_type), "INFO", sprintf("%s Summary", tail_type))
  cli_log(sprintf("Total %s reads processed: %d", tail_type, processing_stats$total_reads_processed), bullet = TRUE)
  cli_log(sprintf("Decorated reads: %d", processing_stats$decorated_reads), bullet = TRUE)
  cli_log(sprintf("Blank reads: %d", processing_stats$blank_reads), bullet = TRUE)
  cli_log(sprintf("Unclassified reads: %d", processing_stats$unclassified_reads), bullet = TRUE)
  cli_log(sprintf("Total modifications detected: %d", processing_stats$total_modifications), bullet = TRUE)

  result <- list(
    read_classes = outputs$read_classes,
    nonadenosine_residues = outputs$nonadenosine_residues,
    processing_stats = processing_stats,
    tail_type = tail_type
  )

  return(result)
} # .process_reads_cdna


################################################################################
#  PROCESSING WITH POLY(A) MODEL
###############################################################################
#' Process poly(A) reads with the cDNA poly(A) model
#'
#' Classifies the tail signals of reads sequenced in the poly(A) orientation
#' (tail at the 3' end of the read) with the dedicated cDNA model and creates
#' the ninetails output tables for them.
#'
#' @details
#' Signals are classified with \code{\link{predict_nonA_cdna}}, one signal
#' file at a time. The calls of each file are written to disk, so memory usage
#' does not grow with the size of the sequencing run. The tables are then
#' built from these files by \code{\link{create_outputs_dorado_cdna}}. Reads with invalid
#' coordinates, tails shorter than 10 nt or (with \code{qc = TRUE}) a
#' Lachesis QC tag other than \code{PASS} are not classified.
#'
#' @param polya_sequences Data frame/tibble. Poly(A) reads with their tail
#' coordinates: \code{read_id}, \code{tail_start}, \code{tail_end},
#' \code{polya_length} and optionally \code{contig} and \code{qc_tag}
#' (see \code{\link{read_tail_coordinates_cdna}}).
#' @param signal_files Character vector. Paths to the \code{.rds} files with
#' raw flanked tail signals (see \code{\link{extract_tail_signals_cdna}}).
#' @param num_cores Integer. Number of CPU cores (currently unused; the
#' classification is parallelised by TensorFlow).
#' @param qc Logical. Whether to apply quality control filtering.
#' @param save_dir Character. Output directory. Intermediate files are
#' written to its \code{polya_processing_dir} (or \code{polyt_processing_dir})
#' subdirectory.
#' @param prefix Character. Optional prefix for intermediate file names.
#' @param cli_log Function for logging messages and progress.
#'
#' @return List containing poly(A) processing results:
#'   \describe{
#'     \item{read_classes}{Data frame with read classification results}
#'     \item{nonadenosine_residues}{Data frame with predicted non-adenosines}
#'     \item{processing_stats}{Summary statistics for poly(A) processing}
#'     \item{tail_type}{Character indicating this is "polyA" data}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' polya_results <- ninetails::process_polya_reads_cdna(
#'   polya_sequences = polya_coordinates,
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

  result <- .process_reads_cdna(sequences = polya_sequences,
                                signal_files = signal_files,
                                tail_type = "polyA",
                                sequences_arg = "polya_sequences",
                                num_cores = num_cores,
                                qc = qc,
                                save_dir = save_dir,
                                prefix = prefix,
                                cli_log = cli_log)

  return(result)
}

################################################################################
#  PROCESSING WITH POLY(T) MODEL
###############################################################################
#' Process poly(T) reads with the cDNA poly(T) model
#'
#' Classifies the tail signals of reads sequenced in the poly(T) orientation
#' (tail at the 5' end of the read, read as its complement) with the
#' dedicated cDNA model and creates the ninetails output tables for them.
#' Predictions refer to the residues of the original RNA tail (C, G, U), so
#' they are directly comparable with the poly(A) orientation.
#'
#' @inherit process_polya_reads_cdna details
#'
#' @param polyt_sequences Data frame/tibble. Poly(T) reads with their tail
#' coordinates: \code{read_id}, \code{tail_start}, \code{tail_end},
#' \code{polya_length} and optionally \code{contig} and \code{qc_tag}
#' (see \code{\link{read_tail_coordinates_cdna}}).
#' @inheritParams process_polya_reads_cdna
#'
#' @return List containing poly(T) processing results:
#'   \describe{
#'     \item{read_classes}{Data frame with read classification results}
#'     \item{nonadenosine_residues}{Data frame with predicted non-adenosines}
#'     \item{processing_stats}{Summary statistics for poly(T) processing}
#'     \item{tail_type}{Character indicating this is "polyT" data}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' polyt_results <- ninetails::process_polyt_reads_cdna(
#'   polyt_sequences = polyt_coordinates,
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
    stop("Poly(T) sequences are missing. Please provide valid polyt_sequences argument.",
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

  result <- .process_reads_cdna(sequences = polyt_sequences,
                                signal_files = signal_files,
                                tail_type = "polyT",
                                sequences_arg = "polyt_sequences",
                                num_cores = num_cores,
                                qc = qc,
                                save_dir = save_dir,
                                prefix = prefix,
                                cli_log = cli_log)

  return(result)
}



################################################################################
# cDNA-SPECIFIC OUTPUT CREATION FUNCTION
################################################################################

#' Create ninetails output tables for the cDNA pipeline
#'
#' Combines the tail coordinates with the non-adenosine calls of
#' \code{\link{predict_nonA_cdna}} into the two standard ninetails outputs:
#' per-read classes and non-adenosine residues with their estimated positions.
#' Inputs are read from the intermediate directories written by
#' \code{\link{process_polya_reads_cdna}} or
#' \code{\link{process_polyt_reads_cdna}}.
#'
#' @section Read classes:
#' \itemize{
#'   \item \code{decorated}, \code{YAY}: at least one non-adenosine residue
#'   \item \code{blank}, \code{MAU}: no non-adenosine residue detected
#'   \item \code{blank}, \code{MPU}: only terminal calls, discarded by
#'     \code{qc} (closer than 2 nt to either end of the tail)
#'   \item \code{unclassified}, \code{IRL}: tail shorter than 10 nt
#'   \item \code{unclassified}, \code{BAC}: invalid tail coordinates or no
#'     signal available
#'   \item \code{unclassified}, \code{QCF}: Lachesis QC tag other than
#'     \code{PASS} (only with \code{qc = TRUE})
#' }
#'
#' @param dorado_summary_dir Character string. Path to a directory with the
#' tail coordinate tables (.txt, .tsv or .csv) of the reads of one
#' orientation. Required columns: \code{read_id}, \code{contig},
#' \code{polya_length}, \code{qc_tag} and \code{tail_type}.
#' @param nonA_temp_dir Character string. Path to a directory with the
#' per-part RDS files of non-adenosine calls. Each file holds a list with
#' \code{nonA_calls} (output of \code{\link{predict_nonA_cdna}}) and
#' \code{classified_reads} (IDs of the reads whose signals were classified).
#' @param polya_chunks_dir Character string. Kept for compatibility with the
#' DRS pipeline (currently unused; the cDNA model does not produce chunk
#' files).
#' @param num_cores Integer [1]. Number of cores used for loading the files.
#' @param qc Logical [TRUE]. Whether to discard terminal calls and reads
#' failing the Lachesis QC.
#'
#' @return A named list:
#' \describe{
#'   \item{read_classes}{Data frame with columns \code{readname},
#'     \code{contig}, \code{polya_length}, \code{qc_tag}, \code{class},
#'     \code{comments}, \code{tail_type}}
#'   \item{nonadenosine_residues}{Data frame with columns \code{readname},
#'     \code{contig}, \code{prediction} (residue in the RNA tail),
#'     \code{est_nonA_pos} (position from the 3' end of the tail),
#'     \code{polya_length}, \code{qc_tag}, \code{tail_type}}
#' }
#'
#' @seealso \code{\link{predict_nonA_cdna}},
#'   \code{\link{merge_cdna_results}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' outputs <- ninetails::create_outputs_dorado_cdna(
#'   dorado_summary_dir = "polya_processing_dir/dorado_summary_dir",
#'   nonA_temp_dir = "polya_processing_dir/nonA_temp_dir",
#'   polya_chunks_dir = "polya_processing_dir/polya_chunks_dir",
#'   num_cores = 4,
#'   qc = TRUE
#' )
#' }
create_outputs_dorado_cdna <- function(dorado_summary_dir,
                                       nonA_temp_dir,
                                       polya_chunks_dir,
                                       num_cores = 1,
                                       qc = TRUE) {

  # Assertions
  if (missing(dorado_summary_dir)) stop("Dorado summary directory is missing.", call. = FALSE)
  if (missing(nonA_temp_dir)) stop("Non-A predictions directory is missing.", call. = FALSE)
  if (missing(polya_chunks_dir)) stop("Poly(A) chunks directory is missing.", call. = FALSE)

  assert_condition(is.numeric(num_cores) && num_cores > 0,
                   "Number of cores must be a positive integer")
  assert_condition(is.logical(qc),
                   "QC must be a logical value")
  assert_dir_exists(dorado_summary_dir, "Dorado summary")
  assert_dir_exists(nonA_temp_dir, "Non-A predictions")

  ################################################################################
  # LOAD DATA FILES
  ################################################################################

  dorado_summary_files <- list.files(dorado_summary_dir, pattern = "\\.(txt|tsv|csv)$", full.names = TRUE)

  if (length(dorado_summary_files) == 0) {
    stop("No summary files found in dorado_summary_dir", call. = FALSE)
  }

  tail_coordinates <- do.call(rbind, parallel::mclapply(dorado_summary_files, function(file) {
    as.data.frame(vroom::vroom(file, show_col_types = FALSE,
                               col_types = vroom::cols(read_id = "c", contig = "c",
                                                       qc_tag = "c", tail_type = "c")))
  }, mc.cores = num_cores))

  required_cols <- c("read_id", "contig", "polya_length", "qc_tag", "tail_type")
  missing_cols <- setdiff(required_cols, colnames(tail_coordinates))
  if (length(missing_cols) > 0) {
    stop(sprintf("Required columns missing: %s", paste(missing_cols, collapse = ", ")),
         call. = FALSE)
  }

  # per-part calls and IDs of the classified reads
  nonA_files <- list.files(nonA_temp_dir, pattern = "\\.rds$", full.names = TRUE)
  nonA_parts <- parallel::mclapply(nonA_files, readRDS, mc.cores = num_cores)

  nonA_calls <- do.call(rbind, lapply(nonA_parts, function(x) x$nonA_calls))
  classified_reads <- unlist(lapply(nonA_parts, function(x) x$classified_reads))
  rm(nonA_parts)

  if (is.null(nonA_calls)) {
    nonA_calls <- data.frame()
  }

  read_classes <- as.data.frame(tail_coordinates)[, required_cols]
  read_classes <- read_classes[!duplicated(read_classes$read_id), ]

  # Non-adenosine residues
  ##############################################################################
  nonadenosine_residues <- data.frame(
    readname = character(), contig = character(), prediction = character(),
    est_nonA_pos = numeric(), polya_length = numeric(), qc_tag = character(), tail_type = character(),
    stringsAsFactors = FALSE
  )
  moved_blank_readnames <- character(0)

  if (nrow(nonA_calls) > 0) {
    nonadenosine_residues <- merge(nonA_calls[, c("read_id", "prediction", "est_nonA_pos")],
                                   read_classes,
                                   by = "read_id")

    # terminal calls are unreliable (tail boundaries are estimates)
    if (qc == TRUE) {
      terminal_mask <- (nonadenosine_residues$est_nonA_pos < 2) |
        (nonadenosine_residues$est_nonA_pos > nonadenosine_residues$polya_length - 2)
      discarded_read_ids <- unique(nonadenosine_residues$read_id[terminal_mask])
      nonadenosine_residues <- nonadenosine_residues[!terminal_mask, ]
      # reads left without any call after the terminal filter
      moved_blank_readnames <- setdiff(discarded_read_ids, nonadenosine_residues$read_id)
    }

    names(nonadenosine_residues)[names(nonadenosine_residues) == "read_id"] <- "readname"
    nonadenosine_residues <- nonadenosine_residues[, c("readname", "contig", "prediction", "est_nonA_pos",
                                                       "polya_length", "qc_tag", "tail_type")]
  }

  # Read classes
  ##############################################################################
  read_classes$class <- "blank"
  read_classes$comments <- "MAU"

  read_classes$comments[read_classes$read_id %in% moved_blank_readnames] <- "MPU"

  decorated <- read_classes$read_id %in% nonadenosine_residues$readname
  read_classes$class[decorated] <- "decorated"
  read_classes$comments[decorated] <- "YAY"

  # reads which were not classified
  not_classified <- !read_classes$read_id %in% classified_reads
  read_classes$class[not_classified] <- "unclassified"
  read_classes$comments[not_classified] <- "BAC"

  if (qc == TRUE) {
    qc_failed <- not_classified & .is_qc_failed_cdna(read_classes$qc_tag)
    read_classes$comments[qc_failed] <- "QCF"
  }

  short_tail <- not_classified & !is.na(read_classes$polya_length) & read_classes$polya_length < 10
  read_classes$comments[short_tail] <- "IRL"

  names(read_classes)[names(read_classes) == "read_id"] <- "readname"
  read_classes <- read_classes[, c("readname", "contig", "polya_length", "qc_tag", "class", "comments", "tail_type")]

  ninetails_output <- list(
    read_classes = read_classes,
    nonadenosine_residues = nonadenosine_residues
  )

  return(ninetails_output)
}



################################################################################
# RESULTS MERGING FUNCTION FOR cDNA (UPDATED FOR SIMPLIFIED APPROACH)
################################################################################

#' Merge polyA and polyT processing results for cDNA analysis
#'
#' This function combines the results from separate poly(A) and poly(T) processing
#' paths into a unified output in standard ninetails format. Updated to work
#' with the simplified approach using tibbles instead of files.
#'
#' @param polya_results List. Results from poly(A) processing (from process_polya_reads_cdna).
#' Can be NULL if no polyA reads were found.
#'
#' @param polyt_results List. Results from poly(T) processing (from process_polyt_reads_cdna).
#' Can be NULL if no poly(T) reads were found.
#'
#' @param unidentified_reads Data frame/tibble. Reads that could not be classified
#' as polyA or polyT (with tail_type = "unidentified").
#'
#' @param save_dir Character. Directory where merged results will be saved.
#'
#' @param prefix Character. Optional prefix for output file names.
#'
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

  assert_condition(is.character(save_dir),
                   "save_dir must be a character string")

  cli_log("Starting cDNA results merging", "INFO", "Merging cDNA Results")

  # Check if we have any results to merge
  if (is.null(polya_results) && is.null(polyt_results)) {
    cli_log("WARNING: No poly(A) or poly(T) results provided for merging", "WARNING")
  }

  ################################################################################
  # MERGE READ CLASSIFICATIONS
  ################################################################################

  cli_log("Merging read classifications...", "INFO", bullet = TRUE)

  # Initialize empty data frame with expected structure
  merged_read_classes <- data.frame()

  # Add poly(A) read classes
  if (!is.null(polya_results) && !is.null(polya_results$read_classes)) {
    polya_classes <- polya_results$read_classes

    # Ensure tail_type column exists
    if (!"tail_type" %in% colnames(polya_classes) && nrow(polya_classes) > 0) {
      polya_classes$tail_type <- "polyA"
    }

    merged_read_classes <- dplyr::bind_rows(merged_read_classes, polya_classes)
    cli_log(sprintf("Added %d poly(A) read classifications", nrow(polya_classes)), "INFO")
  } else {
    cli_log("No poly(A) read classifications to merge", "INFO")
  }

  # Add poly(T0 read classes
  if (!is.null(polyt_results) && !is.null(polyt_results$read_classes)) {
    polyt_classes <- polyt_results$read_classes

    # Make sure that tail_type column exists
    if (!"tail_type" %in% colnames(polyt_classes) && nrow(polyt_classes) > 0) {
      polyt_classes$tail_type <- "polyT"
    }

    merged_read_classes <- dplyr::bind_rows(merged_read_classes, polyt_classes)
    cli_log(sprintf("Added %d poly(T) read classifications", nrow(polyt_classes)), "INFO")
  } else {
    cli_log("No poly(T) read classifications to merge", "INFO")
  }

  ################################################################################
  # MERGE NONADENOSINE RESIDUE PREDICTIONS
  ################################################################################

  cli_log("Merging nonA residue predictions...", "INFO", bullet = TRUE)

  # Initialize empty data frame for nonAs
  merged_modifications <- data.frame()

  # Add polyA nonAs
  if (!is.null(polya_results) && !is.null(polya_results$nonadenosine_residues)) {
    polya_mods <- polya_results$nonadenosine_residues

    # Ensure tail_type column exists
    if (!"tail_type" %in% colnames(polya_mods) && nrow(polya_mods) > 0) {
      polya_mods$tail_type <- "polyA"
    }

    merged_modifications <- dplyr::bind_rows(merged_modifications, polya_mods)
    cli_log(sprintf("Added %d poly(A) nonA residues", nrow(polya_mods)), "INFO")
  } else {
    cli_log("No poly(A) nonA residues to merge", "INFO")
  }

  # Add poly(T) nonAs
  if (!is.null(polyt_results) && !is.null(polyt_results$nonadenosine_residues)) {
    polyt_mods <- polyt_results$nonadenosine_residues

    # Make sure that tail_type column exists
    if (!"tail_type" %in% colnames(polyt_mods) && nrow(polyt_mods) > 0) {
      polyt_mods$tail_type <- "polyT"
    }

    merged_modifications <- dplyr::bind_rows(merged_modifications, polyt_mods)
    cli_log(sprintf("Added %d poly(T) nonA residues", nrow(polyt_mods)), "INFO")
  } else {
    cli_log("No poly(T) nonA residues to merge", "INFO")
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
  cli_log(sprintf("- Poly(A) reads processed: %d (%.2f%%)",
                  polya_count, if (grand_total > 0) polya_count / grand_total * 100 else 0), bullet = TRUE)
  cli_log(sprintf("- Poly(T) reads processed: %d (%.2f%%)",
                  polyt_count, if (grand_total > 0) polyt_count / grand_total * 100 else 0), bullet = TRUE)
  cli_log(sprintf("- Unidentified reads: %d (%.2f%%)",
                  unidentified_count, if (grand_total > 0) unidentified_count / grand_total * 100 else 0), bullet = TRUE)

  cli_log(sprintf("Total nonAs detected: %d", total_modifications), bullet = TRUE)

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
save_cdna_outputs <- function(outputs,
                              save_dir,
                              prefix = "") {

  # Input validation
  if (missing(outputs)) {
    stop("Outputs are missing. Please provide valid outputs argument.", call. = FALSE)
  }

  if (missing(save_dir)) {
    stop("Output directory is missing. Please provide a valid save_dir argument.", call. = FALSE)
  }

  assert_condition(is.list(outputs),
                   "outputs must be a list from merge_cdna_results()")

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
