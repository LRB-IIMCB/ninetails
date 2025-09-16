#' Extract poly(A) information from BAM file
#'
#' Extracts poly(A) tail coordinates and related information from a BAM file using a summary file for pod5 file mapping.
#'
#' @param bam_file Character path to the BAM file containing aligned reads with poly(A) information.
#'
#' @param summary_file Character path to the summary file containing read_id and filename mapping.
#'
#' @param cli_log Function for logging messages, defaults to message.
#'
#' @return A data frame containing poly(A) information with columns:
#'   \describe{
#'     \item{read_id}{Unique read identifier}
#'     \item{filename}{Corresponding pod5 file name}
#'     \item{reference}{Reference sequence name}
#'     \item{ref_start}{Alignment start position}
#'     \item{ref_end}{Alignment end position}
#'     \item{mapq}{Mapping quality score}
#'     \item{poly_tail_length}{Length of poly(A) tail}
#'     \item{poly_tail_start}{Poly(A) start position}
#'     \item{poly_tail_end}{Poly(A) end position}
#'     \item{poly_tail2_start}{Secondary poly(A) start position}
#'     \item{poly_tail2_end}{Secondary poly(A) end position}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' polya_data <- extract_polya_from_bam(
#'   bam_file = "path/to/aligned.bam",
#'   summary_file = "path/to/summary.txt"
#' )
#' }
extract_polya_from_bam <- function(bam_file, summary_file, cli_log = message) {

  # Check for required Bioconductor packages
  if (!requireNamespace("Rsamtools", quietly = TRUE)) {
    stop("Package 'Rsamtools' is required for extract_polya_from_bam(). Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("S4Vectors", quietly = TRUE)) {
    stop("Package 'S4Vectors' is required for extract_polya_from_bam(). Please install it.",
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
    tag = c("pt", "pa"),  # tags to extract
    what = c("qname", "rname", "pos", "mapq", "flag", "cigar"))

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
  ref_starts <- integer(initial_count)
  ref_ends <- integer(initial_count)
  mapqs <- integer(initial_count)
  poly_tail_lengths <- numeric(initial_count)
  #anchor_positions <- integer(initial_count)
  poly_tail_starts <- integer(initial_count)
  poly_tail_ends <- integer(initial_count)
  poly_tail2_starts <- integer(initial_count)
  poly_tail2_ends <- integer(initial_count)
  filenames <- character(initial_count)  # New vector for pod5 files

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
    poly_tail_start <- -1L
    poly_tail_end <- -1L
    poly_tail2_start <- -1L
    poly_tail2_end <- -1L

    # Add pa tag information if it exists and is valid
    if (!is.null(bam_data$tag$pa[[i]]) && length(bam_data$tag$pa[[i]]) == 5) {
      pa_values <- bam_data$tag$pa[[i]]
      if (!any(is.na(pa_values))) {
        #anchor_pos <- pa_values[1]
        poly_tail_start <- pa_values[2]
        poly_tail_end <- pa_values[3]
        poly_tail2_start <- pa_values[4]
        poly_tail2_end <- pa_values[5]
      }
    }

    # Get pod5 file from lookup
    current_read_id <- bam_data$qname[i]
    filename <- pod5_lookup[current_read_id]
    if (is.na(filename)) {
      next  # Skip if no matching pod5 file found
    }

    # Increment counter and store values
    valid_count <- valid_count + 1
    read_ids[valid_count] <- current_read_id
    references[valid_count] <- as.character(bam_data$rname[i])
    ref_starts[valid_count] <- bam_data$pos[i]
    ref_ends[valid_count] <- bam_data$pos[i] + nchar(bam_data$cigar[i])
    mapqs[valid_count] <- bam_data$mapq[i]
    poly_tail_lengths[valid_count] <- poly_tail_length
    #anchor_positions[valid_count] <- anchor_pos
    poly_tail_starts[valid_count] <- poly_tail_start
    poly_tail_ends[valid_count] <- poly_tail_end
    poly_tail2_starts[valid_count] <- poly_tail2_start
    poly_tail2_ends[valid_count] <- poly_tail2_end
    filenames[valid_count] <- filename
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
    ref_starts <- ref_starts[1:valid_count]
    ref_ends <- ref_ends[1:valid_count]
    mapqs <- mapqs[1:valid_count]
    poly_tail_lengths <- poly_tail_lengths[1:valid_count]
    #anchor_positions <- anchor_positions[1:valid_count]
    poly_tail_starts <- poly_tail_starts[1:valid_count]
    poly_tail_ends <- poly_tail_ends[1:valid_count]
    poly_tail2_starts <- poly_tail2_starts[1:valid_count]
    poly_tail2_ends <- poly_tail2_ends[1:valid_count]
    filenames <- filenames[1:valid_count]

    # Check for duplicates
    duplicate_reads <- duplicated(read_ids)
    if (any(duplicate_reads)) {
      cli_log(sprintf("WARNING: Found %d duplicate read names, keeping first occurrence only",
                      sum(duplicate_reads)), "WARNING")
      # Keep only unique entries
      keep_idx <- !duplicate_reads
      read_ids <- read_ids[keep_idx]
      references <- references[keep_idx]
      ref_starts <- ref_starts[keep_idx]
      ref_ends <- ref_ends[keep_idx]
      mapqs <- mapqs[keep_idx]
      poly_tail_lengths <- poly_tail_lengths[keep_idx]
      #anchor_positions <- anchor_positions[keep_idx]
      poly_tail_starts <- poly_tail_starts[keep_idx]
      poly_tail_ends <- poly_tail_ends[keep_idx]
      poly_tail2_starts <- poly_tail2_starts[keep_idx]
      poly_tail2_ends <- poly_tail2_ends[keep_idx]
      filenames <- filenames[keep_idx]
    }

    # Create the final data frame
    polya_df <- tibble::tibble(
      read_id = read_ids,
      filename = filenames,  # Added pod5 file information
      reference = references,
      ref_start = ref_starts,
      ref_end = ref_ends,
      mapq = mapqs,
      poly_tail_length = poly_tail_lengths,
      #anchor_pos = anchor_positions,
      poly_tail_start = poly_tail_starts,
      poly_tail_end = poly_tail_ends,
      poly_tail2_start = poly_tail2_starts,
      poly_tail2_end = poly_tail2_ends
    )

    cli_log("Processing complete", "SUCCESS")
    return(polya_df)
  } else {
    cli_log("No reads with poly(A) information found after filtering", "WARNING")
    return(tibble::tibble())
  }
}
