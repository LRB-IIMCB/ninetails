#' Extract poly(A) tail information from BAM file
#'
#' This function extracts poly(A) tail information from a BAM file produced by
#' Dorado basecaller.
#'
#' - Requires BAM files produced by dorado version >= 1.0.0 for poly(A) coordinate information
#' - The pa tag should contain 5 values: anchor position and poly(A) region coordinates
#' - Returns empty data frame if no valid reads are found
#' - Uses memory-efficient vector pre-allocation for large files
#'
#' The function performs the following steps:
#' 1. Validates input BAM and summary files
#' 2. Reads pod5 file mappings from summary file
#' 3. Extracts BAM tags (pt for poly(A) tail length, pa for positions)
#' 4. Filters out:
#'    - Unmapped reads
#'    - Non-primary alignments
#'    - Reads without poly(A) tags
#' 5. Combines BAM information with pod5 file mappings
#' 6. Handles duplicate reads by keeping first occurrence
#'
#' Progress and filtering statistics are reported during execution.
#'
#' @param bam_file character string. Full path to the BAM file containing Dorado
#' basecalling results. The BAM file must contain pt and pa tags for poly(A) tail
#' information.
#'
#' @param summary_file Character string. Path to the corresponding dorado
#' summary file containing read to pod5 file mappings.
#'
#' @param cli_log Function for logging. This function is encoded in main
#' pipeline wrapper. Its purpose is to provide neatly formatted & informative log file.
#'
#' @return A data frame containing poly(A) information with the following columns:
#' \describe{
#'   \item{read_id}{Character. Read identifier}
#'   \item{pod5_file}{Character. Path to the corresponding pod5 file}
#'   \item{reference}{Character. Reference sequence name}
#'   \item{ref_start}{Integer. Start position on reference}
#'   \item{ref_end}{Integer. End position on reference}
#'   \item{mapq}{Integer. Mapping quality score}
#'   \item{polya_length}{Numeric. Length of poly(A) tail}
#'   \item{anchor_pos}{Integer. Position of poly(A) anchor (-1 if not available)}
#'   \item{polya_start}{Integer. Start position of poly(A) region (-1 if not available)}
#'   \item{polya_end}{Integer. End position of poly(A) region (-1 if not available)}
#'   \item{secondary_polya_start}{Integer. Start position of secondary poly(A) region (-1 if not available)}
#'   \item{secondary_polya_end}{Integer. End position of secondary poly(A) region (-1 if not available)}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Extract poly(A) information from BAM file
#' polya_data <- extract_polya_from_bam(
#'   bam_file = "path/to/dorado_calls.bam",
#'   summary_file = "path/to/summary_part1.txt"
#' )
#'
#' # Check results
#' head(polya_data)
#' }
#'
extract_polya_from_bam <- function(bam_file, summary_file, cli_log) {

  # Assertions
  if (!file.exists(bam_file)) {
    stop("BAM file does not exist: ", bam_file)}
  if (!file.exists(summary_file)) {
    stop("Summary file does not exist: ", summary_file)}

  # Read summary file efficiently (only needed columns)
  cli_log(sprintf("Reading summary file: %s", base::basename(summary_file)), "INFO")
  summary_data <- data.table::fread(summary_file, select = c("filename", "read_id"))
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
    cli_log("ERROR: poly(A) coordinate information (pa tag) not found in BAM file.", "ERROR")
    cli_log("Please ensure data were basecalled with dorado version >= 1.0.0", "ERROR")
    stop("poly(A) coordinate information (pa tag) not found in BAM file")
  }

  # Pre-allocate vectors with maximum possible size
  read_ids <- character(initial_count)
  references <- character(initial_count)
  ref_starts <- integer(initial_count)
  ref_ends <- integer(initial_count)
  mapqs <- integer(initial_count)
  polya_lengths <- numeric(initial_count)
  anchor_positions <- integer(initial_count)
  polya_starts <- integer(initial_count)
  polya_ends <- integer(initial_count)
  secondary_polya_starts <- integer(initial_count)
  secondary_polya_ends <- integer(initial_count)
  pod5_files <- character(initial_count)  # New vector for pod5 files

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
    polya_length <- bam_data$tag$pt[[i]]
    pt_tag_count <- pt_tag_count + 1

    # Initialize pa tag values
    anchor_pos <- -1L
    polya_start <- -1L
    polya_end <- -1L
    secondary_polya_start <- -1L
    secondary_polya_end <- -1L

    # Add pa tag information if it exists and is valid
    if (!is.null(bam_data$tag$pa[[i]]) && length(bam_data$tag$pa[[i]]) == 5) {
      pa_values <- bam_data$tag$pa[[i]]
      if (!any(is.na(pa_values))) {
        anchor_pos <- pa_values[1]
        polya_start <- pa_values[2]
        polya_end <- pa_values[3]
        secondary_polya_start <- pa_values[4]
        secondary_polya_end <- pa_values[5]
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
    ref_starts[valid_count] <- bam_data$pos[i]
    ref_ends[valid_count] <- bam_data$pos[i] + nchar(bam_data$cigar[i])
    mapqs[valid_count] <- bam_data$mapq[i]
    polya_lengths[valid_count] <- polya_length
    anchor_positions[valid_count] <- anchor_pos
    polya_starts[valid_count] <- polya_start
    polya_ends[valid_count] <- polya_end
    secondary_polya_starts[valid_count] <- secondary_polya_start
    secondary_polya_ends[valid_count] <- secondary_polya_end
    pod5_files[valid_count] <- pod5_file
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
    polya_lengths <- polya_lengths[1:valid_count]
    anchor_positions <- anchor_positions[1:valid_count]
    polya_starts <- polya_starts[1:valid_count]
    polya_ends <- polya_ends[1:valid_count]
    secondary_polya_starts <- secondary_polya_starts[1:valid_count]
    secondary_polya_ends <- secondary_polya_ends[1:valid_count]
    pod5_files <- pod5_files[1:valid_count]

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
      polya_lengths <- polya_lengths[keep_idx]
      anchor_positions <- anchor_positions[keep_idx]
      polya_starts <- polya_starts[keep_idx]
      polya_ends <- polya_ends[keep_idx]
      secondary_polya_starts <- secondary_polya_starts[keep_idx]
      secondary_polya_ends <- secondary_polya_ends[keep_idx]
      pod5_files <- pod5_files[keep_idx]
    }

    # Create the final data frame
    polya_df <- data.frame(
      read_id = read_ids,
      pod5_file = pod5_files,  # Added pod5 file information
      reference = references,
      ref_start = ref_starts,
      ref_end = ref_ends,
      mapq = mapqs,
      polya_length = polya_lengths,
      anchor_pos = anchor_positions,
      polya_start = polya_starts,
      polya_end = polya_ends,
      secondary_polya_start = secondary_polya_starts,
      secondary_polya_end = secondary_polya_ends,
      stringsAsFactors = FALSE
    )

    cli_log("Processing complete", "SUCCESS")
    return(polya_df)
  } else {
    cli_log("No reads with poly(A) information found after filtering", "WARNING")
    return(data.frame())
  }
}
