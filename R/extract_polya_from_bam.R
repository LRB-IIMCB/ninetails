#' Extract poly(A) tail information from BAM file
#'
#' This function extracts poly(A) tail information from a BAM file produced by
#' Dorado basecaller. The function processes BAM files containing pt (poly(A) tail length)
#' and pa (poly(A) region coordinates) tags, filtering out unmapped and non-primary
#' alignments. Only reads containing valid pt tags are included in the output.
#'
#' The BAM file must be produced with Dorado version 1.0.0 or higher
#' and contain pt and pa tags for poly(A) tail information. Files produced with earlier
#' versions of Dorado will not generate correct output due to differences in tags.
#'
#' The function applies the following filtering criteria:\itemize{
#' \item Excludes reads without pt tag
#' \item Excludes unmapped reads
#' \item Excludes non-primary alignments
#' \item Removes duplicate read names (keeps first occurrence)
#' }
#'
#' Progress and filtering statistics are reported during execution.
#'
#' @param bam_file character string. Full path to the BAM file containing Dorado
#' basecalling results. The BAM file must contain pt and pa tags for poly(A) tail
#' information.
#'
#' @return A nested list organized by read names, where each element corresponds to a read and contains:\itemize{
#' \item reference - Reference sequence name
#' \item ref_start - Start position on the reference
#' \item ref_end - End position on the reference
#' \item mapq - Mapping quality score
#' \item polya_length - Poly(A) tail length from pt tag
#' \item anchor_pos - Anchor position from pa tag
#' \item polya_start - Poly(A) start position from pa tag
#' \item polya_end - Poly(A) end position from pa tag
#' \item secondary_polya_start - Secondary poly(A) start position from pa tag
#' \item secondary_polya_end - Secondary poly(A) end position from pa tag
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' polya_data <- extract_polya_from_bam("/path/to/basecalled.bam")
#' }
#'
extract_polya_from_bam <- function(bam_file) {

  # Assertions
  if (!file.exists(bam_file)) {
    stop("BAM file does not exist: ", bam_file)}

  # Create parameter settings for scanning BAM file
  param <- Rsamtools::ScanBamParam(
    tag = c("pt", "pa"),  # tags to extract
    what = c("qname", "rname", "pos", "mapq", "flag", "cigar"))

  # messages for logging
  cat(paste0('[', as.character(Sys.time()), '] ', 'Reading BAM file: ', base::basename(bam_file), '\n', sep=''))
  bam_data <- Rsamtools::scanBam(bam_file, param = param)[[1]]
  initial_count <- length(bam_data$qname)
  cat(paste0('[', as.character(Sys.time()), '] ', sprintf('Total read count: %d', initial_count), '\n', sep=''))

  # Create output list structure
  polya_reads <- list()
  mapped_count <- 0
  primary_count <- 0
  pt_tag_count <- 0

  # Process reads directly from BAM data
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

    read_info <- list(
      reference = as.character(bam_data$rname[i]),
      ref_start = bam_data$pos[i],
      ref_end = bam_data$pos[i] + nchar(bam_data$cigar[i]),
      mapq = bam_data$mapq[i],
      polya_length = polya_length,
      anchor_pos = -1L,
      polya_start = -1L,
      polya_end = -1L,
      secondary_polya_start = -1L,
      secondary_polya_end = -1L)

    # Add pa tag information if it exists and is valid
    #
    if (!is.null(bam_data$tag$pa[[i]]) && length(bam_data$tag$pa[[i]]) == 5) {
      pa_values <- bam_data$tag$pa[[i]]
      if (!any(is.na(pa_values))) {
        read_info$anchor_pos <- pa_values[1]
        read_info$polya_start <- pa_values[2]
        read_info$polya_end <- pa_values[3]
        read_info$secondary_polya_start <- pa_values[4]
        read_info$secondary_polya_end <- pa_values[5]
      }
    }

    polya_reads[[bam_data$qname[i]]] <- read_info
  }

  # Output filtering statistics for logging
  cat(paste0('[', as.character(Sys.time()), '] ', sprintf('Mapped reads filtered: %d', mapped_count), '\n', sep=''))
  cat(paste0('[', as.character(Sys.time()), '] ', sprintf('Primary alignment filtered: %d', primary_count), '\n', sep=''))
  cat(paste0('[', as.character(Sys.time()), '] ', sprintf('Reads with pt tag filtered: %d', pt_tag_count), '\n', sep=''))

  # Check for duplicates in the final output (debuging purpose)
  if (length(polya_reads) > 0) {
    duplicate_reads <- duplicated(names(polya_reads))
    if (any(duplicate_reads)) {
      cat(paste0('[', as.character(Sys.time()), '] ',
                 sprintf('WARNING: Found %d duplicate read names, keeping first occurrence only', sum(duplicate_reads)),
                 '\n', sep=''))
      polya_reads <- polya_reads[!duplicate_reads]
    }

    cat(paste0('[', as.character(Sys.time()), '] ', 'Processing complete', '\n', sep=''))
    return(polya_reads)
  } else {
    cat(paste0('[', as.character(Sys.time()), '] ', 'No reads with poly(A) information found after filtering', '\n', sep=''))
    return(list())
  }
}
