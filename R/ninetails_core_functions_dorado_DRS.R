#' Process and split Dorado summary file into smaller parts
#'
#' Splits a Dorado summary file or data frame into multiple smaller files for downstream analysis.
#'
#' @param dorado_summary Character path to Dorado summary file, or a data frame containing summary information.
#'
#' @param save_dir Character path to directory where split summary files will be saved.
#'
#' @param part_size Integer. Number of reads per file part when splitting the summary file.
#'
#' @param cli_log Function for logging messages and progress.
#'
#' @return Character vector containing paths to the split summary files.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' summary_files <- process_dorado_summary(
#'   dorado_summary = "path/to/summary.txt",
#'   save_dir = "path/to/output/",
#'   part_size = 40000,
#'   cli_log = message
#' )
#' }
process_dorado_summary <- function(dorado_summary,
                                   save_dir,
                                   part_size,
                                   cli_log) {

  # Assertions
  if (missing(part_size)) {
    stop("Number of reads per file part (part_size) is missing. Please provide a valid part_size argument.", call. = FALSE)
  }

  assertthat::assert_that(is.numeric(part_size), part_size > 0,
                          msg = paste0("Reads per part must be numeric and positive. Please provide a valid argument."))

  # Handle input data and get prefix
  if (is.character(dorado_summary)) {
    cli_log(sprintf("Reading summary file: %s", basename(dorado_summary)), "INFO")

    # First, let's peek at the file structure
    cli_log("Checking file structure...", "INFO")
    tryCatch({
      # Read first few lines to check structure
      headers <- names(vroom::vroom(dorado_summary, n_max = 1, show_col_types = FALSE))
      cli_log(sprintf("Found columns: %s", paste(headers, collapse=", ")), "INFO")

      # Now read the actual data
      summary_data <- vroom::vroom(dorado_summary, show_col_types = FALSE)

      if (!"read_id" %in% headers) {
        cli_log("WARNING: 'read_id' column not found. Available columns are:", "WARNING")
        cli_log(paste(headers, collapse=", "), "WARNING")
        stop("Required column 'read_id' not found in summary file")
      }

    }, error = function(e) {
      cli_log(sprintf("Error reading summary file: %s", e$message), "ERROR")
      stop(e$message)
    })

    # Extract prefix from original file name
    summary_prefix <- tools::file_path_sans_ext(basename(dorado_summary))
  } else if (is.data.frame(dorado_summary)) {
    summary_data <- tibble::as_tibble(dorado_summary)
    # Check for required columns in data frame
    if (!"read_id" %in% names(summary_data)) {
      cli_log("Required column 'read_id' not found in data frame", "ERROR")
      stop("Required column 'read_id' not found in data frame")
    }
    # Use default prefix for data frame input
    summary_prefix <- "summary"
  } else {
    stop("Invalid dorado_summary format. Must be either a file path or a data frame.")
  }

  # Create output directory if it doesn't exist
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }

  # Calculate number of parts needed
  total_reads <- nrow(summary_data)
  num_parts <- ceiling(total_reads / part_size)

  cli_log(sprintf("Total reads: %d", total_reads), "INFO")
  cli_log(sprintf("Splitting into %d parts", num_parts), "INFO")

  # Initialize vector to store output file paths
  output_files <- character(num_parts)

  # Split and save summary data
  for (i in seq_len(num_parts)) {
    start_idx <- ((i - 1) * part_size) + 1
    end_idx <- min(i * part_size, total_reads)

    # Extract subset of data
    part_data <- summary_data[start_idx:end_idx, ]

    # Create output filename using prefix from original file
    output_file <- file.path(save_dir,
                             sprintf("%s_part%d.txt", summary_prefix, i))

    # Save data
    vroom::vroom_write(part_data,
                       file = output_file,
                       delim = "\t")

    output_files[i] <- output_file

    cli_log(sprintf("Saved part %d/%d: %s", i, num_parts, basename(output_file)), "INFO")
  }

  cli_log(sprintf("Summary split into %d parts", num_parts), "SUCCESS")

  return(output_files)
}



#' Split BAM file into smaller parts based on read IDs from summary file
#'
#' Splits a BAM file into multiple smaller BAM files, each containing reads matching a subset of read IDs from a summary file.
#'
#' @param bam_file Character path to the BAM file to be split.
#'
#' @param dorado_summary Character path to the Dorado summary file containing read IDs for filtering.
#'
#' @param part_size Integer. Number of reads per output file part (default: 100000).
#'
#' @param save_dir Character path to directory where split BAM files will be saved.
#'
#' @param part_number Integer. Numerical identifier for the current part being processed.
#'
#' @param cli_log Function for logging messages and progress (default: message).
#'
#' @return Character vector containing paths to the split BAM files.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' bam_files <- split_bam_file(
#'   bam_file = "path/to/aligned.bam",
#'   dorado_summary = "path/to/summary_part1.txt",
#'   save_dir = "path/to/output/",
#'   part_number = 1
#' )
#' }
split_bam_file <- function(bam_file,
                           dorado_summary,
                           part_size = 100000,
                           save_dir,
                           part_number,
                           cli_log = message) {

  # Variable binding (suppressing R CMD check from throwing an error)
  read_id <- NULL

  # Check for required Bioconductor packages
  if (!requireNamespace("Rsamtools", quietly = TRUE)) {
    stop("Package 'Rsamtools' is required for split_bam_file(). Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("S4Vectors", quietly = TRUE)) {
    stop("Package 'S4Vectors' is required for split_bam_file(). Please install it.",
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
  cli_log(sprintf("[INFO] Reading summary file: %s", basename(dorado_summary)))

  t1 <- Sys.time()
  # Get read IDs efficiently with vroom
  read_ids <- vroom::vroom(dorado_summary,
                           col_select = "read_id",
                           show_col_types = FALSE,
                           altrep = TRUE)

  total_alignments <- nrow(read_ids)
  unique_reads <- length(unique(read_ids$read_id))

  cli_log(sprintf("[INFO] Summary processed in %.2f seconds",
                  as.numeric(difftime(Sys.time(), t1, units = "secs"))))

  cli_log(sprintf("[INFO] Found %d alignments for %d unique reads",
                  total_alignments, unique_reads))

  if(total_alignments > unique_reads) {
    cli_log(sprintf("[INFO] Average alignments per read: %.2f",
                    total_alignments/unique_reads))
  }

  # Get unique read IDs and clean up
  read_ids <- unique(read_ids$read_id)
  total_reads <- length(read_ids)
  gc()

  # Calculate how many output files we'll create
  num_parts <- ceiling(total_reads / part_size)

  # Calculate number of digits needed for file numbering
  num_digits <- nchar(as.character(num_parts))

  cli_log(sprintf("[INFO] Splitting into %d files with ~%d reads each",
                  num_parts, part_size))

  # Prepare output file paths with dynamic padding
  output_files <- sprintf(paste0("%s_%0", num_digits, "d.bam"),
                          output_prefix, seq_len(num_parts))

  # Split read IDs into parts of desired size
  read_parts <- split(read_ids, ceiling(seq_along(read_ids) / part_size))

  # Set up BAM reading parameters - we only need read names for filtering
  param <- Rsamtools::ScanBamParam(what = "qname")

  # Process each part
  for (i in seq_along(read_parts)) {
    # Log progress
    cli_log(sprintf("[INFO] Processing part %d/%d (expecting %d reads)",
                    i, num_parts, length(read_parts[[i]])))

    # Prepare filter function for the current part of read IDs
    current_part <- read_parts[[i]]
    filter <- function(x) {
      matches <- x$qname %in% current_part
      return(matches)
    }

    # Filter BAM file to create part
    Rsamtools::filterBam(bam_file, output_files[i],
                         filter = S4Vectors::FilterRules(filter),
                         param = param)

    # Count actual number of reads in output file
    actual_reads <- Rsamtools::countBam(output_files[i])$records

    cli_log(sprintf("[INFO] Created: %s (contains %d alignments)",
                    basename(output_files[i]), actual_reads))

    if(actual_reads > length(current_part)) {
      cli_log(sprintf("[INFO] Note: %d additional alignments included for these reads",
                      actual_reads - length(current_part)))
    }
  }

  cli_log(sprintf("[SUCCESS] Splitting complete"))

  return(output_files)
}


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



#' Extract poly(A) tail signal segments from POD5 files
#'
#' Uses reticulate to import the Python library `pod5`. For each pod5 file referenced in
#' `polya_data`, opens a `pod5.Reader`, matches `read_id`, and extracts the numeric
#' slice `signal[poly_tail_start:poly_tail_end]` if indices are valid. Applies winsorization and interpolation.
#'
#' Only keeps reads if:
#' \itemize{
#'   \item `poly_tail_length > 10` (if the column exists), and
#'   \item `poly_tail_start > 0` (to exclude artifacts, e.g. caused by clogged pores).
#' }
#'
#' @param polya_data Data frame with columns:
#'   \describe{
#'     \item{read_id}{Character: ONT read identifiers}
#'     \item{filename}{Character: corresponding POD5 filename}
#'     \item{poly_tail_start}{Integer/numeric: start index of the poly(A) region}
#'     \item{poly_tail_end}{Integer/numeric: end index of the poly(A) region}
#'   }
#' @param pod5_dir Character. Path to the directory containing POD5 files.
#'
#' @return List of numeric vectors, each containing the extracted signal for a read.
#' @export
#'
#' @examples
#' \dontrun{
#' res <- extract_tails_from_pod5(polya_data, "/path/to/pod5_dir")
#' }
extract_tails_from_pod5 <- function(polya_data, pod5_dir) {
  i <- reader <- NULL

  # Check for required columns
  required_cols <- c("read_id", "filename", "poly_tail_start", "poly_tail_end")
  missing_cols <- setdiff(required_cols, colnames(polya_data))
  if (length(missing_cols) > 0) {
    stop(sprintf("polya_data is missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  if (!reticulate::py_module_available("pod5")) {
    stop("Python module 'pod5' is not available in the current environment.")
  }
  pod5 <- reticulate::import("pod5")

  # Keep only valid reads:
  # - poly_tail_length > 10 (if column exists)
  # - poly_tail_start > 0
  if ("poly_tail_length" %in% colnames(polya_data)) {
    polya_data <- polya_data[polya_data$poly_tail_length > 10 & polya_data$poly_tail_start > 0, ]
    if (nrow(polya_data) == 0) {
      stop("No valid reads found: require poly_tail_length > 10 and poly_tail_start > 0.")
    }
  } else {
    polya_data <- polya_data[polya_data$poly_tail_start > 0, ]
    if (nrow(polya_data) == 0) {
      stop("No valid reads found: require poly_tail_start > 0.")
    }
  }

  polya_by_file <- split(polya_data, polya_data$filename)

  signals_list <- list()
  for (filename in names(polya_by_file)) {
    current_data <- polya_by_file[[filename]]
    pod5_path <- file.path(pod5_dir, filename)
    if (!file.exists(pod5_path)) {
      warning(sprintf("POD5 file not found: %s", pod5_path))
      next
    }
    tryCatch({
      reader <- pod5$Reader(pod5_path)
      all_reads <- reticulate::iterate(reader$reads())
      available_ids <- sapply(all_reads, function(r) trimws(r$read_id))
      reads_by_id <- stats::setNames(all_reads, available_ids)
      for (i in seq_len(nrow(current_data))) {
        read_id <- trimws(current_data$read_id[i])
        if (!read_id %in% names(reads_by_id)) next
        signal <- reticulate::py_to_r(reads_by_id[[read_id]]$signal)
        start_idx <- current_data$poly_tail_start[i]
        end_idx <- current_data$poly_tail_end[i]
        if (is.numeric(start_idx) && is.numeric(end_idx) &&
            start_idx > 0 && end_idx > 0 &&
            start_idx < end_idx && end_idx <= length(signal)) {
          polya_signal <- signal[start_idx:end_idx]
          # Winsorize and interpolate
          polya_signal <- ninetails::winsorize_signal(polya_signal)
          polya_signal <- round(stats::approx(polya_signal, method = "linear", n = ceiling(0.2 * length(polya_signal)))[[2]], digits = 0)
        } else {
          polya_signal <- numeric(0)
        }
        signals_list[[read_id]] <- as.numeric(polya_signal)
      }
      reader$close()
    }, error = function(e) {
      warning(sprintf("Error processing %s: %s", pod5_path, e$message))
    })
  }
  return(signals_list)
}


#' Preprocess Dorado input files for poly(A) analysis
#'
#' Validates and splits input files, extracts poly(A) information and signals, and prepares all files for downstream analysis.
#'
#' @param bam_file Character. Path to the BAM file containing aligned reads.
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
#'     \item{bam_files}{Paths to split BAM files (if BAM processing performed)}
#'     \item{polya_files}{Paths to extracted poly(A) data files}
#'     \item{polya_signal_files}{Paths to extracted poly(A) signal files}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' processed_files <- preprocess_inputs(
#'   bam_file = "path/to/aligned.bam",
#'   dorado_summary = "path/to/summary.txt",
#'   pod5_dir = "path/to/pod5/",
#'   num_cores = 2,
#'   qc = TRUE,
#'   save_dir = "path/to/output/",
#'   prefix = "experiment1",
#'   part_size = 40000,
#'   cli_log = message
#' )
#' }
preprocess_inputs <- function(bam_file,
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
  assertthat::assert_that(is.numeric(num_cores), num_cores > 0,
                          msg = "Number of cores must be a positive numeric value")
  assertthat::assert_that(is.character(pod5_dir), dir.exists(pod5_dir),
                          msg = "Pod5 files directory does not exist or path is invalid")
  assertthat::assert_that(is.character(save_dir),
                          msg = "Output directory path must be a character string")
  assertthat::assert_that(is.logical(qc),
                          msg = "qc must be logical [TRUE/FALSE]")
  if (nchar(prefix) > 0) {
    assertthat::assert_that(is.character(prefix),
                            msg = "File name prefix must be a character string")
  }
  if (checkmate::test_string(bam_file)) {
    assertthat::assert_that(file.exists(bam_file),
                            msg = "BAM file does not exist")
  }
  if (checkmate::test_string(dorado_summary)) {
    assertthat::assert_that(file.exists(dorado_summary),
                            msg = "Dorado summary file does not exist")
  } else {
    assertthat::assert_that(is.data.frame(dorado_summary) && nrow(dorado_summary) > 0,
                            msg = "Dorado summary must be a non-empty data frame or valid file path")
  }
  cli_log("Provided input files/paths are in correct format", "SUCCESS")

  ################################################################################
  # CHECK FOR POLYA COLUMNS IN SUMMARY
  ################################################################################
  if (is.character(dorado_summary)) {
    summary_preview <- vroom::vroom(dorado_summary, n_max = 1, show_col_types = FALSE)
    summary_cols <- colnames(summary_preview)
  } else if (is.data.frame(dorado_summary)) {
    summary_cols <- colnames(dorado_summary)
  } else {
    stop("dorado_summary must be a file path or a data frame")
  }
  has_polya_cols <- all(c("poly_tail_length", "poly_tail_start", "poly_tail_end") %in% summary_cols)

  ################################################################################
  # DORADO SUMMARY PROCESSING
  ################################################################################
  cli_log("Processing dorado summary...", "INFO", "Dorado Summary Processing", bullet = TRUE)
  summary_dir <- file.path(save_dir, "dorado_summary_dir")
  if (!dir.exists(summary_dir)) {
    dir.create(summary_dir, recursive = TRUE)
  }
  part_files <- ninetails::process_dorado_summary(
    dorado_summary = dorado_summary,
    save_dir = summary_dir,
    part_size = part_size,
    cli_log = cli_log
  )
  if (length(part_files) == 0) {
    stop("No summary parts were created")
  }

  ################################################################################
  # ROUTE SELECTION: DIRECT SUMMARY OR BAM PROCESSING
  ################################################################################
  if (has_polya_cols) {
    cli_log("Summary contains poly(A) columns. Skipping BAM processing.", "INFO", "Route Selection", bullet = TRUE)
    bam_files <- NULL
    polya_files <- part_files
  } else {
    cli_log("Summary does not contain poly(A) columns. Performing BAM processing.", "INFO", "Route Selection", bullet = TRUE)
    ################################################################################
    # BAM PRE-PROCESSING
    ################################################################################
    bam_dir <- file.path(save_dir, "bam_dir")
    if (!dir.exists(bam_dir)) {
      dir.create(bam_dir, recursive = TRUE)
    }
    cli_log("Processing BAM file...", "INFO", "BAM Processing", bullet = TRUE)
    if (length(part_files) > 1) {
      cli_log(sprintf("Splitting BAM file to match %d summary parts", length(part_files)), "INFO")
      bam_files <- vector("character", length(part_files))
      for (i in seq_along(part_files)) {
        cli_log(sprintf("Processing BAM part %d/%d", i, length(part_files)), "INFO")
        current_summary <- part_files[i]
        if (!file.exists(current_summary)) {
          stop(sprintf("Summary file does not exist: %s", current_summary))
        }
        bam_files[i] <- ninetails::split_bam_file(
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
      new_bam_path <- file.path(bam_dir, basename(bam_file))
      file.copy(from = bam_file, to = new_bam_path)
      bam_files <- new_bam_path
      cli_log("BAM file copied as single file", "SUCCESS")
    }
    ################################################################################
    # EXTRACTING POLYA FROM BAM
    ################################################################################
    cli_log("Extracting poly(A) info from BAM file...", "INFO", "Poly(A) info Extracting", bullet = TRUE)
    polya_dir <- file.path(save_dir, "polya_data_dir")
    if (!dir.exists(polya_dir)) {
      dir.create(polya_dir, recursive = TRUE)
    }
    polya_files <- character(length(bam_files))
    for (i in seq_along(bam_files)) {
      current_bam <- bam_files[i]
      current_summary <- part_files[i]
      cli_log(sprintf("Processing BAM file %d/%d for poly(A) extraction", i, length(bam_files)), "INFO")
      polya_data <- ninetails::extract_polya_from_bam(
        bam_file = current_bam,
        summary_file = current_summary,
        cli_log = cli_log
      )
      if (nrow(polya_data) > 0) {
        output_file <- file.path(polya_dir,
                                 sprintf("%spolya_data_part%d.tsv",
                                         ifelse(nchar(prefix) > 0, paste0(prefix, "_"), ""),
                                         i))
        vroom::vroom_write(polya_data, output_file, delim = "\t")
        polya_files[i] <- output_file
        cli_log(sprintf("Saved poly(A) data for BAM %d to %s (%d reads)",
                        i, output_file, nrow(polya_data)), "INFO")
      } else {
        cli_log(sprintf("No poly(A) data extracted from BAM %d", i), "INFO")
      }
      gc()
    }
    polya_files <- polya_files[file.exists(polya_files)]
    if (length(polya_files) > 0) {
      cli_log(sprintf("Successfully processed %d BAM files", length(polya_files)), "SUCCESS")
    } else {
      cli_log("No poly(A) data was extracted from any BAM file", "WARNING")
    }
  }

  ################################################################################
  # EXTRACTING POLYA SIGNALS FROM POD5
  ################################################################################
  cli_log("Signal extracting", "INFO", "Signal Extracting")
  cli_log("Extracting poly(A) signals from pod5 file...", "INFO", bullet = TRUE)
  polya_signal_dir <- file.path(save_dir, "polya_signal_dir")
  if (!dir.exists(polya_signal_dir)) {
    dir.create(polya_signal_dir, recursive = TRUE)
  }
  polya_signal_files <- character(length(polya_files))
  for (i in seq_along(polya_files)) {
    polya_data <- vroom::vroom(polya_files[i], delim = "\t", show_col_types = FALSE)
    signal_list <- ninetails::extract_tails_from_pod5(polya_data, pod5_dir)
    output_file <- file.path(polya_signal_dir,
                             sprintf("%spolya_signal_part%d.rds",
                                     ifelse(nchar(prefix) > 0, paste0(prefix, "_"), ""),
                                     i))
    saveRDS(signal_list, output_file)
    polya_signal_files[i] <- output_file
    cli_log(sprintf("Saved extracted signals for poly(A) data part %d to %s", i, output_file), "INFO")
  }

  ################################################################################
  # RETURN OUTPUT FILE PATHS
  ################################################################################
  return(list(
    summary_files = part_files,
    bam_files = bam_files,
    polya_files = polya_files,
    polya_signal_files = polya_signal_files
  ))
}



#' Creates a nested list of Dorado tail features (raw signal + pseudomoves).
#'
#' This function processes raw poly(A) tail signal traces from Dorado and computes
#' pseudomoves using a threshold-based filter. Each read is stored in a nested
#' list structure, keyed by its read ID, with both the raw signal and the
#' computed pseudomove vector.
#'
#' Unlike Guppy-based pipelines, Dorado does not provide moves directly in BAM
#' files, so pseudomoves are computed empirically from the raw signal traces.
#'
#' Parallel execution is handled with \pkg{foreach} and \pkg{doSNOW}, and a
#' progress bar is displayed during processing.
#'
#' @param signal_list list of numeric vectors. Each element must represent
#' a raw poly(A) tail signal trace, with list names corresponding to read IDs.
#' @param num_cores numeric [1]. Number of physical cores to use in processing.
#' Do not exceed 1 less than the total available cores on your system.
#'
#' @return A nested list of tail features, organized by read IDs. Each read entry
#' contains:
#'   \itemize{
#'     \item \code{tail_signal}: numeric vector of the raw poly(A) tail signal.
#'     \item \code{tail_pseudomoves}: numeric vector of pseudomove states
#'           computed by \code{ninetails::filter_signal_by_threshold}.
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' features <- ninetails::create_tail_features_list_dorado(
#'   signal_list = signals,
#'   num_cores = 4
#' )
#' }
create_tail_features_list_dorado <- function(signal_list,
                                             num_cores) {

  # variable binding (suppressing R CMD check from throwing an error)
  x <- NULL

  # Assertions
  if (missing(signal_list)) {
    stop("Signal list is missing. Please provide a valid signal_list argument.", call. = FALSE)
  }
  if (missing(num_cores)) {
    stop("Number of declared cores is missing. Please provide a valid num_cores argument.", call. = FALSE)
  }

  assertthat::assert_that(is.list(signal_list),
                          msg = "Provided signal_list is not a list. Please provide a valid list of raw signal vectors.")
  assertthat::assert_that(is.numeric(num_cores),
                          msg = "Declared num_cores must be numeric. Please provide a valid value.")
  assertthat::assert_that(all(sapply(signal_list, is.numeric)),
                          msg = "All elements of signal_list must be numeric vectors (raw tail signal traces).")


  # creating cluster for parallel computing
  my_cluster <- parallel::makeCluster(num_cores)
  on.exit(parallel::stopCluster(my_cluster))
  doSNOW::registerDoSNOW(my_cluster)
  `%dopar%` <- foreach::`%dopar%`
  `%do%` <- foreach::`%do%`
  # Options for multicore execution
  mc_options <- list(preschedule = TRUE, set.seed = FALSE, cleanup = FALSE)

  #creating list for outputs

  output_nested_list <- list()

  #progressbar header
  cat(paste0('[', as.character(Sys.time()), '] ','Computing pseudomoves...', '\n', sep=''))

  # progress bar
  pb <- utils::txtProgressBar(min = 0,
                              max = length(signal_list),
                              style = 3,
                              width = 50,
                              char = "=",
                              file = stderr())
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  output_nested_list <- list()

  # Main parallel loop
  # For each read (named element of signal_list):
  #   - store the raw signal trace
  #   - compute pseudomoves with ninetails::filter_signal_by_threshold()
  #   - return as a named sublist under the read ID
  output_nested_list <- foreach::foreach(x = names(signal_list),
                                         .combine = c,
                                         .inorder = TRUE,
                                         .errorhandling = "pass",
                                         .options.snow = opts,
                                         .options.multicore = mc_options) %dopar% (function(x) {
                                           stats::setNames(list(list(
                                             tail_signal = signal_list[[x]],
                                             tail_pseudomoves = ninetails::filter_signal_by_threshold(signal_list[[x]]))
                                           ), x)
                                         })(x)

  #close(pb)
  # Done comm
  cat(paste0('[', as.character(Sys.time()), '] ','Done!', '\n', sep=''))

  return(output_nested_list)
}


#' Extracts fragments of poly(A) tail signal (Dorado mode) containing potential modifications
#' along with their delimitation (positional indices; coordinates) within the tail.
#'
#' This function processes raw poly(A) tail signal and Dorado-derived pseudomoves
#' to identify and extract signal segments (chunks) potentially corresponding
#' to modified positions (e.g., non-A residues). Each extracted chunk spans
#' 100 signal points, centered on the midpoint of a pseudomove run.
#'
#' In the Dorado pipeline, moves are not used:
#'   * retrieving them from BAM files is computationally expensive
#'   * processing is non-intuitive
#'
#' Instead, only pseudomoves are considered. As a safeguard against Dorado’s
#' tendency to extend poly(A) boundaries into the transcript body, the last
#' 3 pseudomove values are forced to 0. This prevents misclassification of
#' transcript nucleotides as part of the tail.
#'
#' Candidate modification regions are detected by:
#'   * run-length encoding (RLE) of the pseudomove vector
#'   * filtering runs of pseudomoves with length ≥ 5
#'
#' Extracted fragments are padded/imputed if they extend beyond signal boundaries:
#'   * upstream/downstream missing values (NAs) are replaced
#'   * imputation is based on random draws from the 5 most frequent signal values
#'
#' The function returns a list object (nested), where each element represents
#' one candidate modification region, containing:
#'   * `chunk_sequence`: the raw signal subsequence (length = 100, imputed if needed)
#'   * `chunk_start_pos`: starting index of the subsequence
#'   * `chunk_end_pos`: ending index of the subsequence
#'
#' @param readname character string. Name of the given read within the
#' analyzed dataset.
#'
#' @param tail_feature_list list object produced by
#'   \code{\link{create_tail_features_list_dorado}} Dorado-tail feature
#'   extraction function. Must contain
#'   \code{$tail_signal} (numeric vector) and \code{$tail_pseudomoves} (integer vector)
#'   for each read.
#'
#' @return a nested list where each element corresponds to a signal fragment.
#' Each fragment is itself a list with three entries:
#'   \itemize{
#'     \item \code{chunk_sequence}: numeric vector of raw signal values
#'     \item \code{chunk_start_pos}: integer, start index of the chunk
#'     \item \code{chunk_end_pos}: integer, end index of the chunk
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ninetails::split_tail_centered_dorado(
#'   readname = "1234-anexample-r3adn4m3",
#'   tail_feature_list = tail_feature_list
#' )
#' }
split_tail_centered_dorado <- function(readname, tail_feature_list) {

  # variable binding (suppressing R CMD check from throwing an error)
  i <- NULL

  #assertions
  if (missing(readname)) {
    stop("Readname is missing. Please provide a valid readname argument.", call. =FALSE)
  }

  if (missing(tail_feature_list)) {
    stop("List of tail features is missing. Please provide a valid tail_feature_list argument.", call. =FALSE)
  }

  assertthat::assert_that(is.character(readname),
                          msg = paste0("Given readname is not a character string. Please provide a valid readname."))
  assertthat::assert_that(is.list(tail_feature_list),
                          msg = paste0("Given tail_feature_list is not a list (class). Please provide valid file format."))

  # Extract required data
  # In the Dorado pipeline, moves are not used:
  #   - retrieving moves from BAM is time-consuming
  #   - processing moves is non-intuitive
  signal <- tail_feature_list[[readname]]$tail_signal
  pseudomoves <- tail_feature_list[[readname]]$tail_pseudomoves

  # Force the last 3 pseudomove values to 0 if they are not already 0
  # Dorado often extends tail boundaries slightly into the transcript body,
  # which can cause false positives by misclassifying transcript nucleotides
  # as part of the poly(A) tail.
  # Instead of adjusting tail boundaries directly, we zero these trailing positions.
  # This strategy is chosen with the expectation that ONT will eventually refine
  # their poly(A) detection algorithm.
  if (length(pseudomoves) >= 3) {
    idx <- (length(pseudomoves)-2):length(pseudomoves)
    pseudomoves[idx][pseudomoves[idx] != 0] <- 0
  }

  # mod-centered chunk extraction
  # Strategy:
  # - Compute run-length encoding (RLE) of pseudomoves
  # - Select runs that:
  #     * are "positive" pseudomove stretches (value != 0)
  #     * are at least 5 bases long
  # - These runs are likely to correspond to poly(A) tail signal
  mod_rle <- rle(pseudomoves)
  # pseudomoves filtered by condition (potentially decorated - empyrical!)
  condition <- mod_rle$lengths >= 5 & mod_rle$values != 0
  if (!any(condition)) return(NULL)  # No valid chunks, return NULL
  # For each valid run:
  # - Compute its starting position
  # - Compute chunk boundaries centered on the middle of the run
  # - Each chunk has a fixed window size of 100 bases (±50 around center)
  first_filtered_positions <- cumsum(c(1, utils::head(mod_rle$lengths, -1)))[condition]
  # length of pseudomoves satisfying condition
  filtered_length <- mod_rle$lengths[condition]
  # extracted coordinates (indices)
  start_positions <- first_filtered_positions + floor(filtered_length / 2) - 50
  end_positions <- start_positions + 99

  # extract signal chunks centered on potential modification
  # - If the start position is < 1, left-pad with NAs
  # - Later, NAs are imputed with frequent signal values
  list_1 <- lapply(seq_along(start_positions), function(i) {
    if (start_positions[i] > 0) signal[start_positions[i]:end_positions[i]]
    else c(rep(NA, abs(start_positions[i] - 1)), signal[1:end_positions[i]])
  })

  # NAs introduced by left-padding are replaced by random draws
  # from the 5 most frequent observed signal values.
  # This ensures consistent chunk length without introducing
  # artificial bias from fixed-value padding.
  most_freq_vals <- as.numeric(names(sort(table(signal), decreasing = TRUE)[1:5]))
  list_1 <- lapply(list_1, function(n) replace(n, is.na(n), sample(most_freq_vals, sum(is.na(n)), TRUE)))

  # naming chunks
  #   <readname>_<chunk_index>
  chunk_names <- paste0(rep(readname, length(list_1)), '_', seq_along(start_positions))
  names(list_1) <- chunk_names

  # retrieve coordinates as list_2 and list_3:
  list_2 <- as.list(start_positions)
  list_3 <- as.list(end_positions)

  # Output is a list of named lists:
  #   - chunk_sequence   (raw signal segment, length = 100)
  #   - chunk_start_pos  (starting index in original signal)
  #   - chunk_end_pos    (ending index in original signal)
  out <- purrr::transpose(purrr::set_names(list(list_1, list_2, list_3), c('chunk_sequence', 'chunk_start_pos', 'chunk_end_pos')))
  return(out)
}


#' Creates list of poly(A) tail chunks (Dorado mode) centered on significant signal deviations.
#'
#' Processes raw poly(A) tail signals and pseudomoves (as generated by Dorado)
#' in parallel to extract candidate signal fragments potentially containing
#' non-A nucleotides. Each fragment is 100 signal points long, centered on
#' pseudomove runs of sufficient length, and is returned with its positional
#' coordinates. The resulting data are organized into a nested list keyed
#' by read IDs.
#'
#' This Dorado-specific function differs from the Guppy-based version:
#'   * moves are not used (to avoid costly BAM parsing and processing)
#'   * pseudomoves are corrected at the tail ends (last 3 values forced to 0)
#'
#' Parallelization is handled with \pkg{foreach} and \pkg{doSNOW}, allowing
#' efficient scaling across multiple CPU cores. A progress bar is displayed
#' to monitor job completion.
#'
#' @param tail_feature_list list object produced by \code{\link{create_tail_feature_list}}
#' or an equivalent Dorado-tail feature extraction function. Must contain
#' per-read entries with \code{$tail_signal} and \code{$tail_pseudomoves}.
#'
#' @param num_cores numeric [1]. Number of physical cores to use in processing.
#' Do not exceed 1 less than the number of available cores on your machine.
#'
#' @return A nested list containing the segmented tail data (chunks and coordinates),
#' organized by read IDs. Each read entry contains one or more fragments, where
#' each fragment is a list with:
#'   \itemize{
#'     \item \code{chunk_sequence}: numeric vector of raw signal values (length 100)
#'     \item \code{chunk_start_pos}: integer, starting index of the chunk
#'     \item \code{chunk_end_pos}: integer, ending index of the chunk
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' tcl_dorado <- ninetails::create_tail_chunk_list_dorado(
#'   tail_feature_list = tfl,
#'   num_cores = 3
#' )
#' }
create_tail_chunk_list_dorado <- function(tail_feature_list, num_cores) {
  # variable binding (suppressing R CMD check from throwing an error)
  i <- NULL

  # initial assertions
  if (missing(num_cores)) {
    stop("Number of declared cores is missing. Please provide a valid num_cores argument.", call. =FALSE)
  }

  if (missing(tail_feature_list)) {
    stop("List of features is missing. Please provide a valid tail_feature_list argument.", call. =FALSE)
  }

  assertthat::assert_that(is.numeric(num_cores),
                          msg = paste0("Declared core number must be numeric. Please provide a valid argument."))
  assertthat::assert_that(is.list(tail_feature_list),
                          msg = paste0("Given tail_feature_list is not a list (class). Please provide valid file format."))

  # creating cluster for parallel computing
  my_cluster <- parallel::makeCluster(num_cores)
  on.exit(parallel::stopCluster(my_cluster))
  doSNOW::registerDoSNOW(my_cluster)
  `%dopar%` <- foreach::`%dopar%`
  mc_options <- list(preschedule = TRUE, set.seed = FALSE, cleanup = FALSE)

  # header for progress bar
  cat(paste0('[', as.character(Sys.time()), '] ','Creating tail segmentation data...', '\n'))

  # progress bar
  pb <- utils::txtProgressBar(min = 0,
                              max = length(tail_feature_list),
                              style = 3,
                              width = 50,
                              char = "=",
                              file = stderr())
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  #output list
  tail_chunk_list <- list()

  # For each read (element of tail_feature_list):
  #   - Call split_tail_centered_dorado() on the read
  #   - Collect resulting signal chunks into a nested list
  tail_chunk_list <- foreach::foreach(i = seq_along(tail_feature_list),
                                      .combine = c,
                                      .inorder = TRUE,
                                      .errorhandling = 'pass',
                                      .options.snow = opts,
                                      .options.multicore = mc_options) %dopar% {
                                        lapply(names(tail_feature_list[i]), function(x) ninetails::split_tail_centered_dorado(x, tail_feature_list))
                                      }

  names(tail_chunk_list) <- names(tail_feature_list)
  # Remove NULLs (reads with no valid chunks)
  tail_chunk_list <- rrapply::rrapply(tail_chunk_list,
                                      condition = Negate(is.null),
                                      how = "prune")
  # Done comment
  cat(paste0('[', as.character(Sys.time()), '] ','Done!', '\n'))

  return(tail_chunk_list)
}


