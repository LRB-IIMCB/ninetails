################################################################################
# FUNCTIONS DEVELOPED TO PRODUCE TRAINING SETS FROM DORADO cDNA DATA
################################################################################
# The cDNA training-set branch reuses as much of the Guppy training-set
# machinery as possible. Only the signal source (POD5 + BAM instead of
# multi-fast5 + nanopolish) and the orientation split (polyA vs polyT)
# are new. The feature list produced here keeps the four-slot per-read
# layout of extract_tail_data_trainingset() so that
# split_tail_centered_trainingset(), split_with_overlaps(),
# create_tail_chunk_list_trainingset(), create_tail_chunk_list_A(),
# filter_nonA_chunks_trainingset(), create_gaf_list() and
# create_gaf_list_A() can be called unchanged.
#
# Orientation conventions (Dorado cDNA, no signal reversal is applied):
#   polyA reads: transcript body | poly(A) tail | VNP (rc); the tail signal
#                runs from the transcript body towards the 3' end, opposite
#                to the DRS convention (3' end first).
#   polyT reads: VNP | poly(T) tail | transcript body (rc); the tail signal
#                runs from the 3' end of the original tail towards the body,
#                and every residue is observed as its complement (C -> G,
#                G -> C, U -> A, A -> T).
# Because the two orientations are trained as separate models, the signals
# are used as extracted. Class labels always refer to the residue in the
# original RNA tail (known from the experimental design), not to the
# complementary base read in polyT orientation.
################################################################################

#' Extracts poly(A)/poly(T) tail signals of Dorado cDNA reads for
#' training-set preparation, split by read orientation.
#'
#' This is the cDNA counterpart of the fast5/nanopolish-based read
#' extraction performed inside \code{\link{create_tail_feature_list_trainingset}}.
#' It filters the Dorado summary with \code{\link{filter_dorado_summary}},
#' optionally restricts the reads to the contigs of interest (the spike-in
#' constructs whose tail composition is known from the experimental design),
#' extracts basecalled sequences from the BAM file with
#' \code{\link{extract_data_from_bam}}, classifies the read orientation with
#' \code{\link{detect_orientation_single}}, and extracts the winsorized and
#' downsampled tail signals from POD5 files with
#' \code{\link{extract_tails_from_pod5}}.
#'
#' @details
#' The signal extraction is delegated to the same Python helper as the
#' production cDNA pipeline, so the resulting signals are winsorized
#' (0.5\% and 99.5\% percentiles) and interpolated to 20\% of their
#' original length, exactly as in \code{\link{extract_tail_data_trainingset}}.
#' Reads whose signal could not be extracted (empty vectors) are dropped.
#'
#' No signal reversal is applied. In polyA orientation the tail signal
#' runs from the transcript body towards the 3' end, in polyT orientation
#' it runs from the 3' end (right after the VNP primer) towards the body
#' and every residue is observed as its complement. Keep this in mind when
#' inspecting the signals; the downstream models are orientation-specific.
#'
#' @param dorado_summary Character string or data frame. Full path of the
#'   Dorado summary file (\code{dorado summary} run on the aligned BAM with
#'   poly(A) estimation enabled) or the already loaded data frame. Must
#'   contain \code{read_id}, \code{filename} (or \code{input_filename}),
#'   \code{poly_tail_length}, \code{poly_tail_start}, \code{poly_tail_end},
#'   \code{alignment_genome}, \code{alignment_direction} and
#'   \code{alignment_mapq}.
#'
#' @param bam_file Character string. Full path of the aligned BAM file with
#'   basecalled sequences and Dorado \code{pt}/\code{pa} tags.
#'
#' @param pod5_dir Character string. Full path of the directory containing
#'   POD5 files referenced in the summary.
#'
#' @param num_cores Numeric \code{[1]}. Number of physical cores to use for
#'   the POD5 extraction. Do not exceed 1 less than the number of cores at
#'   your disposal.
#'
#' @param contig Character vector \code{[NA]}. Names of the reference
#'   contigs (\code{alignment_genome}) to retain. Use it to select the reads
#'   mapping to the construct carrying the residue of interest. If
#'   \code{NA}, all mapped reads are retained.
#'
#' @return A named list with three elements:
#' \describe{
#'   \item{polya_signals}{Named list of numeric vectors. Tail signals of
#'     reads classified as polyA, named by read ID.}
#'   \item{polyt_signals}{Named list of numeric vectors. Tail signals of
#'     reads classified as polyT, named by read ID.}
#'   \item{read_annotation}{Data frame with one row per extracted read:
#'     \code{read_id}, \code{contig}, \code{tail_type} (polyA, polyT or
#'     unidentified), \code{poly_tail_length}, \code{poly_tail_start},
#'     \code{poly_tail_end} and \code{signal_length} (data points after
#'     downsampling). Use it as a lookup to verify that the reads come
#'     from the expected construct and orientation.}
#' }
#' Reads classified as unidentified are listed in \code{read_annotation}
#' but their signals are not returned. Always assign this returned list to
#' a variable; printing the full list to the console may crash the R session.
#'
#' @seealso \code{\link{create_tail_feature_list_trainingset_cdna}} for the
#'   next pipeline step,
#'   \code{\link{extract_tails_from_pod5}} for the signal extraction,
#'   \code{\link{detect_orientation_single}} for the orientation call,
#'   \code{\link{prepare_trainingset_cdna}} for the top-level wrapper.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' extracted <- ninetails::extract_tail_signals_trainingset_cdna(
#'   dorado_summary = '/path/to/dorado_summary.txt',
#'   bam_file = '/path/to/aligned.bam',
#'   pod5_dir = '/path/to/pod5_dir/',
#'   num_cores = 10,
#'   contig = "RlucB")
#'
#' }
#'
extract_tail_signals_trainingset_cdna <- function(dorado_summary,
                                                  bam_file,
                                                  pod5_dir,
                                                  num_cores = 1,
                                                  contig = NA) {

  #assertions
  if (missing(dorado_summary)) {
    stop(
      "Dorado summary is missing. Please provide a valid dorado_summary argument.",
      call. = FALSE
    )
  }

  if (missing(bam_file)) {
    stop(
      "BAM file is missing. Please provide a valid bam_file argument.",
      call. = FALSE
    )
  }

  if (missing(pod5_dir)) {
    stop(
      "Directory with POD5 files is missing. Please provide a valid pod5_dir argument.",
      call. = FALSE
    )
  }

  assert_condition(
    is.numeric(num_cores),
    "Declared core number must be numeric. Please provide a valid argument."
  )
  assert_condition(
    is_string(bam_file),
    "Path to BAM file is not a character string. Please provide a valid bam_file argument."
  )
  assert_file_exists(bam_file, "BAM")
  assert_condition(
    is_string(pod5_dir),
    "Path to POD5 files is not a character string. Please provide a valid pod5_dir argument."
  )
  assert_dir_exists(pod5_dir, "POD5")
  assert_condition(
    all(is.na(contig)) || is.character(contig),
    "Contig must be a character vector or NA. Please provide a valid contig argument."
  )

  # read summary (path or in-memory data frame)
  if (is_string(dorado_summary)) {
    assert_file_exists(dorado_summary, "Dorado summary")
    dorado_summary <- vroom::vroom(dorado_summary, show_col_types = FALSE)
  } else if (!is.data.frame(dorado_summary) || nrow(dorado_summary) == 0) {
    stop(
      "Empty data frame provided as an input (dorado_summary). Please provide valid input"
    )
  }

  # From dorado 1.4.0, column name for pod5 files is "input_filename"
  if (!"filename" %in% colnames(dorado_summary) && "input_filename" %in% colnames(dorado_summary)) {
    colnames(dorado_summary)[colnames(dorado_summary) == "input_filename"] <- "filename"
  }

  required_cols <- c(
    "read_id",
    "filename",
    "poly_tail_length",
    "poly_tail_start",
    "poly_tail_end",
    "alignment_genome",
    "alignment_direction",
    "alignment_mapq"
  )
  missing_cols <- setdiff(required_cols, colnames(dorado_summary))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Required columns missing: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  # quality filtering (mapped, mapq > 0, valid coordinates, tail >= 10 nt)
  # same criteria as in the production pipelines
  dorado_summary <- ninetails::filter_dorado_summary(dorado_summary)

  # keep only reads mapping to the construct(s) of interest
  if (!all(is.na(contig))) {
    dorado_summary <- dorado_summary[
      dorado_summary$alignment_genome %in% contig,
    ]
  }

  if (nrow(dorado_summary) == 0) {
    stop(
      "No reads left after quality/contig filtering. Please check the dorado_summary and contig arguments."
    )
  }

  # extract basecalled sequences from BAM
  # extract_data_from_bam() requires the summary as a file, so the filtered
  # summary is written to a temporary tsv (only filename and read_id are read)
  cat(paste0(
    '[',
    as.character(Sys.time()),
    '] ',
    'Extracting basecalled sequences from BAM file...',
    '\n',
    sep = ''
  ))

  temp_summary <- tempfile(pattern = "trainingset_summary_", fileext = ".tsv")
  on.exit(unlink(temp_summary), add = TRUE)
  vroom::vroom_write(
    dorado_summary[, c("filename", "read_id")],
    temp_summary,
    delim = "\t"
  )

  bam_data <- ninetails::extract_data_from_bam(
    bam_file = bam_file,
    summary_file = temp_summary,
    seq_only = TRUE
  )

  if (nrow(bam_data) == 0) {
    stop(
      "No reads with poly(A) tags found in the BAM file for the filtered summary. Please check the bam_file argument."
    )
  }

  # classify read orientation (Dorado-style SSP/VNP primer matching)
  cat(paste0(
    '[',
    as.character(Sys.time()),
    '] ',
    'Classifying read orientation...',
    '\n',
    sep = ''
  ))

  # sequences are coerced to character in case the BAM reader returns them
  # as a Biostrings object
  raw_orientation <- sapply(
    as.character(bam_data$sequence),
    ninetails::detect_orientation_single,
    USE.NAMES = FALSE
  )
  orientation_dict <- c("A" = "polyA", "T" = "polyT", "unknown" = "unidentified")
  bam_data$tail_type <- unname(orientation_dict[raw_orientation])

  # extract tail signals from POD5 (winsorized + downsampled to 20%)
  cat(paste0(
    '[',
    as.character(Sys.time()),
    '] ',
    'Extracting tail signals from POD5 files...',
    '\n',
    sep = ''
  ))

  signal_list <- ninetails::extract_tails_from_pod5(
    polya_data = dorado_summary,
    pod5_dir = pod5_dir,
    num_cores = num_cores
  )

  # drop reads for which extraction failed (empty vectors)
  signal_list <- Filter(function(x) length(x) > 0, signal_list)

  if (length(signal_list) == 0) {
    stop(
      "No tail signals were extracted from POD5 files. Please check the pod5_dir argument."
    )
  }

  # read annotation (lookup table for data verification)
  read_annotation <- dorado_summary[, c(
    "read_id",
    "alignment_genome",
    "poly_tail_length",
    "poly_tail_start",
    "poly_tail_end"
  )]
  names(read_annotation)[names(read_annotation) == "alignment_genome"] <- "contig"
  read_annotation <- dplyr::inner_join(
    read_annotation,
    bam_data[, c("read_id", "tail_type")],
    by = "read_id"
  )
  read_annotation <- read_annotation[
    read_annotation$read_id %in% names(signal_list),
  ]
  read_annotation$signal_length <- unname(sapply(
    signal_list[read_annotation$read_id],
    length
  ))
  #coerce tibble to df
  read_annotation <- as.data.frame(read_annotation)

  # split signals by orientation
  polya_readnames <- read_annotation$read_id[read_annotation$tail_type == "polyA"]
  polyt_readnames <- read_annotation$read_id[read_annotation$tail_type == "polyT"]

  #create final output
  extracted_signals <- list()

  extracted_signals[["polya_signals"]] <- signal_list[polya_readnames]
  extracted_signals[["polyt_signals"]] <- signal_list[polyt_readnames]
  extracted_signals[["read_annotation"]] <- read_annotation

  cat(sprintf(
    "Extracted %d reads: %d polyA, %d polyT, %d unidentified\n",
    nrow(read_annotation),
    length(polya_readnames),
    length(polyt_readnames),
    sum(read_annotation$tail_type == "unidentified")
  ))

  # Done comm
  cat(paste0('[', as.character(Sys.time()), '] ', 'Done!', '\n', sep = ''))

  return(extracted_signals)
}


#' Creates the training-set feature list (signal + pseudomoves) from Dorado
#' cDNA tail signals of a single orientation.
#'
#' This is the cDNA counterpart of
#' \code{\link{create_tail_feature_list_trainingset}} (for C, G, U) and
#' \code{\link{create_tail_feature_list_A}} (for A). It computes
#' pseudomoves for every tail signal in parallel with
#' \code{\link{filter_signal_by_threshold_trainingset}} and then applies the
#' nucleotide-specific retention criterion:
#' \itemize{
#'   \item \code{"C"}, \code{"G"}, \code{"U"}: keep reads whose pseudomove
#'         vector contains at least one non-zero run of length >= 4
#'         (potential modification present).
#'   \item \code{"A"}: keep reads whose pseudomove vector contains
#'         \emph{no} non-zero run of length >= 4 (pure homopolymer tail).
#' }
#'
#' @details
#' Dorado does not provide a move table, so the per-read feature sublists
#' carry \code{NA} in the \code{pod5_filename} and \code{tail_moves} slots.
#' The four-slot layout (\code{pod5_filename}, \code{tail_signal},
#' \code{tail_moves}, \code{tail_pseudomoves}) is preserved on purpose, so
#' that the Guppy training-set chunkers
#' (\code{\link{create_tail_chunk_list_trainingset}},
#' \code{\link{create_tail_chunk_list_A}}) and everything downstream can be
#' reused without modification. Consequently, no zero-moved read category
#' exists in this variant.
#'
#' Provide signals of a single orientation only (either
#' \code{polya_signals} or \code{polyt_signals} from
#' \code{\link{extract_tail_signals_trainingset_cdna}}); mixing orientations
#' in one feature list would produce a mixed training set.
#'
#' @param signal_list Named list of numeric vectors. Tail signals of a
#'   single orientation, named by read ID (as produced by
#'   \code{\link{extract_tail_signals_trainingset_cdna}}).
#'
#' @param num_cores Numeric \code{[1]}. Number of physical cores to use.
#'   Do not exceed 1 less than the number of cores at your disposal.
#'
#' @param nucleotide Character. One of \code{"A"}, \code{"C"}, \code{"G"} or
#'   \code{"U"}. The residue inserted into the tails of the analysed
#'   construct (known from the experimental design). Selects the retention
#'   criterion described above.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{tail_feature_list}{Named list of per-read feature lists with
#'     four slots: \code{pod5_filename} (\code{NA}), \code{tail_signal},
#'     \code{tail_moves} (\code{NA}) and \code{tail_pseudomoves}.}
#'   \item{discarded_readnames}{Character vector. Read IDs discarded by the
#'     nucleotide-specific criterion (no qualifying pseudomove run for
#'     C/G/U; a qualifying run present for A).}
#' }
#' Always assign this returned list to a variable; printing the full list
#' to the console may crash the R session.
#'
#' @seealso \code{\link{extract_tail_signals_trainingset_cdna}} for the
#'   preceding pipeline step,
#'   \code{\link{count_pseudomove_runs_trainingset_cdna}} and
#'   \code{\link{plot_tail_features_trainingset_cdna}} for data inspection,
#'   \code{\link{create_tail_chunk_list_trainingset}} and
#'   \code{\link{create_tail_chunk_list_A}} for the next pipeline step,
#'   \code{\link{prepare_trainingset_cdna}} for the top-level wrapper.
#'
#' @importFrom foreach %dopar%
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' tfl <- ninetails::create_tail_feature_list_trainingset_cdna(
#'   signal_list = extracted$polya_signals,
#'   num_cores = 10,
#'   nucleotide = "C")
#'
#' }
#'
create_tail_feature_list_trainingset_cdna <- function(signal_list,
                                                      num_cores,
                                                      nucleotide) {

  # Assertions
  if (missing(signal_list)) {
    stop(
      "Signal list is missing. Please provide a valid signal_list argument.",
      call. = FALSE
    )
  }

  if (missing(num_cores)) {
    stop(
      "Number of declared cores is missing. Please provide a valid num_cores argument.",
      call. = FALSE
    )
  }

  if (missing(nucleotide)) {
    stop(
      "Nucleotide is missing. Please provide a valid nucleotide argument.",
      call. = FALSE
    )
  }

  assert_condition(
    is.list(signal_list),
    "Provided signal_list is not a list. Please provide a valid list of tail signal vectors."
  )
  assert_condition(
    length(signal_list) > 0,
    "Provided signal_list is empty. Please provide a valid list of tail signal vectors."
  )
  assert_condition(
    !is.null(names(signal_list)) && all(nchar(names(signal_list)) > 0),
    "Elements of signal_list must be named by read IDs. Please provide a valid argument."
  )
  assert_condition(
    all(sapply(signal_list, is.numeric)),
    "All elements of signal_list must be numeric vectors (tail signal traces)."
  )
  assert_condition(
    is.numeric(num_cores),
    "Declared core number must be numeric. Please provide a valid argument."
  )

  if (!is.character(nucleotide) || length(nucleotide) != 1 || !nucleotide %in% c("A", "C", "G", "U")) {
    stop(
      "Wrong nucleotide selected. Please provide either 'A', 'C', 'G' or 'U'.",
      call. = FALSE
    )
  }

  #create empty list for extracted data
  tail_features_list = list()

  # creating cluster for parallel computing
  my_cluster <- parallel::makeCluster(num_cores)
  on.exit(parallel::stopCluster(my_cluster))
  doSNOW::registerDoSNOW(my_cluster)
  `%dopar%` <- foreach::`%dopar%`
  mc_options <- list(preschedule = TRUE, set.seed = FALSE, cleanup = FALSE)

  # header for progress bar
  cat(paste0(
    '[',
    as.character(Sys.time()),
    '] ',
    'Computing pseudomoves of provided reads...',
    '\n',
    sep = ''
  ))

  # progress bar
  pb <- utils::txtProgressBar(
    min = 0,
    max = length(signal_list),
    style = 3,
    width = 50,
    char = "=",
    file = stderr()
  )
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  # parallel extraction
  # four-slot layout kept for compatibility with the Guppy training chunkers;
  # dorado provides neither the fast5 name nor the move table, hence NAs
  tail_features_list <- foreach::foreach(
    i = seq_along(signal_list),
    .combine = c,
    .inorder = TRUE,
    .errorhandling = 'pass',
    .options.snow = opts,
    .options.multicore = mc_options
  ) %dopar%
    {
      lapply(signal_list[i], function(x) {
        extracted_data_single_list = list()
        extracted_data_single_list[["pod5_filename"]] <- NA_character_
        extracted_data_single_list[["tail_signal"]] <- x
        extracted_data_single_list[["tail_moves"]] <- NA_real_
        extracted_data_single_list[["tail_pseudomoves"]] <- ninetails::filter_signal_by_threshold_trainingset(x)
        extracted_data_single_list
      })
    }

  #label each signal according to corresponding read name to avoid confusion
  squiggle_names <- names(signal_list)
  names(tail_features_list) <- squiggle_names

  if (nucleotide == "A") {
    # prevent from running on reads which fulfill the pseudomove condition
    # (same inverse criterion as in create_tail_feature_list_A)
    find_nonzeros <- function(x, y) rowSums(stats::embed(x != 0, y))
    tail_features_list <- Filter(
      function(x) {
        length(x$tail_pseudomoves) >= 4 &&
          all(find_nonzeros(x$tail_pseudomoves, 4) < 4)
      },
      tail_features_list
    )
  } else {
    # prevent from running on reads which do not fulfill the pseudomove condition
    # (same criterion as in create_tail_feature_list_trainingset)
    tail_features_list <- Filter(
      function(x) any(with(rle(x$tail_pseudomoves), lengths[values != 0] >= 4)),
      tail_features_list
    )
  }

  # reads discarded by the nucleotide-specific criterion
  discarded_readnames <- squiggle_names[
    !(squiggle_names %in% names(tail_features_list))
  ]

  #create final output
  tail_feature_list <- list()

  tail_feature_list[["tail_feature_list"]] <- tail_features_list
  tail_feature_list[["discarded_readnames"]] <- discarded_readnames

  # Done comm
  cat(paste0('[', as.character(Sys.time()), '] ', 'Done!', '\n', sep = ''))

  return(tail_feature_list)
}


#' Counts peak and valley pseudomove runs per read in a training-set
#' feature list (lookup for the cDNA training data).
#'
#' In the Guppy/DRS training routine the polarity of the signal deviation
#' caused by a given residue was known (G produces a peak, C and U produce
#' valleys). The polarity for cDNA (DNA chemistry, two orientations, and
#' complementary bases in polyT reads) has to be established empirically
#' before \code{\link{filter_nonA_chunks_trainingset}} is called with the
#' correct \code{value}. This function tabulates, for every read, how many
#' qualifying peak (+1) and valley (-1) runs the pseudomove vector contains,
#' so that the dominant polarity of a labelled dataset can be read off the
#' column sums.
#'
#' @param tail_feature_list List object produced by
#'   \code{\link{create_tail_feature_list_trainingset_cdna}}.
#'
#' @param min_run_length Numeric \code{[4]}. Minimum length of a non-zero
#'   pseudomove run to be counted. The default matches the chunking
#'   criterion of \code{\link{split_tail_centered_trainingset}}.
#'
#' @return A data frame with one row per read and columns:
#' \describe{
#'   \item{readname}{Character. Read ID.}
#'   \item{signal_length}{Integer. Length of the downsampled tail signal.}
#'   \item{peak_runs}{Integer. Number of +1 runs of length >=
#'     \code{min_run_length}.}
#'   \item{valley_runs}{Integer. Number of -1 runs of length >=
#'     \code{min_run_length}.}
#' }
#'
#' @seealso \code{\link{create_tail_feature_list_trainingset_cdna}} for the
#'   input,
#'   \code{\link{plot_tail_features_trainingset_cdna}} for visual
#'   inspection of single reads,
#'   \code{\link{filter_nonA_chunks_trainingset}} where the polarity is
#'   used.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' run_counts <- ninetails::count_pseudomove_runs_trainingset_cdna(
#'   tail_feature_list = tfl)
#' colSums(run_counts[, c("peak_runs", "valley_runs")])
#'
#' }
#'
count_pseudomove_runs_trainingset_cdna <- function(tail_feature_list,
                                                   min_run_length = 4) {

  #assertions
  if (missing(tail_feature_list)) {
    stop(
      "List of tail features is missing. Please provide a valid tail_feature_list argument.",
      call. = FALSE
    )
  }

  assert_condition(
    is.list(tail_feature_list),
    "Given tail_feature_list is not a list (class). Please provide valid file format."
  )
  assert_condition(
    is.numeric(min_run_length) && min_run_length >= 1,
    "Minimum run length must be a positive number. Please provide a valid argument."
  )

  features <- tail_feature_list[[1]]

  if (length(features) == 0) {
    return(data.frame(
      readname = character(0),
      signal_length = integer(0),
      peak_runs = integer(0),
      valley_runs = integer(0),
      stringsAsFactors = FALSE
    ))
  } # No reads, return empty table

  run_counts <- lapply(names(features), function(x) {
    pseudomove_rle <- rle(features[[x]][["tail_pseudomoves"]])
    qualifying <- pseudomove_rle$lengths >= min_run_length
    data.frame(
      readname = x,
      signal_length = length(features[[x]][["tail_signal"]]),
      peak_runs = sum(qualifying & pseudomove_rle$values == 1),
      valley_runs = sum(qualifying & pseudomove_rle$values == -1),
      stringsAsFactors = FALSE
    )
  })
  run_counts <- do.call(rbind, run_counts)

  return(run_counts)
}


#' Draws the downsampled tail signal of a single training-set read with
#' pseudomove runs highlighted.
#'
#' Visual lookup for the cDNA training data. The shaded regions mark the
#' non-zero pseudomove runs (length >= 4) which will become the centres of
#' the extracted chunks, labelled as peak (+1) or valley (-1). Inspect a
#' handful of reads per labelled dataset and orientation to confirm that
#' the deviations look like genuine residue-induced distortions and to
#' determine their polarity.
#'
#' @param readname Character string. Name (UUID) of the given read within
#'   the feature list.
#'
#' @param tail_feature_list List object produced by
#'   \code{\link{create_tail_feature_list_trainingset_cdna}}.
#'
#' @return ggplot2 object with the tail signal and highlighted pseudomove
#'   runs.
#'
#' @seealso \code{\link{count_pseudomove_runs_trainingset_cdna}} for the
#'   tabular lookup,
#'   \code{\link{plot_tail_chunk}} for plotting extracted chunks,
#'   \code{\link{plot_gaf}} for plotting the resulting GAFs.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' example <- ninetails::plot_tail_features_trainingset_cdna(
#'   readname = "5c2386e6-32e9-4e15-a5c7-2831f4750b2b",
#'   tail_feature_list = tfl)
#'
#' print(example)
#'
#' }
#'
plot_tail_features_trainingset_cdna <- function(readname, tail_feature_list) {

  #assertions
  if (missing(readname)) {
    stop(
      "Readname is missing. Please provide a valid readname argument.",
      call. = FALSE
    )
  }

  if (missing(tail_feature_list)) {
    stop(
      "List of tail features is missing. Please provide a valid tail_feature_list argument.",
      call. = FALSE
    )
  }

  assert_condition(
    is.character(readname),
    "Given readname is not a character string. Please provide a valid readname."
  )
  assert_condition(
    is.list(tail_feature_list),
    "Given tail_feature_list is not a list (class). Please provide valid file format."
  )
  assert_condition(
    readname %in% names(tail_feature_list[[1]]),
    "Given readname is not present in the tail_feature_list. Please provide a valid readname."
  )

  #extract required data
  signal <- tail_feature_list[[1]][[readname]][["tail_signal"]]
  pseudomoves <- tail_feature_list[[1]][[readname]][["tail_pseudomoves"]]

  #create signal dataframe for plotting
  signal_df <- data.frame(
    position = seq_along(signal),
    signal = signal
  )

  # locate qualifying pseudomove runs (same criterion as chunk extraction)
  pseudomove_rle <- rle(pseudomoves)
  run_ends <- cumsum(pseudomove_rle$lengths)
  run_starts <- run_ends - pseudomove_rle$lengths + 1
  qualifying <- pseudomove_rle$values != 0 & pseudomove_rle$lengths >= 4

  run_palette <- c(
    "peak"   = "#50a675",
    "valley" = "#3a424f"
  )

  # highlight runs behind the signal line
  run_layers <- list()
  for (i in which(qualifying)) {
    run_type <- ifelse(pseudomove_rle$values[i] == 1, "peak", "valley")
    run_layers[[length(run_layers) + 1]] <- ggplot2::annotate(
      "rect",
      xmin = run_starts[i],
      xmax = run_ends[i],
      ymin = -Inf,
      ymax = Inf,
      fill = run_palette[[run_type]],
      alpha = .15
    )
    run_layers[[length(run_layers) + 1]] <- ggplot2::annotate(
      "text",
      x = mean(c(run_starts[i], run_ends[i])),
      y = Inf,
      label = run_type,
      vjust = 1.5,
      fontface = "bold",
      size = 3.5,
      color = run_palette[[run_type]]
    )
  }

  g.line <- ggplot2::geom_line(ggplot2::aes(y = signal), color = "#3271a8")
  g.labs <- ggplot2::labs(
    title = paste0("Read ", readname),
    subtitle = "Shaded areas: pseudomove runs of at least 4 data points (peak = +1, valley = -1)",
    x = "\nposition [data points]",
    y = "signal [raw]\n"
  )

  #plot
  plot_squiggle <- ggplot2::ggplot(
    data = signal_df,
    ggplot2::aes(x = position)
  ) +
    run_layers +
    g.line +
    g.labs +
    ggplot2::theme_bw()

  return(plot_squiggle)
}


#' Produces GAF training data of a given nucleotide and read orientation
#' from Dorado cDNA data.
#'
#' Top-level convenience wrapper mirroring \code{\link{prepare_trainingset}}
#' for the cDNA pipeline. It chains
#' \code{\link{extract_tail_signals_trainingset_cdna}},
#' \code{\link{create_tail_feature_list_trainingset_cdna}} and the reused
#' Guppy training-set chunking, filtering and GAF functions for a single
#' nucleotide and a single orientation.
#'
#' @details
#' The internal pipeline differs by nucleotide:
#' \describe{
#'   \item{\code{"A"}}{\code{\link{create_tail_chunk_list_A}}
#'     \eqn{\rightarrow} \code{\link{create_gaf_list_A}} (overlapping
#'     windows, data augmentation).}
#'   \item{\code{"C"}, \code{"G"}, \code{"U"}}{\code{\link{create_tail_chunk_list_trainingset}}
#'     \eqn{\rightarrow} \code{\link{filter_nonA_chunks_trainingset}}
#'     (with the supplied \code{value}) \eqn{\rightarrow}
#'     \code{\link{create_gaf_list}}.}
#' }
#'
#' Unlike the Guppy/DRS routine, the pseudomove polarity of each residue
#' is not hardcoded, because it has not been established for the DNA
#' chemistry and for both orientations. Determine it first with
#' \code{\link{count_pseudomove_runs_trainingset_cdna}} and
#' \code{\link{plot_tail_features_trainingset_cdna}}, then pass it as
#' \code{value}. Run the wrapper once per nucleotide and per orientation;
#' the polyA and polyT sets train two separate models.
#'
#' @param nucleotide Character. One of \code{"A"}, \code{"C"}, \code{"G"}
#'   or \code{"U"}. The residue inserted into the tails of the analysed
#'   construct (known from the experimental design).
#'
#' @param tail_type Character. Either \code{"polyA"} or \code{"polyT"}.
#'   Read orientation for which the training data are produced.
#'
#' @param dorado_summary Character string or data frame. Dorado summary
#'   (see \code{\link{extract_tail_signals_trainingset_cdna}}).
#'
#' @param bam_file Character string. Full path of the aligned BAM file.
#'
#' @param pod5_dir Character string. Full path of the directory containing
#'   POD5 files.
#'
#' @param num_cores Numeric \code{[1]}. Number of physical cores to use.
#'   Do not exceed 1 less than the number of cores at your disposal.
#'
#' @param contig Character vector \code{[NA]}. Reference contigs to retain
#'   (see \code{\link{extract_tail_signals_trainingset_cdna}}).
#'
#' @param value Numeric \code{[NA]}. Pseudomove polarity of the residue of
#'   interest in the given orientation: \code{1} to retain chunks with
#'   peaks, \code{-1} to retain chunks with valleys. Required for
#'   \code{"C"}, \code{"G"} and \code{"U"}; ignored for \code{"A"}.
#'
#' @return A named list of GAF arrays (100, 100, 2) organised by
#'   \code{<read_ID>_<index>}. Always assign this returned list to a
#'   variable; printing the full list to the console may crash the R session.
#'
#' @seealso \code{\link{prepare_trainingset}} for the Guppy counterpart,
#'   \code{\link{extract_tail_signals_trainingset_cdna}},
#'   \code{\link{create_tail_feature_list_trainingset_cdna}},
#'   \code{\link{count_pseudomove_runs_trainingset_cdna}},
#'   \code{\link{filter_nonA_chunks_trainingset}},
#'   \code{\link{create_gaf_list}},
#'   \code{\link{create_gaf_list_A}}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # C-containing construct, polyT orientation, valleys established beforehand
#' C_polyt_gafs <- ninetails::prepare_trainingset_cdna(
#'   nucleotide = "C",
#'   tail_type = "polyT",
#'   dorado_summary = '/path/to/dorado_summary.txt',
#'   bam_file = '/path/to/aligned.bam',
#'   pod5_dir = '/path/to/pod5_dir/',
#'   num_cores = 10,
#'   contig = "RlucB",
#'   value = -1)
#'
#' }
#'
prepare_trainingset_cdna <- function(nucleotide,
                                     tail_type,
                                     dorado_summary,
                                     bam_file,
                                     pod5_dir,
                                     num_cores = 1,
                                     contig = NA,
                                     value = NA) {

  #assertions
  if (missing(nucleotide)) {
    stop(
      "Nucleotide is missing. Please provide a valid nucleotide argument.",
      call. = FALSE
    )
  }

  if (missing(tail_type)) {
    stop(
      "Tail type is missing. Please provide a valid tail_type argument.",
      call. = FALSE
    )
  }

  if (!is.character(nucleotide) || length(nucleotide) != 1 || !nucleotide %in% c("A", "C", "G", "U")) {
    stop(
      "Wrong nucleotide selected. Please provide either 'A', 'C', 'G' or 'U'.",
      call. = FALSE
    )
  }

  if (!is.character(tail_type) || length(tail_type) != 1 || !tail_type %in% c("polyA", "polyT")) {
    stop(
      "Wrong tail type selected. Please provide either 'polyA' or 'polyT'.",
      call. = FALSE
    )
  }

  # polarity is required for non-A residues (not hardcoded for cDNA)
  if (nucleotide != "A") {
    assert_condition(
      is.numeric(value) && length(value) == 1 && value %in% c(-1, 1),
      "Pseudomove polarity (value) must be either 1 or -1 for C, G and U. Establish it with count_pseudomove_runs_trainingset_cdna() first."
    )
  }

  # signal extraction and orientation split
  extracted_signals <- ninetails::extract_tail_signals_trainingset_cdna(
    dorado_summary = dorado_summary,
    bam_file = bam_file,
    pod5_dir = pod5_dir,
    num_cores = num_cores,
    contig = contig
  )

  if (tail_type == "polyA") {
    signal_list <- extracted_signals[["polya_signals"]]
  } else {
    signal_list <- extracted_signals[["polyt_signals"]]
  }

  if (length(signal_list) == 0) {
    stop(
      sprintf("No %s reads found in the provided data. Please check the inputs.", tail_type),
      call. = FALSE
    )
  }

  # pseudomoves and nucleotide-specific retention
  tail_feature_list <- ninetails::create_tail_feature_list_trainingset_cdna(
    signal_list = signal_list,
    num_cores = num_cores,
    nucleotide = nucleotide
  )

  # chunking, filtering and GAF creation reuse the Guppy training-set branch
  if (nucleotide == "A") {
    # process the A-containing data; split tails with overlaps
    tail_chunk_list <- ninetails::create_tail_chunk_list_A(
      tail_feature_list,
      num_cores
    )
    gafs_list <- ninetails::create_gaf_list_A(tail_chunk_list, num_cores)
  } else {
    # process the non-A data to filter only internal, residue-containing chunks
    tail_chunk_list <- ninetails::create_tail_chunk_list_trainingset(
      tail_feature_list,
      num_cores
    )
    filtered_chunk_list <- ninetails::filter_nonA_chunks_trainingset(
      tail_chunk_list,
      value = value,
      num_cores
    )
    gafs_list <- ninetails::create_gaf_list(filtered_chunk_list, num_cores)
  }

  return(gafs_list)
}
