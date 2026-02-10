################################################################################
# CORE NINETAILS FUNCTIONS SHARED ACROSS THE PIPELINES
################################################################################
#' Extract poly(A) data from nanopolish output and sequencing summary
#'
#' Extracts features of poly(A) tails of selected RNA reads from the output
#' table provided by nanopolish polya function and the sequencing summary
#' provided by the sequencer. Filenames are taken from the sequencing summary
#' file. Only reads with tail lengths estimated as >= 10 nt by nanopolish polya
#' function are taken into account.
#'
#' @param nanopolish Character string or data frame. Either the full path of
#'   the \code{.tsv} file produced by nanopolish polya function or an in-memory
#'   data frame containing nanopolish data.
#'
#' @param sequencing_summary Character string or data frame. Either the full
#'   path of the \code{.txt} file with sequencing summary or an in-memory data
#'   frame containing sequencing summary data.
#'
#' @param pass_only Logical. If \code{TRUE} (default), only reads tagged by
#'   nanopolish as \code{"PASS"} are taken into consideration. If \code{FALSE},
#'   reads tagged as \code{"PASS"} and \code{"SUFFCLIP"} are both included in
#'   the analysis.
#'
#' @details
#' The function performs the following operations:
#' \enumerate{
#'   \item Reads and validates nanopolish and sequencing summary inputs
#'         (accepts both file paths and in-memory data frames)
#'   \item Filters reads by QC tag (\code{pass_only} parameter)
#'   \item Joins nanopolish poly(A) data with sequencing summary by read name
#'   \item Filters reads with poly(A) tail length >= 10 nt
#'   \item Removes duplicate entries from secondary alignments
#' }
#'
#' @return A data frame containing read information organized by the read ID.
#'   Columns include:
#'   \describe{
#'     \item{readname}{Character. Read identifier}
#'     \item{polya_start}{Integer. Start position of the poly(A) tail in the
#'       raw signal}
#'     \item{transcript_start}{Integer. Start position of the transcript in the
#'       raw signal}
#'     \item{polya_length}{Numeric. Estimated poly(A) tail length in nucleotides}
#'     \item{qc_tag}{Character. Nanopolish quality control tag}
#'     \item{filename}{Character. Name of the source Fast5 file}
#'   }
#'   Always assign the returned data frame to a variable. Printing the full
#'   output to the console may crash your R session.
#'
#' @seealso
#' \code{\link{extract_tail_data}} for extracting tail features from individual
#' reads, \code{\link{create_tail_feature_list}} for batch feature extraction
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' ninetails::extract_polya_data(
#'   nanopolish = '/path/to/nanopolish/polya/output.tsv',
#'   sequencing_summary = '/path/to/sequencing_summary.txt',
#'   pass_only = TRUE
#' )
#'
#' }
extract_polya_data <- function(nanopolish,
                               sequencing_summary,
                               pass_only = TRUE) {

  if (missing(nanopolish)) {
    stop(
      "Nanopolish polya output is missing. Please provide a valid nanopolish argument.",
      .call = FALSE
    )
  }

  if (missing(sequencing_summary)) {
    stop(
      "Sequencing summary file is missing. Please provide a valid sequencing_summary argument.",
      .call = FALSE
    )
  }

  assert_condition(
    is.logical(pass_only),
    "Please provide TRUE/FALSE values for pass_only parameter"
  )

  # Accept either path to file or in-memory file - PK's GH issue
  if (is_string(nanopolish)) {
    # if string provided as an argument, read from file
    assert_file_exists(nanopolish, "Nanopolish")
    nanopolish_polya_table <- vroom::vroom(
      nanopolish,
      col_select = c(
        readname,
        polya_start,
        transcript_start,
        polya_length,
        qc_tag
      ),
      show_col_types = FALSE
    )
  } else {
    # make sure that nanopolish is an object with rows
    if (!is.data.frame(nanopolish) || nrow(nanopolish) == 0) {
      stop(
        "Empty data frame provided as an input (nanopolish). Please provide valid input"
      )
    }

    nanopolish_polya_table <- nanopolish[, c(
      "readname",
      "polya_start",
      "transcript_start",
      "polya_length",
      "qc_tag"
    )]
  }

  # Accept either path to file or in-memory file - PK's GH issue
  if (is_string(sequencing_summary)) {
    # if string provided as an argument, read from file
    assert_file_exists(sequencing_summary, "Sequencing summary")
    sequencing_summary_table <- vroom::vroom(
      sequencing_summary,
      col_select = c(filename, read_id),
      show_col_types = FALSE
    )
  } else {
    if (!is.data.frame(sequencing_summary) || nrow(sequencing_summary) == 0) {
      stop(
        "Empty data frame provided as an input (sequencing_summary). Please provide valid input"
      )
    }

    sequencing_summary_table <- sequencing_summary
  }

  #rename read id column
  names(sequencing_summary_table)[
    names(sequencing_summary_table) == "read_id"
  ] <- "readname"

  # Add filtering criterion: select only pass or pass $ suffclip
  if (pass_only == TRUE) {
    polya_summary <- dplyr::left_join(
      nanopolish_polya_table[which(nanopolish_polya_table$qc_tag == "PASS"), ],
      sequencing_summary_table,
      by = "readname"
    )
  } else {
    polya_summary <- dplyr::left_join(
      nanopolish_polya_table[
        which(nanopolish_polya_table$qc_tag %in% c("PASS", "SUFFCLIP")),
      ],
      sequencing_summary_table,
      by = "readname"
    )
  }

  # Add filtering criterion: tail length >= 10 nt
  polya_summary <- dplyr::filter(polya_summary, polya_length >= 10)

  # Prevent bugs from custom models used for mapping (secondary alignments)
  polya_summary <- unique(polya_summary)

  names(polya_summary$filename) <- polya_summary$readname #named vec of filenames (names = readnames)
  attr(polya_summary, 'spec') <- NULL #drop attributes left by vroom

  return(polya_summary)
}


#' Extract tail features of a single RNA read from a multi-Fast5 file
#'
#' Extracts metadata and signal features of a single RNA read from a
#' multi-Fast5 file basecalled by Guppy. The tail signal, as delimited by
#' nanopolish polya function, is extracted, winsorized (to remove signal
#' cliffs), and downsampled to 20\% of its original length to facilitate
#' further analysis. Pseudomoves are computed from the processed signal
#' using \code{\link{filter_signal_by_threshold}}.
#'
#' @param readname Character string. Name of the given read (UUID) within
#'   the analyzed dataset.
#'
#' @param polya_summary Data frame. The table containing data extracted from
#'   nanopolish and sequencing summary, as produced by
#'   \code{\link{extract_polya_data}}.
#'
#' @param workspace Character string. Full path of the directory containing
#'   the basecalled multi-Fast5 files.
#'
#' @param basecall_group Character string. Name of the level in the Fast5
#'   file hierarchy from which data should be extracted (e.g.,
#'   \code{"Basecall_1D_000"}).
#'
#' @details
#' The function performs the following operations:
#' \enumerate{
#'   \item Reads raw signal from the multi-Fast5 file via \pkg{rhdf5}
#'   \item Retrieves basecaller moves and stride information
#'   \item Extracts the poly(A) tail region of the signal based on nanopolish
#'         coordinates
#'   \item Applies winsorization to the tail signal
#'         (\code{\link{winsorize_signal}})
#'   \item Downsamples both signal and moves to 20\% of original length via
#'         linear interpolation
#'   \item Computes pseudomoves using
#'         \code{\link{filter_signal_by_threshold}}
#' }
#'
#' @return A named list containing per-read tail features:
#'   \describe{
#'     \item{fast5_filename}{Character. Name of the source Fast5 file}
#'     \item{tail_signal}{Numeric vector. Winsorized and downsampled poly(A)
#'       tail signal}
#'     \item{tail_moves}{Numeric vector. Downsampled basecaller moves for
#'       the tail region}
#'     \item{tail_pseudomoves}{Numeric vector. Pseudomove states computed
#'       from the tail signal (\{-1, 0, 1\})}
#'   }
#'   Always assign the returned list to a variable. Printing the full output
#'   to the console may crash your R session.
#'
#' @seealso
#' \code{\link{extract_polya_data}} for preparing the \code{polya_summary}
#' input, \code{\link{create_tail_feature_list}} for batch extraction,
#' \code{\link{filter_signal_by_threshold}} for pseudomove computation
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' ninetails::extract_tail_data(
#'   readname = 'abc123de-fg45-6789-0987-6543hijk2109',
#'   polya_summary = polya_summary_table,
#'   workspace = '/path/to/folder/containing/multifast5s',
#'   basecall_group = 'Basecall_1D_000'
#' )
#'
#' }
extract_tail_data <- function(readname,
                              polya_summary,
                              workspace,
                              basecall_group) {

  #Assertions
  if (missing(readname)) {
    stop(
      "Readname is missing. Please provide a valid readname argument.",
      call. = FALSE
    )
  }

  if (missing(workspace)) {
    stop(
      "Directory with basecalled fast5s is missing. Please provide a valid workspace argument.",
      call. = FALSE
    )
  }

  if (missing(basecall_group)) {
    stop(
      "Basecall group is missing. Please provide a valid basecall_group argument.",
      call. = FALSE
    )
  }

  if (missing(polya_summary)) {
    stop(
      "Polya_summary is missing. Please provide a valid polya_summary argument.",
      call. = FALSE
    )
  }

  if (!is.data.frame(polya_summary) || nrow(polya_summary) == 0) {
    stop(
      "Empty data frame provided as an input (polya_summary). Please provide a valid input table.",
      call. = FALSE
    )
  }

  assert_condition(
    is.character(workspace),
    "Path to basecalled fast5 files is not a character string. Please provide a valid path to basecalled fast5 files."
  )
  assert_condition(
    is_string(workspace, null.ok = FALSE, min.chars = 1),
    "Empty string provided as an input. Please provide a valid path to basecalled fast5 files."
  )
  assert_condition(
    is.character(readname),
    "Given readname is not a character string. Please provide a valid readname argument."
  )

  # Extract data from fast5 file
  fast5_filenames <- polya_summary$filename
  fast5_readname <- paste0("read_", readname) # suffix "read_"
  fast5_file <- fast5_filenames[readname] # particular file name (eg. FAO12345.fast5)
  fast5_file_path <- file.path(workspace, fast5_file) #path to browsed fast5 file

  signal <- rhdf5::h5read(
    file.path(fast5_file_path),
    paste0(fast5_readname, "/Raw/Signal")
  )
  signal <- as.vector(signal) #vectorized sig (natively: array)

  # retrieving moves
  BaseCalled_template <- rhdf5::h5read(
    fast5_file_path,
    paste0(fast5_readname, "/Analyses/", basecall_group, "/BaseCalled_template")
  )
  move <- BaseCalled_template$Move #how the basecaller "moves" through the called sequence, and allows for a mapping from basecall to raw data
  move <- as.numeric(move) # change data type as moves are originally stored as raw

  #read parameter stored in channel_id group
  channel_id <- rhdf5::h5readAttributes(
    fast5_file_path,
    paste0(fast5_readname, "/channel_id")
  ) # parent dir for attributes (within fast5 file)
  sampling_rate <- channel_id$sampling_rate # number of data points collected per second

  #read parameters (attrs) stored in basecall_1d_template
  basecall_1d_template <- rhdf5::h5readAttributes(
    fast5_file_path,
    paste0(
      fast5_readname,
      "/Analyses/",
      basecall_group,
      "/Summary/basecall_1d_template"
    )
  ) # parent dir for attributes (within fast5); fixed hardcoding bug
  stride <- basecall_1d_template$block_stride #  this parameter allows to sample data elements along a dimension
  called_events <- basecall_1d_template$called_events # number of events (nanopore translocations) recorded by device for given read
  number_of_events <- called_events * stride # number of events expanded for whole signal vec (this is estimation of signal length, however keep in mind that decimal values are ignored)

  # close all handled instances (groups, attrs) of fast5 file
  rhdf5::h5closeAll()

  # extract features from polya summary
  read_idx <- which(polya_summary$readname == readname) # this is to retrieve row number because tibbles cannot have rownames
  polya_start_position <- polya_summary$polya_start[read_idx]
  transcript_start_position <- polya_summary$transcript_start[read_idx]
  #define polya end position
  polya_end_position <- transcript_start_position - 1

  # handle move data
  #this is to expand moves along raw signal
  moves_sample_wise_vector <- c(
    rep(move, each = stride),
    rep(NA, length(signal) - number_of_events)
  )
  moves_tail_range <- moves_sample_wise_vector[
    polya_start_position:polya_end_position
  ]

  # extract polya tail region of the signal
  # signal is winsorized here!
  signal <- ninetails::winsorize_signal(signal[
    polya_start_position:polya_end_position
  ])

  # downsample (interpolate signal and moves) so the computation would be faster/easier to handle.
  signal <- round(
    stats::approx(
      signal,
      method = "linear",
      n = ceiling(0.2 * length(signal))
    )[[2]],
    digits = 0
  )
  moves_tail_range <- round(
    stats::approx(
      moves_tail_range,
      method = "linear",
      n = ceiling(0.2 * length(moves_tail_range))
    )[[2]],
    digits = 0
  )

  # filter signal to find local minima & maxima corresponding to potential C, G, U modifications
  pseudomoves <- ninetails::filter_signal_by_threshold(signal)

  extracted_data_single_list = list() # creating empty list for the extracted fast5 data

  extracted_data_single_list[["fast5_filename"]] <- polya_summary$filename[
    read_idx
  ]
  extracted_data_single_list[["tail_signal"]] <- signal
  extracted_data_single_list[["tail_moves"]] <- moves_tail_range
  #added for robustness
  extracted_data_single_list[["tail_pseudomoves"]] <- pseudomoves

  return(extracted_data_single_list)
}


#' Create list of poly(A) tail features from multi-Fast5 files
#'
#' Extracts tail features of RNA reads from multi-Fast5 files basecalled by
#' Guppy and poly(A) tail characteristics (coordinates) produced by nanopolish
#' polya function. Processing is parallelized across reads using
#' \pkg{foreach} and \pkg{doSNOW}. A progress bar is displayed during
#' extraction.
#'
#' After extraction, reads with zero-moved tails and reads that do not satisfy
#' the pseudomove condition (minimum run length of 5) are filtered out and
#' their identifiers are stored separately for downstream classification.
#'
#' @param nanopolish Character string or data frame. Full path of the
#'   \code{.tsv} file produced by nanopolish polya function, or an in-memory
#'   data frame.
#'
#' @param sequencing_summary Character string or data frame. Full path of the
#'   \code{.txt} file with sequencing summary, or an in-memory data frame.
#'
#' @param workspace Character string. Full path of the directory containing
#'   the basecalled multi-Fast5 files.
#'
#' @param num_cores Numeric. Number of physical cores to use in processing.
#'   Do not exceed 1 less than the number of cores at your disposal.
#'
#' @param basecall_group Character string. Name of the level in the Fast5
#'   file hierarchy from which data should be extracted (e.g.,
#'   \code{"Basecall_1D_000"}).
#'
#' @param pass_only Logical. If \code{TRUE} (default), only reads tagged by
#'   nanopolish as \code{"PASS"} are taken into consideration. If \code{FALSE},
#'   reads tagged as \code{"PASS"} and \code{"SUFFCLIP"} are both included.
#'
#' @return A named list with three elements:
#'   \describe{
#'     \item{tail_feature_list}{Named list of per-read tail features. Each
#'       element contains \code{fast5_filename}, \code{tail_signal},
#'       \code{tail_moves}, and \code{tail_pseudomoves} (see
#'       \code{\link{extract_tail_data}}).}
#'     \item{zeromoved_readnames}{Character vector. Read IDs discarded because
#'       all basecaller moves in their tail region were zero.}
#'     \item{nonpseudomoved_readnames}{Character vector. Read IDs discarded
#'       because their pseudomove chain was too short (< 5 consecutive
#'       positions) to indicate a potential modification.}
#'   }
#'   Always assign the returned list to a variable. Printing the full output
#'   to the console may crash your R session.
#'
#' @seealso
#' \code{\link{extract_polya_data}} for reading nanopolish and sequencing
#' summary data, \code{\link{extract_tail_data}} for single-read extraction,
#' \code{\link{create_tail_chunk_list}} for downstream chunk segmentation
#'
#' @importFrom foreach %dopar%
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' tfl <- ninetails::create_tail_feature_list(
#'   nanopolish = system.file('extdata',
#'                            'test_data',
#'                            'nanopolish_output.tsv',
#'                            package = 'ninetails'),
#'   sequencing_summary = system.file('extdata',
#'                                    'test_data',
#'                                    'sequencing_summary.txt',
#'                                    package = 'ninetails'),
#'   workspace = system.file('extdata',
#'                           'test_data',
#'                           'basecalled_fast5',
#'                           package = 'ninetails'),
#'   num_cores = 2,
#'   basecall_group = 'Basecall_1D_000',
#'   pass_only = TRUE
#' )
#'
#' }
#'
create_tail_feature_list <- function(nanopolish,
                                     sequencing_summary,
                                     workspace,
                                     num_cores,
                                     basecall_group,
                                     pass_only = TRUE) {

  # Assertions
  if (missing(num_cores)) {
    stop(
      "Number of declared cores is missing. Please provide a valid num_cores argument.",
      call. = FALSE
    )
  }

  if (missing(basecall_group)) {
    stop(
      "Basecall group is missing. Please provide a valid basecall_group argument.",
      call. = FALSE
    )
  }

  if (missing(workspace)) {
    stop(
      "Directory with basecalled fast5s (guppy workspace) is missing. Please provide a valid workspace argument.",
      call. = FALSE
    )
  }

  assert_condition(
    is.numeric(num_cores),
    "Declared core number must be numeric. Please provide a valid argument."
  )

  # Extracting and processing polya & sequencing summary data
  polya_summary <- ninetails::extract_polya_data(
    nanopolish,
    sequencing_summary,
    pass_only
  )

  #create empty list for extracted fast5 data
  tail_features_list = list()

  # creating cluster for parallel computing
  my_cluster <- parallel::makeCluster(num_cores)
  on.exit(parallel::stopCluster(my_cluster))
  doSNOW::registerDoSNOW(my_cluster)
  `%dopar%` <- foreach::`%dopar%`
  `%do%` <- foreach::`%do%`
  mc_options <- list(preschedule = TRUE, set.seed = FALSE, cleanup = FALSE)

  # header for progress bar
  cat(paste0(
    '[',
    as.character(Sys.time()),
    '] ',
    'Extracting features of provided reads...',
    '\n',
    sep = ''
  ))

  # progress bar
  pb <- utils::txtProgressBar(
    min = 0,
    max = length(polya_summary$filename),
    style = 3,
    width = 50,
    char = "=",
    file = stderr()
  )
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  # parallel extraction
  tail_features_list <- foreach::foreach(
    i = seq_along(polya_summary$readname),
    .combine = c,
    .inorder = TRUE,
    .errorhandling = 'pass',
    .options.snow = opts,
    .options.multicore = mc_options
  ) %dopar%
    {
      lapply(polya_summary$readname[i], function(x) {
        ninetails::extract_tail_data(
          x,
          polya_summary,
          workspace,
          basecall_group
        )
      })
    }

  #label each signal according to corresponding read name to avoid confusion
  # this deals with issue#5 (nanopolish records do not exactly match seqsummary file)
  squiggle_names <- as.vector(sapply(tail_features_list, function(x) {
    attributes(x[[1]])$names
  }))
  tail_features_list <- stats::setNames(tail_features_list, squiggle_names)

  # remove reads with only zero moved tails
  tail_features_list <- Filter(
    function(x) sum(x$tail_moves) != 0,
    tail_features_list
  )
  zeromoved_readnames <- squiggle_names[
    !(squiggle_names %in% names(tail_features_list))
  ]

  # prevent from running on reads which do not fulfill the pseudomove condition
  tail_features_list <- Filter(
    function(x) any(with(rle(x$tail_pseudomoves), lengths[values != 0] >= 5)),
    tail_features_list
  )

  # reads discarded because of unmet pseudomove condition
  #in this reads reported pseudomove chain is too short to be considered as potential modification
  nonpseudomoved_readnames <- squiggle_names[
    !(squiggle_names %in% c(zeromoved_readnames, names(tail_features_list)))
  ]

  #create final output
  tail_feature_list <- list()

  tail_feature_list[["tail_feature_list"]] <- tail_features_list
  tail_feature_list[["zeromoved_readnames"]] <- zeromoved_readnames
  tail_feature_list[["nonpseudomoved_readnames"]] <- nonpseudomoved_readnames

  # Done comm
  cat(paste0('[', as.character(Sys.time()), '] ', 'Done!', '\n', sep = ''))

  return(tail_feature_list)
}


#' Detect outliers (peaks and valleys) in ONT signal using z-scores
#'
#' Detects areas in the poly(A) tail signal where values significantly deviate
#' from those typical of an adenosine homopolymer. The function produces a
#' vector of pseudomoves with values in the range \{-1, 0, 1\}, where -1
#' corresponds to signals significantly below the local mean (potential C/U),
#' 1 corresponds to signals significantly above the local mean (potential G),
#' and 0 corresponds to typical adenosine homopolymer values.
#'
#' The pseudomoves vector allows more accurate calibration of nucleotide
#' positions of potential non-adenosine residues than the moves produced
#' by the Guppy basecaller.
#'
#' @param signal Numeric vector. An ONT read fragment corresponding to the
#'   poly(A) tail region of the read of interest, as delimited by nanopolish
#'   polya function. Fragments are stored in
#'   \code{tail_feature_list[[1]][[readname]][[2]]} produced by
#'   \code{\link{create_tail_feature_list}}.
#'
#' @details
#' The algorithm implements a sliding-window approach based on adaptive z-scores:
#' \enumerate{
#'   \item A calibration phase (100 synthetic data points drawn from the most
#'     frequent signal values) establishes baseline signal characteristics
#'   \item Rolling mean and standard deviation are computed over a 100-point
#'     adaptive window
#'   \item Outliers are detected when signal deviates > 3.5 standard deviations
#'     from the local mean
#'   \item Direction of deviation determines the pseudomove value
#' }
#'
#' Terminal positions (first 5) are masked to prevent false positives at the
#' adapter-tail boundary. This is a temporary safeguard until the
#' resegmentation model is refined.
#'
#' @section References:
#' Based on: Brakel, J.P.G. van (2014). "Robust peak detection algorithm
#' using z-scores". Stack Overflow. Available at:
#' \url{https://stackoverflow.com/questions/22583391/peak-signal-detection-in-realtime-timeseries-data/22640362#22640362}
#' (version: 2020-11-08).
#'
#' @return Numeric vector of pseudomoves corresponding to the analyzed tail
#'   region, with the same length as the input signal. Values are integers
#'   in the range \{-1, 0, 1\}:
#'   \describe{
#'     \item{-1}{Signal valleys (potential C/U nucleotides)}
#'     \item{0}{Normal signal (likely A nucleotides)}
#'     \item{1}{Signal peaks (potential G nucleotides)}
#'   }
#'
#' @seealso
#' \code{\link{filter_signal_by_threshold_vectorized}} for the optimized
#' Dorado version, \code{\link{substitute_gaps}} for gap handling,
#' \code{\link{create_tail_feature_list}} for the complete feature extraction
#' pipeline
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' ninetails::filter_signal_by_threshold(
#'   signal = tail_feature_list[[1]][["readname"]][[4]]
#' )
#'
#' }
filter_signal_by_threshold <- function(signal) {

  # Initialize vectors for incremental subset assignment
  baseline <- std_cutoff <- NULL

  #assertions
  if (missing(signal)) {
    stop(
      "Signal is missing. Please provide a valid signal argument [numeric vec].",
      call. = FALSE
    )
  }

  assert_condition(
    is.numeric(signal),
    "Signal must be numeric. Please provide a valid argument."
  )

  # reproducibility
  set.seed(123)

  # calibrate the algo on the sampled vals
  start_vals <- signal[1:10]
  most_freq_vals <- as.numeric(names(sort(table(signal), decreasing = TRUE)[
    1:20
  ]))
  adjusted_signal <- c(
    sample(c(most_freq_vals, start_vals), 100, replace = TRUE),
    signal
  )

  # Empyrical parameters:
  adaptive_sampling_window <- 100 # datapoints window for adjusting algo
  SD_threshold <- 3.5 # how many SD thresholds from avg signal pseudomove should be reported

  pseudomoves <- rep(0, length(adjusted_signal))
  filtered_signal <- adjusted_signal

  baseline[adaptive_sampling_window] <- mean(
    adjusted_signal[1:adaptive_sampling_window],
    na.rm = TRUE
  )
  std_cutoff[adaptive_sampling_window] <- stats::sd(
    adjusted_signal[1:adaptive_sampling_window],
    na.rm = TRUE
  )
  for (i in (adaptive_sampling_window + 1):length(adjusted_signal)) {
    if (
      abs(adjusted_signal[i] - baseline[i - 1]) >
        SD_threshold * std_cutoff[i - 1]
    ) {
      if (adjusted_signal[i] > baseline[i - 1]) {
        pseudomoves[i] <- 1 #if they go up
      } else {
        pseudomoves[i] <- -1 # if they go down
      }
      filtered_signal[i] <- filtered_signal[i - 1] #update
    } else {
      pseudomoves[i] <- 0 # uniform distr
      filtered_signal[i] <- adjusted_signal[i] #update
    }

    baseline[i] <- mean(
      filtered_signal[(i - adaptive_sampling_window):i],
      na.rm = TRUE
    )
    std_cutoff[i] <- stats::sd(
      filtered_signal[(i - adaptive_sampling_window):i],
      na.rm = TRUE
    )
  }

  # trim pseudomoves so they fit actual signal/moves
  adjusted_pseudomoves <- pseudomoves[101:length(pseudomoves)]

  #temporary hotfix for terminal false-positives (resegmentation model
  #- work in progress); current model was not trained on the terminal position,
  #therefore to avoid the increased rate of false-positives reported on the
  #adaptor-tail border, data at the very beginning are currently going to be
  #ignored; this is just a temporary solution, introduced to reduce
  #the bias which would be handled in more elegant way in the future releases
  adjusted_pseudomoves[1:5] <- 0

  adjusted_pseudomoves <- ninetails::substitute_gaps(adjusted_pseudomoves)

  return(adjusted_pseudomoves)
}

#' Extract modification-centered signal fragments from a poly(A) tail
#'
#' Finds areas in the poly(A) tail signal containing potential non-adenosine
#' residues and extracts 100-point signal fragments where the potential
#' modification is always at the center of a given extracted fragment.
#'
#' Candidate modification regions are identified based on two assumptions:
#' the presence of significant raw signal distortion (recorded as a
#' pseudomove by the thresholding algorithm) and the transition of state
#' (\code{move == 1}) recorded by Guppy. If only \code{move == 0} values
#' are present within a given signal chunk, then that chunk is dropped
#' from the analysis (the distortion is most likely caused by a sequencing
#' artifact, not a non-A residue itself).
#'
#' If the data indicating the presence of modifications are near the signal
#' ends (3' or 5'), missing upstream or downstream data are imputed based on
#' the most frequent values in the entire signal.
#'
#' @param readname Character string. Name of the given read (UUID) within
#'   the analyzed dataset.
#'
#' @param tail_feature_list List object produced by
#'   \code{\link{create_tail_feature_list}}. Must contain per-read entries
#'   with \code{tail_signal}, \code{tail_moves}, and
#'   \code{tail_pseudomoves}.
#'
#' @details
#' The extraction procedure is as follows:
#' \enumerate{
#'   \item Run-length encoding (RLE) of the pseudomove vector
#'   \item Filtering runs of pseudomoves with length >= 5
#'   \item Centering a 100-point window on the midpoint of each qualifying run
#'   \item Imputing NAs at boundaries with draws from the 5 most frequent
#'         signal values
#'   \item Removing chunks where basecaller moves are all zero (likely
#'         artifacts)
#' }
#'
#' @return A nested list where each element represents one candidate
#'   modification region. Each element is itself a list with:
#'   \describe{
#'     \item{chunk_sequence}{Numeric vector of raw signal values (length 100)}
#'     \item{chunk_start_pos}{Integer. Start index of the chunk in the
#'       original signal}
#'     \item{chunk_end_pos}{Integer. End index of the chunk in the
#'       original signal}
#'     \item{chunk_moves}{Numeric vector. Basecaller moves for the chunk
#'       region}
#'   }
#'   Element names follow the pattern \code{<readname>_<chunk_index>}.
#'   Positions are indexed from the 3' end.
#'
#' @seealso
#' \code{\link{create_tail_feature_list}} for preparing the input,
#' \code{\link{create_tail_chunk_list}} for batch segmentation,
#' \code{\link{split_tail_centered_dorado}} for the Dorado-specific version
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' ninetails::split_tail_centered(
#'   readname = "1234-anexample-r3adn4m3",
#'   tail_feature_list = tail_feature_list
#' )
#'
#' }
split_tail_centered <- function(readname, tail_feature_list) {

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

  #extract required data
  signal <- tail_feature_list[[1]][[readname]][[2]]
  moves <- tail_feature_list[[1]][[readname]][[3]]
  pseudomoves <- tail_feature_list[[1]][[readname]][[4]]

  # mod-centered chunk extraction
  # recompute rle
  mod_rle <- rle(pseudomoves)
  # pseudomoves filtered by condition (potentially decorated - empyrical!)
  condition <- mod_rle$lengths >= 5 & mod_rle$values
  # beginning positions of filtered pseudomoves which satisfy conditions
  first_filtered_positions <- cumsum(c(1, utils::head(mod_rle$lengths, -1)))[
    condition
  ]
  # length of pseudomoves satisfying condition
  filtered_length <- mod_rle$lengths[condition]
  # extracted coordinates (indices)
  start_positions <- first_filtered_positions + floor(filtered_length / 2) - 50
  end_positions <- start_positions + 99

  # extract signal chunks centered on potential modification
  list_1 <- lapply(1:length(start_positions), function(i) {
    if (start_positions[i] > 0) {
      signal[start_positions[i]:end_positions[i]]
    } else {
      c(rep(NA, abs(start_positions[i] - 1)), signal[1:end_positions[i]])
    }
  })

  # replace NAs with 5 most freq values (impute the data for gafs)
  most_freq_vals <- as.numeric(names(sort(table(signal), decreasing = TRUE)[
    1:5
  ]))
  list_1 <- lapply(list_1, function(n) {
    replace(n, is.na(n), sample(most_freq_vals, sum(is.na(n)), TRUE))
  })

  #add indices
  chunks_indices <- c(1:length(start_positions))

  # naming chunks based on names & indices
  chunk_names <- paste0(
    rep(readname, length(list_1)),
    '_',
    unlist(chunks_indices)
  )
  names(list_1) <- chunk_names

  # retrieve coordinates as list_2 and list_3:
  list_2 <- as.list(start_positions)
  #names(list_2) <- chunk_names
  list_3 <- as.list(end_positions)
  #names(list_3) <- chunk_names

  # retrieve move vector fragments based on the same positional info:
  list_4 <- lapply(1:length(start_positions), function(i) {
    if (start_positions[i] > 0) {
      moves[start_positions[i]:end_positions[i]]
    } else {
      c(rep(NA, abs(start_positions[i] - 1)), moves[1:end_positions[i]])
    }
  })

  out <- mget(ls(pattern = '^list.*\\d$')) %>%
    split(sub("_\\d+$", '', names(.))) %>%
    purrr::map(
      ~ purrr::transpose(purrr::set_names(
        .,
        c('chunk_sequence', 'chunk_start_pos', 'chunk_end_pos', 'chunk_moves')
      ))
    ) %>%
    purrr::flatten(.)

  # remove chunks with discordance between moves & pseudomoves reported
  # if pseudomoves registered while moves not - remove chunk
  # this usually causes the loss of terminal positions, but this is to avoid
  # inheritance of nanopolish error as well as the warn labeling
  out <- Filter(function(x) sum(x$chunk_moves) != 0, out)

  return(out)
}


#' Create list of poly(A) tail chunks centered on significant signal deviations
#'
#' Extracts fragments of poly(A) tails of ONT RNA reads potentially containing
#' non-A nucleotides, along with their coordinates, and appends the data to a
#' nested list organized by read IDs.
#'
#' Parallelization is handled with \pkg{foreach} and \pkg{doSNOW}. A progress
#' bar is displayed during processing. After extraction, empty list elements
#' are pruned and chunk moves are dropped to minimize memory footprint.
#'
#' @param tail_feature_list List object produced by
#'   \code{\link{create_tail_feature_list}}.
#'
#' @param num_cores Numeric. Number of physical cores to use in processing.
#'   Do not exceed 1 less than the number of cores at your disposal.
#'
#' @return A nested list containing the segmented tail data (chunks and
#'   coordinates) organized by read IDs. Each read entry contains one or more
#'   fragments, where each fragment is a list with:
#'   \describe{
#'     \item{chunk_sequence}{Numeric vector. Raw signal values (length 100)}
#'     \item{chunk_start_pos}{Integer. Starting index of the chunk in the
#'       original signal}
#'     \item{chunk_end_pos}{Integer. Ending index of the chunk in the
#'       original signal}
#'   }
#'
#' @seealso
#' \code{\link{create_tail_feature_list}} for preparing the input,
#' \code{\link{split_tail_centered}} for the per-read segmentation logic,
#' \code{\link{create_gaf_list}} for downstream GAF transformation
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' tcl <- ninetails::create_tail_chunk_list(
#'   tail_feature_list = tfl,
#'   num_cores = 3
#' )
#'
#' }
create_tail_chunk_list <- function(tail_feature_list, num_cores) {

  # initial assertions
  if (missing(num_cores)) {
    stop(
      "Number of declared cores is missing. Please provide a valid num_cores argument.",
      call. = FALSE
    )
  }

  if (missing(tail_feature_list)) {
    stop(
      "List of features is missing. Please provide a valid tail_feature_list argument.",
      call. = FALSE
    )
  }

  assert_condition(
    is.numeric(num_cores),
    "Declared core number must be numeric. Please provide a valid argument."
  )
  assert_condition(
    is.list(tail_feature_list),
    "Given tail_feature_list is not a list (class). Please provide valid file format."
  )

  # creating cluster for parallel computing
  my_cluster <- parallel::makeCluster(num_cores)
  on.exit(parallel::stopCluster(my_cluster))
  doSNOW::registerDoSNOW(my_cluster)
  `%dopar%` <- foreach::`%dopar%`
  `%do%` <- foreach::`%do%`
  mc_options <- list(preschedule = TRUE, set.seed = FALSE, cleanup = FALSE)

  # header for progress bar
  cat(paste0(
    '[',
    as.character(Sys.time()),
    '] ',
    'Creating tail segmentation data...',
    '\n',
    sep = ''
  ))

  # progress bar
  pb <- utils::txtProgressBar(
    min = 0,
    max = length(tail_feature_list[[1]]),
    style = 3,
    width = 50,
    char = "=",
    file = stderr()
  )
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  #output list
  tail_chunk_list <- list()

  #parallel extraction
  tail_chunk_list <- foreach::foreach(
    i = seq_along(tail_feature_list[[1]]),
    .combine = c,
    .inorder = TRUE,
    .errorhandling = 'pass',
    .options.snow = opts,
    .options.multicore = mc_options
  ) %dopar%
    {
      lapply(names(tail_feature_list[[1]][i]), function(x) {
        ninetails::split_tail_centered(x, tail_feature_list)
      })
    }

  #rename first level of nested list accordingly
  names(tail_chunk_list) <- names(tail_feature_list[[1]])

  #drop moves variable - recursive fn for prunning moves vector from the output
  # this is to save the memory space and minimize the size of vars
  .prune_moves <- function(i) {
    lapply(i, function(x) {
      if (is.list(x)) {
        if (!is.null(names(x))) {
          .prune_moves(x[names(x) != "chunk_moves"])
        } else {
          .prune_moves(x)
        }
      } else {
        x
      }
    })
  }

  tail_chunk_list <- .prune_moves(tail_chunk_list)

  #remove empty elements from the list (clean the output)
  tail_chunk_list <- rrapply::rrapply(
    tail_chunk_list,
    condition = Negate(is.null),
    how = "prune"
  )

  # Done comm
  cat(paste0('[', as.character(Sys.time()), '] ', 'Done!', '\n', sep = ''))

  return(tail_chunk_list)
}


#' Convert ONT signal to Gramian Angular Field
#'
#' Represents time series data (ONT squiggle) in a polar coordinate system
#' instead of the typical Cartesian coordinates. This is a pure-R equivalent
#' of the Python pyts.image implementation.
#'
#' Two methods of transformation are available: Gramian Angular Summation
#' Field (GASF) and Gramian Angular Difference Field (GADF).
#'
#' @param tail_chunk Numeric vector. A 100-element signal chunk representing
#'   a fragment of the poly(A) tail within the analyzed dataset.
#'
#' @param method Character string specifying the type of Gramian Angular Field.
#'   \code{"s"} (default) produces a summation field (GASF), \code{"d"}
#'   produces a difference field (GADF).
#'
#' @details
#' The transformation proceeds as follows:
#' \enumerate{
#'   \item Rescale signal values to the interval [-1, 1]
#'   \item Compute the inverse trigonometric function (\code{acos} for GASF,
#'         \code{asin} for GADF) of each value
#'   \item Form a pairwise matrix by summing (GASF) or differencing (GADF)
#'         the angle vectors
#'   \item Apply the corresponding trigonometric function (\code{cos} or
#'         \code{sin}) to the matrix
#'   \item Reshape to a 100x100x1 array
#'   \item Rescale result to the interval [0, 1]
#' }
#'
#' @return An array of dimensions (100, 100, 1) with values representing
#'   the Gramian Angular Field transformation of the input signal chunk.
#'   Values are in the range [0, 1].
#'
#' @seealso
#' \code{\link{combine_gafs}} for combining GASF and GADF into a single
#' array, \code{\link{create_gaf_list}} for batch GAF computation
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' ninetails::create_gaf(
#'   tail_chunk = tail_chunk,
#'   method = "s"
#' )
#'
#' }
create_gaf <- function(tail_chunk, method = "s") {

  #assertions
  if (missing(tail_chunk)) {
    stop(
      "Tail_chunk is missing. Please provide a valid tail_chunk argument.",
      call. = FALSE
    )
  }
  if (missing(method)) {
    stop(
      "Transformation method is missing. Please provide a valid method argument.",
      call. = FALSE
    )
  }
  assert_condition(
    is.numeric(tail_chunk),
    "Provided tail_chunk must be numeric. Please provide a valid argument."
  )
  assert_condition(
    is.character(method),
    "Provided method must be character string. Please provide a valid argument."
  )
  assert_condition(
    length(tail_chunk) == 100,
    "Provided chunks of wrong length. The chunk length should be equal to 100. Please provide a valid tail_chunk."
  )

  # rescale values so that all of them fall in the interval [-1, 1]:
  tail_chunk <- (tail_chunk -
    max(tail_chunk) +
    (tail_chunk - min(tail_chunk))) /
    (max(tail_chunk) - min(tail_chunk))

  if (method == "s") {
    # calculate phi coefficient for interpolation to polar coordinates
    # first by computing arc cosine
    # then by converting the data to matrix by replicating the vec
    # calculating sum of phi
    # calculating cosine
    #and reshaping the data into new dimensions
    tail_chunk <- acos(tail_chunk)
    tail_chunk <- cbind(replicate(length(tail_chunk), tail_chunk))
    tail_chunk <- tail_chunk + t(tail_chunk)
    tail_chunk <- cos(tail_chunk)
    tail_chunk <- array(t(tail_chunk), c(100, 100, 1))
  } else if (method == "d") {
    # computing arc sine instead of arc cosine
    # calculating sine instead of cosine
    tail_chunk <- asin(tail_chunk)
    tail_chunk <- cbind(replicate(length(tail_chunk), tail_chunk))
    tail_chunk <- tail_chunk + t(tail_chunk)
    tail_chunk <- sin(tail_chunk)
    tail_chunk <- array(t(tail_chunk), c(100, 100, 1))
  } else {
    stop(
      "Wrong GAF method definition. Please provide 's' for summation or 'd' for difference GAF."
    )
  }

  # rescale values so that all of them fall in the interval [0, 1]:
  tail_chunk <- round(
    (tail_chunk - min(tail_chunk)) / (max(tail_chunk) - min(tail_chunk)),
    4
  )

  return(tail_chunk)
}

#' Combine GASF and GADF into a two-channel array
#'
#' Creates a two-dimensional array containing both the Gramian Angular
#' Summation Field (GASF) and the Gramian Angular Difference Field (GADF)
#' produced from a single ONT tail chunk. Using both representations
#' increases the sensitivity of the CNN classification by overcoming the
#' limitations of each method individually.
#'
#' @param tail_chunk Numeric vector. A 100-element signal chunk representing
#'   a fragment of the poly(A) tail.
#'
#' @return An array of dimensions (100, 100, 2) where the first channel
#'   contains the GASF and the second channel contains the GADF, both
#'   produced by \code{\link{create_gaf}}.
#'
#' @seealso
#' \code{\link{create_gaf}} for individual GAF computation,
#' \code{\link{create_gaf_list}} for batch processing,
#' \code{\link{predict_gaf_classes}} for downstream CNN classification
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' ninetails::combine_gafs(tail_chunk = tail_chunk)
#'
#' }
combine_gafs <- function(tail_chunk) {

  #assertions
  if (missing(tail_chunk)) {
    stop(
      "Tail_chunk is missing. Please provide a valid tail_chunk argument.",
      call. = FALSE
    )
  }

  assert_condition(
    is.numeric(tail_chunk),
    "Provided tail_chunk must be numeric. Please provide a valid argument."
  )

  #create gasf & gaf
  GASF <- ninetails::create_gaf(tail_chunk = tail_chunk, method = "s")
  GADF <- ninetails::create_gaf(tail_chunk = tail_chunk, method = "d")

  #create array of gafs
  combined_gafs <- array(c(GASF, GADF), dim = c(100, 100, 2))

  return(combined_gafs)
}


#' Create list of Gramian Angular Field matrices from tail chunks
#'
#' Computes Gramian Angular Field (GAF) representations for all signal chunks
#' in a tail chunk list. Each chunk is transformed into a combined GASF + GADF
#' array via \code{\link{combine_gafs}}.
#'
#' Parallelization is handled with \pkg{foreach} and \pkg{doSNOW}, and a
#' progress bar is displayed during processing.
#'
#' @param tail_chunk_list List object produced by
#'   \code{\link{create_tail_chunk_list}}.
#'
#' @param num_cores Numeric. Number of physical cores to use in processing.
#'   Do not exceed 1 less than the number of cores at your disposal.
#'
#' @return A named list of GAF arrays (100, 100, 2) organized by read
#'   \code{ID_index}. Always assign the returned list to a variable.
#'   Printing the full output to the console may crash your R session.
#'
#' @seealso
#' \code{\link{create_tail_chunk_list}} for preparing the input,
#' \code{\link{combine_gafs}} for the per-chunk GAF transformation,
#' \code{\link{predict_gaf_classes}} for downstream CNN classification
#'
#' @importFrom foreach %dopar%
#' @importFrom utils head
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' gl <- ninetails::create_gaf_list(
#'   tail_chunk_list = tcl,
#'   num_cores = 2
#' )
#'
#' }
create_gaf_list <- function(tail_chunk_list, num_cores) {

  # Assertions
  if (missing(num_cores)) {
    stop(
      "Number of declared cores is missing. Please provide a valid num_cores argument.",
      call. = FALSE
    )
  }

  if (missing(tail_chunk_list)) {
    stop(
      "List of tail chunks is missing. Please provide a valid tail_chunk_list argument.",
      call. = FALSE
    )
  }

  assert_condition(
    is.list(tail_chunk_list),
    "Given tail_chunk_list is not a list (class). Please provide valid file format."
  )
  assert_condition(
    is.numeric(num_cores),
    "Declared core number must be numeric. Please provide a valid argument."
  )

  # creating cluster for parallel computing
  my_cluster <- parallel::makeCluster(num_cores)
  on.exit(parallel::stopCluster(my_cluster))
  doSNOW::registerDoSNOW(my_cluster)
  `%dopar%` <- foreach::`%dopar%`
  `%do%` <- foreach::`%do%`
  mc_options <- list(preschedule = TRUE, set.seed = FALSE, cleanup = FALSE)

  #create empty list for the data
  gaf_list = list()

  #progressbar header
  cat(paste0(
    '[',
    as.character(Sys.time()),
    '] ',
    'Computing gramian angular fields...',
    '\n',
    sep = ''
  ))

  # progress bar
  pb <- utils::txtProgressBar(
    min = 0,
    max = length(tail_chunk_list),
    style = 3,
    width = 50,
    char = "=",
    file = stderr()
  )
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  gaf_list <- foreach::foreach(
    i = seq_along(tail_chunk_list),
    .combine = c,
    .inorder = TRUE,
    .errorhandling = 'pass',
    .options.snow = opts,
    .options.multicore = mc_options
  ) %dopar%
    {
      lapply(tail_chunk_list[[i]], function(x_ij) {
        ninetails::combine_gafs(x_ij[['chunk_sequence']])
      })
    }

  # Done comm
  cat(paste0('[', as.character(Sys.time()), '] ', 'Done!', '\n', sep = ''))

  return(gaf_list)
}

#' Classify Gramian Angular Field matrices with a pretrained CNN
#'
#' Assigns GAF representations of signal chunks to one of four nucleotide
#' categories (A, C, G, U) using a pretrained convolutional neural network.
#' The model is loaded via \code{\link{load_keras_model}} and inference is
#' performed through the TensorFlow/Keras backend.
#'
#' @param gaf_list List of GAF arrays (100, 100, 2), as produced by
#'   \code{\link{create_gaf_list}}.
#'
#' @details
#' The function reshapes the input GAF list into a 4-D tensor
#' (n_chunks x 100 x 100 x 2) and applies the pretrained CNN to predict
#' the nucleotide class for each chunk. Predictions are returned as integer
#' codes: 0 = A, 1 = C, 2 = G, 3 = U.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{chunkname}{Character vector. Names of the classified signal
#'       chunks (format: \code{<readname>_<index>})}
#'     \item{prediction}{Integer vector. Predicted nucleotide class for
#'       each chunk (0 = A, 1 = C, 2 = G, 3 = U)}
#'   }
#'
#' @seealso
#' \code{\link{create_gaf_list}} for preparing the input,
#' \code{\link{load_keras_model}} for model loading,
#' \code{\link{create_outputs}} for integrating predictions into final output
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' pl <- ninetails::predict_gaf_classes(gl)
#'
#' }
predict_gaf_classes <- function(gaf_list) {

  #assertions
  if (missing(gaf_list)) {
    stop(
      "List of transformed signal chunks is missing. Please provide a valid gaf_list argument.",
      call. = FALSE
    )
  }

  chunknames <- names(gaf_list)
  names(gaf_list) <- NULL

  gaf_list <- simplify2array(gaf_list)
  gaf_list <- aperm(gaf_list, c(4, 1, 2, 3))

  # Output info
  cat(paste0(
    '[',
    as.character(Sys.time()),
    '] ',
    'Classifying gramian angular fields...',
    '\n',
    sep = ''
  ))

  keras_model <- ninetails::load_keras_model()

  #predict chunk class
  predicted_gaf_classes <- keras_model %>%
    stats::predict(gaf_list) %>%
    keras::k_argmax()
  predicted_gaf_classes <- as.numeric(predicted_gaf_classes)

  predicted_list = list() # creating empty list for the extracted  data

  predicted_list[["chunkname"]] <- chunknames
  predicted_list[["prediction"]] <- predicted_gaf_classes

  # Output info
  cat(paste0('[', as.character(Sys.time()), '] ', 'Done!', '\n', sep = ''))

  return(predicted_list)
}


#' Create ninetails output tables (Guppy legacy pipeline)
#'
#' Integrates nanopolish poly(A) tail data, CNN predictions, and positional
#' information to generate two main outputs: per-read classifications and
#' non-adenosine residue predictions with estimated positions along the
#' poly(A) tail.
#'
#' The function implements a complete read accounting system that classifies
#' all reads from the nanopolish output into biologically meaningful categories
#' based on both quality control metrics and modification detection results.
#'
#' @param tail_feature_list List object produced by
#'   \code{\link{create_tail_feature_list}}.
#'
#' @param tail_chunk_list List object produced by
#'   \code{\link{create_tail_chunk_list}}.
#'
#' @param nanopolish Character string or data frame. Full path of the
#'   \code{.tsv} file produced by nanopolish polya function, or an in-memory
#'   data frame.
#'
#' @param predicted_list List object produced by
#'   \code{\link{predict_gaf_classes}}.
#'
#' @param num_cores Numeric. Number of physical cores to use in processing.
#'   Do not exceed 1 less than the number of cores at your disposal.
#'
#' @param pass_only Logical. If \code{TRUE} (default), only reads tagged by
#'   nanopolish as \code{"PASS"} are taken into consideration. If \code{FALSE},
#'   reads tagged as \code{"PASS"} and \code{"SUFFCLIP"} are both included.
#'
#' @param qc Logical. If \code{TRUE} (default), quality control of the output
#'   predictions is performed. Reads with non-A residue positions in terminal
#'   nucleotides (< 2 nt from either end of the tail) are labeled with
#'   \code{"-WARN"} suffixes, as these are most likely artifacts inherited
#'   from nanopolish segmentation. It is then up to the user whether to
#'   include or discard such reads from downstream analysis.
#'
#' @section Read Classification System:
#'   Reads are assigned into three main categories with specific comment codes:
#'
#'   \strong{decorated} - Reads with detected non-adenosine modifications
#'   \itemize{
#'     \item YAY: Move transition present, nonA residue detected
#'   }
#'
#'   \strong{blank} - Reads without detected modifications
#'   \itemize{
#'     \item MAU: Move transition absent, nonA residue undetected
#'     \item MPU: Move transition present, nonA residue undetected
#'   }
#'
#'   \strong{unclassified} - Reads that failed quality control
#'   \itemize{
#'     \item IRL: Insufficient read length (poly(A) < 10 nt)
#'     \item QCF: Nanopolish QC failed
#'     \item NIN: Not included in the analysis (\code{pass_only = TRUE})
#'   }
#'
#' @section Position Estimation:
#'   Non-adenosine residue positions are estimated using the formula:
#'
#'   \code{est_nonA_pos = polya_length - ((polya_length * centr_signal_pos) /
#'   signal_length)}
#'
#'   Positions are reported as distance from the 3' end of the tail
#'   (position 1 = most 3' nucleotide).
#'
#' @return A named list with two data frames:
#'   \describe{
#'     \item{read_classes}{Data frame with per-read classification results.
#'       Columns: \code{readname}, \code{contig}, \code{polya_length},
#'       \code{qc_tag}, \code{class}, and \code{comments}.}
#'     \item{nonadenosine_residues}{Data frame with per-chunk predictions of
#'       non-adenosine residues for decorated reads only. Columns:
#'       \code{readname}, \code{contig}, \code{prediction}, \code{est_nonA_pos},
#'       \code{polya_length}, \code{qc_tag}. When \code{qc = TRUE}, prediction
#'       values may carry a \code{"-WARN"} suffix for terminal positions.}
#'   }
#'
#' @seealso
#' \code{\link{create_tail_feature_list}} for feature extraction,
#' \code{\link{create_tail_chunk_list}} for chunk segmentation,
#' \code{\link{predict_gaf_classes}} for CNN classification,
#' \code{\link{check_tails_guppy}} for the complete pipeline wrapper
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' create_outputs(
#'   tail_feature_list = tail_feature_list,
#'   tail_chunk_list = tail_chunk_list,
#'   nanopolish = '/path/to/nanopolish_output.tsv',
#'   predicted_list = predicted_list,
#'   num_cores = 2,
#'   pass_only = TRUE,
#'   qc = TRUE
#' )
#'
#' }
#'
#'
create_outputs <- function(tail_feature_list,
                           tail_chunk_list,
                           nanopolish,
                           predicted_list,
                           num_cores,
                           pass_only = TRUE,
                           qc = TRUE) {

  #assertions
  if (missing(tail_feature_list)) {
    stop(
      "List of tail features is missing. Please provide a valid tail_feature_list argument.",
      call. = FALSE
    )
  }

  if (missing(tail_chunk_list)) {
    stop(
      "List of tail chunks is missing. Please provide a valid tail_chunk_list argument.",
      call. = FALSE
    )
  }

  if (missing(nanopolish)) {
    stop(
      "Nanopolish polya output is missing. Please provide a valid nanopolish argument.",
      .call = FALSE
    )
  }

  if (missing(predicted_list)) {
    stop(
      "List of predictions is missing. Please provide a valid predicted_list argument.",
      call. = FALSE
    )
  }

  if (missing(num_cores)) {
    stop(
      "Number of declared cores is missing. Please provide a valid num_cores argument.",
      call. = FALSE
    )
  }

  assert_condition(
    is.list(tail_feature_list),
    "Given tail_feature_list is not a list (class). Please provide valid object."
  )
  assert_condition(
    is.list(tail_chunk_list),
    "Given tail_chunk_list is not a list (class). Please provide valid object."
  )
  assert_condition(
    is.list(predicted_list),
    "Given predicted_list is not a list (class). Please provide valid object."
  )

  if (is_string(nanopolish)) {
    # if string provided as an argument, read from file
    assert_file_exists(nanopolish, "Nanopolish")
    nanopolish_polya_table <- vroom::vroom(
      nanopolish,
      col_select = c(readname, contig, polya_length, qc_tag),
      show_col_types = FALSE
    )
  } else {
    # make sure that nanopolish is an object with rows
    if (!is.data.frame(nanopolish) || nrow(nanopolish) == 0) {
      stop(
        "Empty data frame provided as an input (nanopolish). Please provide valid input"
      )
    }

    nanopolish_polya_table <- nanopolish[, c(
      "readname",
      "contig",
      "polya_length",
      "qc_tag"
    )]
  }

  # HANDLE LENGTH/POSITION CALIBRATING DATA
  # creating cluster for parallel computing
  my_cluster <- parallel::makeCluster(num_cores)
  on.exit(parallel::stopCluster(my_cluster))

  doSNOW::registerDoSNOW(my_cluster)
  `%dopar%` <- foreach::`%dopar%`
  `%do%` <- foreach::`%do%`
  mc_options <- list(preschedule = TRUE, set.seed = FALSE, cleanup = FALSE)

  #create empty list for the data
  tail_length_list = list()

  # header for progress bar
  cat(paste0(
    '[',
    as.character(Sys.time()),
    '] ',
    'Retrieving estimated length data...',
    '\n',
    sep = ''
  ))

  #set progressbar
  pb <- utils::txtProgressBar(
    min = 0,
    max = length(names(tail_feature_list[[1]])),
    style = 3,
    width = 50,
    char = "=",
    file = stderr()
  )
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  tail_length_list <- foreach::foreach(
    i = seq_along(tail_feature_list[[1]]),
    .combine = c,
    .inorder = TRUE,
    .errorhandling = 'pass',
    .options.snow = opts,
    .options.multicore = mc_options
  ) %dopar%
    {
      lapply(names(tail_feature_list[[1]][i]), function(x) {
        length(tail_feature_list[[1]][[x]][[2]])
      })
    }

  #coerce to df, add names
  tail_length_list <- do.call("rbind.data.frame", tail_length_list)
  tail_length_list$readname <- names(tail_feature_list[["tail_feature_list"]])
  colnames(tail_length_list) <- c("signal_length", "readname")

  #merge data from feature list with nanopolish estimations
  tails_tail_feature_list <- dplyr::left_join(
    tail_length_list,
    nanopolish_polya_table,
    by = "readname"
  )

  ### Chunk positional data
  #create empty list for extracted data
  non_a_position_list <- list()

  # header for progress bar
  cat(paste0(
    '[',
    as.character(Sys.time()),
    '] ',
    'Retrieving position calibrating data...',
    '\n',
    sep = ''
  ))

  #set progressbar
  # progress bar
  pb <- utils::txtProgressBar(
    min = 0,
    max = length(tail_chunk_list),
    style = 3,
    width = 50,
    char = "=",
    file = stderr()
  )
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  non_a_position_list <- foreach::foreach(
    i = seq_along(tail_chunk_list),
    .combine = c,
    .inorder = TRUE,
    .errorhandling = 'pass',
    .options.snow = opts,
    .options.multicore = mc_options
  ) %dopar%
    {
      lapply(tail_chunk_list[[i]], function(x) x[['chunk_start_pos']] + 50)
    }
  #close(pb)

  # coerce to df, add names
  non_a_position_list <- do.call("rbind.data.frame", non_a_position_list)
  chunknames <- purrr::map_depth(tail_chunk_list, 1, names) %>%
    unlist(use.names = F)
  non_a_position_list$chunkname <- chunknames
  non_a_position_list <- non_a_position_list %>%
    dplyr::mutate(readname = gsub('_.*', '', chunkname))
  colnames(non_a_position_list)[1] <- c("centr_signal_pos")

  #merge data from feature list with nanopolish estimations
  non_a_position_list <- dplyr::left_join(
    non_a_position_list,
    tails_tail_feature_list,
    by = "readname"
  )

  # HANDLE PREDICTIONS
  moved_chunks_table <- data.frame(t(Reduce(rbind, predicted_list)))
  # rename cols
  colnames(moved_chunks_table) <- c("chunkname", "prediction")
  # extract readnames
  moved_chunks_table$readname <- sub('\\_.*', '', moved_chunks_table$chunkname)

  # vectorized substitution of predictions with letter code
  prediction_dict <- c("0" = "A", "1" = "C", "2" = "G", "3" = "U")
  moved_chunks_table$prediction <- prediction_dict[as.character(
    moved_chunks_table$prediction
  )]

  #extract reads with move==1 and modification absent (not detected)
  moved_blank_readnames <- names(which(with(
    moved_chunks_table,
    tapply(prediction, readname, unique) == 'A'
  )))

  # cleaned chunks_table
  moved_chunks_table <- subset(
    moved_chunks_table,
    !(readname %in% moved_blank_readnames)
  )
  # delete A-containing rows
  moved_chunks_table <- moved_chunks_table[
    !(moved_chunks_table$prediction == "A"),
  ]
  #merge data from feats & predictions
  moved_chunks_table <- dplyr::left_join(
    moved_chunks_table,
    non_a_position_list,
    by = c("readname", "chunkname")
  )

  #estimate non-A nucleotide position
  moved_chunks_table$est_nonA_pos <- round(
    moved_chunks_table$polya_length -
      ((moved_chunks_table$polya_length * moved_chunks_table$centr_signal_pos) /
        moved_chunks_table$signal_length),
    digits = 2
  )

  #clean up the output nonA table:
  moved_chunks_table <- moved_chunks_table[, c(3, 6, 2, 9, 7, 8)]

  # Handle other (discarded) reads:
  discarded_reads <- nanopolish_polya_table[
    !nanopolish_polya_table$readname %in% moved_chunks_table$readname,
  ]

  # Add filtering criterion: select only pass or pass $ suffclip
  if (pass_only == TRUE) {
    discarded_reads <- discarded_reads %>%
      dplyr::filter(!readname %in% moved_chunks_table$readname) %>%
      dplyr::mutate(
        comments = dplyr::case_when(
          polya_length < 10 ~ "IRL",
          qc_tag == "SUFFCLIP" ~ "NIN",
          qc_tag == "ADAPTER" ~ "QCF",
          qc_tag == "NOREGION" ~ "QCF",
          qc_tag == "READ_FAILED_LOAD" ~ "QCF",
          readname %in% moved_blank_readnames ~ "MPU",
          TRUE ~ "MAU"
        ),
        class = dplyr::case_when(
          polya_length < 10 ~ "unclassified",
          readname %in% moved_blank_readnames ~ "blank",
          comments == "MAU" ~ "blank",
          TRUE ~ "unclassified"
        )
      )
  } else {
    discarded_reads <- discarded_reads %>%
      dplyr::filter(!readname %in% moved_chunks_table$readname) %>%
      dplyr::mutate(
        comments = dplyr::case_when(
          polya_length < 10 ~ "IRL",
          qc_tag == "ADAPTER" ~ "QCF",
          qc_tag == "NOREGION" ~ "QCF",
          qc_tag == "READ_FAILED_LOAD" ~ "QCF",
          readname %in% moved_blank_readnames ~ "MPU",
          TRUE ~ "MAU"
        ),
        class = dplyr::case_when(
          polya_length < 10 ~ "unclassified",
          readname %in% moved_blank_readnames ~ "blank",
          comments == "MAU" ~ "blank",
          TRUE ~ "unclassified"
        )
      )
  }

  decorated_reads <- nanopolish_polya_table[
    nanopolish_polya_table$readname %in% moved_chunks_table$readname,
  ]

  decorated_reads$class <- "decorated"
  decorated_reads$comments <- "YAY"

  #merge read_classes tabular output:
  nanopolish_polya_table <- rbind(decorated_reads, discarded_reads)
  #corece tibble to df
  nanopolish_polya_table <- data.frame(nanopolish_polya_table)

  #create empty list for the output
  ninetails_output <- list()

  ##############################################################################
  # QUALITY CONTROL OF THE OUTPUTS
  ##############################################################################
  # the model was not trained on the tail termini; also the ninetails inherits
  # the potential segmentation errors from nanopolish (as the tail is delimited
  # based on nanopolish polya function). Thus, to avoid artifact contribution in
  # the final output tables, the potential erroneous positions indicated by the
  # classifier based on the signal shape would be labeled by this module; it
  # then depends on the user whether to include such reads in the analysis or
  # not; additional info would be provided as "WARN" comment to treat those reads
  # with caution.

  if (qc == TRUE) {
    # filter the outermost positions (termini!) from position data:
    # empirically tested constraints! report without terminal data
    moved_chunks_table_trimmed <- moved_chunks_table %>%
      dplyr::filter(!(est_nonA_pos < 2) & !(est_nonA_pos > polya_length - 2))
    moved_chunks_table_discarded <- subset(
      moved_chunks_table,
      !(readname %in% moved_chunks_table_trimmed$readname)
    )

    # subset potential artifacts from read_classes
    potential_artifacts <- subset(
      nanopolish_polya_table,
      readname %in% moved_chunks_table_discarded$readname
    )

    decorated_reads_edited <- nanopolish_polya_table %>%
      dplyr::mutate(
        class = dplyr::case_when(
          readname %in% potential_artifacts$readname ~ paste0(class, "-WARN"),
          TRUE ~ paste0(class)
        )
      )

    # label potential artifacts in nonadenosine residue dataframe
    moved_chunks_table_qc <- moved_chunks_table %>%
      dplyr::mutate(
        prediction = dplyr::case_when(
          est_nonA_pos < 2 ~ paste0(prediction, "-WARN"),
          est_nonA_pos > polya_length - 2 ~ paste0(prediction, "-WARN"),
          TRUE ~ paste0(prediction)
        )
      )

    #CREATE FINAL OUTPUT
    #prevent potential bugs inherited from nanopolish multimapping
    decorated_reads_edited <- unique(decorated_reads_edited)
    moved_chunks_table_qc <- unique(moved_chunks_table_qc)

    ninetails_output[['read_classes']] <- decorated_reads_edited
    ninetails_output[['nonadenosine_residues']] <- moved_chunks_table_qc
  } else {
    #CREATE FINAL OUTPUT
    #prevent potential bugs inherited from nanopolish multimapping
    nanopolish_polya_table <- unique(nanopolish_polya_table)
    moved_chunks_table <- unique(moved_chunks_table)

    ninetails_output[['read_classes']] <- nanopolish_polya_table
    ninetails_output[['nonadenosine_residues']] <- moved_chunks_table
  }

  return(ninetails_output)
}
