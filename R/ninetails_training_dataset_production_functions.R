################################################################################
# FUNCTIONS DEVELOPED TO HANDLA GUPPY DATA
################################################################################
#' Extracts tail features of single RNA read from respective basecalled
#' multi-fast5 file.
#'
#' This is the training-set variant of the single-read extraction function.
#' It reads raw signal, basecaller moves and channel metadata from a
#' multi-Fast5 file, isolates the poly(A) tail region defined by nanopolish,
#' winsorizes the signal, downsamples by linear interpolation (to 20\% of
#' original length), and computes pseudomoves via
#' \code{\link{filter_signal_by_threshold_trainingset}}.
#'
#' @details
#' The function differs from its production counterpart
#' (\code{\link{extract_tail_data}}) in that:
#' \itemize{
#'   \item The signal is downsampled (interpolated) to 20\% of its original
#'         length to speed up downstream training-set preparation.
#'   \item Pseudomoves are computed using the training-set-specific threshold
#'         filter (\code{\link{filter_signal_by_threshold_trainingset}}).
#'   \item The returned list omits production-only fields and instead provides
#'         elements suitable for chunk splitting and GAF creation.
#' }
#'
#' @param readname Character string. Name (UUID) of the given read within
#'   the analyzed dataset.
#'
#' @param polya_summary Data frame. The table containing data extracted
#'   from nanopolish and the sequencing summary (produced by
#'   \code{\link{extract_polya_data}}).
#'
#' @param workspace Character string. Full path of the directory containing
#'   basecalled multi-Fast5 files.
#'
#' @param basecall_group Character string \code{["Basecall_1D_000"]}. Name of
#'   the level in the Fast5 file hierarchy from which the data should be
#'   extracted.
#'
#' @return A named list with four elements:
#' \describe{
#'   \item{fast5_filename}{Character. Name of the source Fast5 file.}
#'   \item{tail_signal}{Numeric vector. Winsorized and downsampled signal
#'     corresponding to the poly(A) tail region.}
#'   \item{tail_moves}{Numeric vector. Downsampled basecaller moves for the
#'     tail region.}
#'   \item{tail_pseudomoves}{Numeric vector. Pseudomove values (-1, 0, 1)
#'     indicating potential non-A deviations in the tail signal.}
#' }
#' Always assign this returned list to a variable; printing the full list
#' to the console may crash the R session.
#'
#' @seealso \code{\link{create_tail_feature_list_trainingset}} for the
#'   parallel wrapper that calls this function,
#'   \code{\link{extract_tail_data}} for the production counterpart,
#'   \code{\link{winsorize_signal}} for the signal clipping step,
#'   \code{\link{filter_signal_by_threshold_trainingset}} for pseudomove
#'   detection.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' extract_tail_data_trainingset(
#'   readname = 'abc123de-fg45-6789-0987-6543hijk2109',
#'   polya_summary = polya_summary_table,
#'   workspace = '/path/to/folder/containing/multifast5s',
#'   basecall_group = 'Basecall_1D_000')
#'
#' }
#'
#'
extract_tail_data_trainingset <- function(readname,
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
      "Empty data frame provided as an input (polya_summary). Please provide valid input"
    )
  }

  assert_condition(
    !is.na(workspace) && nchar(workspace) > 0,
    "Empty string provided as an input. Please provide a valid path to basecalled fast5 files."
  )

  assert_condition(
    is.character(workspace),
    "Path to basecalled fast5 files is not a character string. Please provide a valid path to basecalled fast5 files."
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
  #move <- factor(move, levels=c(0,1))

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
      "/Analyses/Basecall_1D_000/Summary/basecall_1d_template"
    )
  ) # parent dir for attributes (within fast5)
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
  pseudomoves <- filter_signal_by_threshold_trainingset(signal)

  extracted_data_single_list = list() # creating empty list for the extracted fast5 data

  extracted_data_single_list[["fast5_filename"]] <- polya_summary$filename[
    read_idx
  ]
  extracted_data_single_list[["tail_signal"]] <- signal
  extracted_data_single_list[["tail_moves"]] <- moves_tail_range
  extracted_data_single_list[["tail_pseudomoves"]] <- pseudomoves

  return(extracted_data_single_list)
}


#' Extracts features of poly(A) tails of ONT RNA reads required for finding
#' non-A nucleotides within the given tails.
#'
#' This is the training-set variant of the feature extraction wrapper. It
#' processes all reads in parallel, calling
#' \code{\link{extract_tail_data_trainingset}} for each read, then filters
#' out reads with zero-moved tails or pseudomove chains too short to
#' represent potential modifications.
#'
#' @details
#' The function differs from its production counterpart in that it retains
#' reads whose pseudomove chains satisfy a length >= 4 criterion, which is
#' required for subsequent modification-centered chunk splitting. Two
#' categories of discarded reads are tracked (zero-moved and
#' non-pseudomoved) and returned alongside the valid feature list.
#'
#' @param nanopolish Character string. Full path of the \code{.tsv} file
#'   produced by \code{nanopolish polya}.
#'
#' @param sequencing_summary Character string. Full path of the \code{.txt}
#'   file with the sequencing summary.
#'
#' @param workspace Character string. Full path of the directory containing
#'   basecalled multi-Fast5 files.
#'
#' @param num_cores Numeric \code{[1]}. Number of physical cores to use.
#'   Do not exceed 1 less than the number of cores at your disposal.
#'
#' @param basecall_group Character string \code{["Basecall_1D_000"]}. Name of
#'   the level in the Fast5 file hierarchy from which the data should be
#'   extracted.
#'
#' @param pass_only Logical \code{[TRUE]}. If \code{TRUE}, only reads tagged
#'   by nanopolish as \code{"PASS"} are retained. Otherwise, reads tagged as
#'   \code{"PASS"} or \code{"SUFFCLIP"} are included.
#'
#' @return A named list with three elements:
#' \describe{
#'   \item{tail_feature_list}{Named list of per-read feature lists (as
#'     returned by \code{\link{extract_tail_data_trainingset}}).}
#'   \item{zeromoved_readnames}{Character vector. Read IDs discarded because
#'     their tail moves summed to zero.}
#'   \item{nonpseudomoved_readnames}{Character vector. Read IDs discarded
#'     because their pseudomove chains were too short (< 4).}
#' }
#' Always assign this returned list to a variable; printing the full list
#' to the console may crash the R session.
#'
#' @seealso \code{\link{extract_tail_data_trainingset}} for the per-read
#'   extraction step,
#'   \code{\link{create_tail_chunk_list_trainingset}} for the next pipeline
#'   step,
#'   \code{\link{extract_polya_data}} for input data preparation.
#'
#' @importFrom foreach %dopar%
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' create_tail_feature_list_trainingset(
#'   nanopolish = '/path/to/file',
#'   sequencing_summary = '/path/to/file',
#'   workspace = '/path/to/guppy/workspace',
#'   num_cores = 10,
#'   basecall_group = 'Basecall_1D_000')
#'
#' }
#'
create_tail_feature_list_trainingset <- function(nanopolish,
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
        ninetails::extract_tail_data_trainingset(
          x,
          polya_summary,
          workspace,
          basecall_group
        )
      })
    }


  #label each signal according to corresponding read name to avoid confusion
  squiggle_names <- polya_summary$readname
  names(tail_features_list) <- squiggle_names

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
    function(x) any(with(rle(x$tail_pseudomoves), lengths[values != 0] >= 4)),
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


#' Detection of outliers (peaks & valleys) in ONT signal using z-scores.
#'
#' This is the training-set variant of the z-score-based peak/valley
#' detection function. It identifies deviations in the poly(A) tail
#' signal that may correspond to non-adenosine nucleotides (C, G, or U)
#' by applying a moving-average z-score filter with empirically calibrated
#' parameters.
#'
#' @details
#' The algorithm prepends 100 calibration data points (sampled from the
#' 10 most frequent signal values and the first 10 observations) to
#' stabilise the baseline estimate. It then applies a sliding window of
#' 100 data points: any observation exceeding 3.5 standard deviations
#' from the running mean is classified as a peak (+1) or valley (-1);
#' all other positions are scored 0. After trimming the calibration
#' prefix, terminal positions are zeroed out (first 5 and last 5) to
#' prevent edge artefacts from entering downstream chunk extraction.
#' Finally, \code{\link{substitute_gaps}} is applied to fill isolated
#' zero-gaps in the pseudomove vector.
#'
#' @section Acknowledgements:
#' The z-score peak detection algorithm is based on:
#' Brakel, J.P.G. van (2014). \emph{"Robust peak detection algorithm
#' using z-scores"}. Stack Overflow. Available at:
#' \url{https://stackoverflow.com/questions/22583391/peak-signal-detection-in-realtime-timeseries-data/22640362#22640362}
#' (version: 2020-11-08).
#'
#' @param signal Numeric vector. An ONT read fragment corresponding to the
#'   poly(A) tail region as delimited by nanopolish polya (stored in
#'   \code{tail_feature_list[[1]]} produced by
#'   \code{\link{create_tail_feature_list_trainingset}}).
#'
#' @return Numeric vector of pseudomoves with values in \{-1, 0, 1\}:
#' \describe{
#'   \item{1}{Signal exceeds baseline (peak) — potential G nucleotide.}
#'   \item{-1}{Signal falls below baseline (valley) — potential C or U
#'     nucleotide.}
#'   \item{0}{Signal within baseline range — homopolymer A region.}
#' }
#'
#' @seealso \code{\link{extract_tail_data_trainingset}} which calls this
#'   function,
#'   \code{\link{filter_signal_by_threshold}} for the production counterpart,
#'   \code{\link{substitute_gaps}} for the gap-filling post-processing step.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' filter_signal_by_threshold_trainingset(
#'   signal = tail_feature_list[[1]][["readname"]][[4]])
#'
#' }
filter_signal_by_threshold_trainingset <- function(signal) {

  #assertions
  if (missing(signal)) {
    stop(
      "Signal is missing. Please provide a valid signal argument [numeric vec].",
      call. = FALSE
    )
  }
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
  SD_threshold <- 3.5 # how many SD tresholds from avg signal pseudomove should be reported

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

  # this is to prevent from extracting poor quality ones
  adjusted_pseudomoves <- ninetails::substitute_gaps(adjusted_pseudomoves)

  #hardcode
  #substitute original terminal values by zeros, so the terminal chunks
  #would be excluded from the prefiltered dataset; this is to prevent
  #potential terminal segmentation errors to be taken into consideration.
  adjusted_pseudomoves[1:5] <- 0 #beginning
  adjusted_pseudomoves[
    (length(adjusted_pseudomoves) - 5):length(adjusted_pseudomoves)
  ] <- 0 # end

  return(adjusted_pseudomoves)
}


#' Extracts decoration-centered fragments of poly(A) tail signal along
#' with positional coordinates.
#'
#' Splits a single read's poly(A) tail signal into fixed-length (100
#' data-point) chunks, each centered on a potential non-A decoration
#' detected from the pseudomove vector. Chunks whose start position falls
#' before index 1 are left-padded with values randomly sampled from the
#' three most frequent signal values.
#'
#' @details
#' This training-set variant includes an additional \code{pseudomoves}
#' element in each chunk sublist, making the output suitable for
#' supervised training/validation data preparation. Pseudomoves for each
#' extracted chunk are recomputed by calling
#' \code{\link{filter_signal_by_threshold}} on the chunk sequence.
#'
#' The centering procedure:
#' \enumerate{
#'   \item Runs RLE on the pseudomove vector.
#'   \item Selects runs of length >= 4 with non-zero values (empirical
#'         modification threshold).
#'   \item Centers a 100-element window on the midpoint of each selected
#'         run.
#' }
#'
#' Chunk names follow the convention \code{<readname>_<index>}, where
#' \code{index} is the sequential position of the modification within
#' the read (numbered from the 3' end).
#'
#' @param readname Character string. Name (UUID) of the given read.
#'
#' @param tail_feature_list List object produced by
#'   \code{\link{create_tail_feature_list_trainingset}}.
#'
#' @return A named nested list where each element corresponds to one chunk
#'   and contains:
#' \describe{
#'   \item{chunk_sequence}{Numeric vector (length 100). The signal
#'     fragment centered on the potential modification.}
#'   \item{chunk_start_pos}{Integer. Start index of the chunk within the
#'     full tail signal (may be negative for left-padded chunks).}
#'   \item{chunk_end_pos}{Integer. End index of the chunk within the full
#'     tail signal.}
#'   \item{pseudomoves}{Numeric vector (length 100). Recomputed
#'     pseudomoves for the chunk. Coordinates are from the 3' end.}
#' }
#'
#' @seealso \code{\link{create_tail_chunk_list_trainingset}} for the
#'   parallel wrapper that calls this function,
#'   \code{\link{filter_signal_by_threshold}} for pseudomove recomputation
#'   on individual chunks,
#'   \code{\link{split_tail_centered}} for the production counterpart.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' split_tail_centered_trainingset(
#'   readname = "1234-anexample-r3adn4m3",
#'   tail_feature_list = tail_feature_list)
#'
#' }
split_tail_centered_trainingset <- function(readname, tail_feature_list) {

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
  pseudomoves <- tail_feature_list[[1]][[readname]][[4]]

  # mod-centered chunk extraction
  # recompute rle
  mod_rle <- rle(pseudomoves)
  # pseudomoves filtered by condition (potentially modified - empyrical!)
  condition <- mod_rle$lengths >= 4 & mod_rle$values
  #beginning positions of filtered pseudomoves which satisfy conditions
  first_filtered_positions <- cumsum(c(1, utils::head(mod_rle$lengths, -1)))[
    condition
  ]
  #length of pseudomoves satisfying condition
  filtered_length <- mod_rle$lengths[condition]
  #extracted coordinates (indices)
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

  # replace NAs with 3 most freq values (impute the data for gafs)
  most_freq_vals <- as.numeric(names(sort(table(signal), decreasing = TRUE)[
    1:3
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

  # creating final list output (retrieve coordinates as list_2 and list_3; retrieving pseudomoves as list_4):
  list_2 <- as.list(start_positions)
  #names(list_2) <- chunk_names
  list_3 <- as.list(end_positions)
  #names(list_3) <- chunk_names
  #list of pseudomoves -> recreated based on chopped signals
  list_4 <- lapply(list_1, function(i) filter_signal_by_threshold(i))

  out <- mget(ls(pattern = '^list.*\\d$')) %>%
    split(sub("_\\d+$", '', names(.))) %>%
    purrr::map(
      ~ purrr::transpose(purrr::set_names(
        .,
        c('chunk_sequence', 'chunk_start_pos', 'chunk_end_pos', 'pseudomoves')
      ))
    ) %>%
    purrr::flatten(.)

  return(out)
}

#' Extracts decoration-centered fragments of poly(A) tails for all reads
#' and appends positional data to a nested list.
#'
#' Parallel wrapper around \code{\link{split_tail_centered_trainingset}}.
#' For every read in the feature list, it extracts 100-element signal chunks
#' centered on potential non-A modifications and organises them in a nested
#' list keyed by read ID.
#'
#' @details
#' This training-set variant is intended for preparing training and
#' validation datasets. Each chunk sublist contains four fields:
#' \code{chunk_sequence}, \code{chunk_start_pos}, \code{chunk_end_pos},
#' and \code{pseudomoves} (see
#' \code{\link{split_tail_centered_trainingset}}).
#'
#' @param tail_feature_list List object produced by
#'   \code{\link{create_tail_feature_list_trainingset}}.
#'
#' @param num_cores Numeric \code{[1]}. Number of physical cores to use.
#'   Do not exceed 1 less than the number of cores at your disposal.
#'
#' @return A named nested list organised by read IDs, where each element
#'   is a list of chunk sublists as returned by
#'   \code{\link{split_tail_centered_trainingset}}.
#'
#' @seealso \code{\link{create_tail_feature_list_trainingset}} for the
#'   preceding pipeline step,
#'   \code{\link{split_tail_centered_trainingset}} for the per-read
#'   extraction logic,
#'   \code{\link{filter_nonA_chunks_trainingset}} for the next pipeline
#'   step,
#'   \code{\link{create_gaf_list}} for GAF conversion downstream.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' create_tail_chunk_list_trainingset(
#'   tail_feature_list = tail_feature_list,
#'   num_cores = 2)
#'
#' }
create_tail_chunk_list_trainingset <- function(tail_feature_list, num_cores) {

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

  #create empty list for extracted data
  tail_chunk_list = list()

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
        ninetails::split_tail_centered_trainingset(x, tail_feature_list)
      })
    }

  #close(pb)

  #rename first level of nested list accordingly
  names(tail_chunk_list) <- names(tail_feature_list[[1]])

  # Done comm
  cat(paste0('[', as.character(Sys.time()), '] ', 'Done!', '\n', sep = ''))

  return(tail_chunk_list)
}


#' Splits signal to overlapping fragments of equal length.
#'
#' Divides a single read's poly(A) tail signal into fixed-length,
#' overlapping segments for data augmentation during training-set
#' preparation. If the signal is not evenly divisible by the segment
#' length, trailing \code{NA} values are imputed with values randomly
#' sampled from the 10 most frequent signal values.
#'
#' @details
#' This function is used exclusively in the A-only training pipeline
#' (\code{\link{create_tail_chunk_list_A}}) to produce overlapping
#' windows from reads that do \emph{not} contain modification signals.
#' The overlap acts as data augmentation, multiplying the number of
#' training examples per read.
#'
#' @param readname Character string. Name (UUID) of the given ONT signal.
#'
#' @param tail_feature_list List object produced by
#'   \code{\link{create_tail_feature_list_A}}.
#'
#' @param segment Numeric \code{[1]}. Length (in data points) of each
#'   chunk to be created.
#'
#' @param overlap Numeric \code{[1]}. Length (in data points) of the
#'   overlap between consecutive chunks.
#'
#' @return A list of numeric vectors, each of length \code{segment},
#'   representing the overlapping signal fragments. Trailing \code{NA}
#'   values are replaced by imputed values.
#'
#' @seealso \code{\link{create_tail_chunk_list_A}} which calls this
#'   function,
#'   \code{\link{create_tail_feature_list_A}} for the preceding step,
#'   \code{\link{split_tail_centered_trainingset}} for the
#'   modification-centered splitting alternative.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' split_with_overlaps(
#'   readname = '12fdcb3-ewfd543-34552-1234ddta345',
#'   tail_feature_list = tail_feature_list,
#'   segment = 300,
#'   overlap = 100)
#'
#' }
split_with_overlaps <- function(readname,
                                tail_feature_list,
                                segment,
                                overlap) {

  #extract required data
  signal <- tail_feature_list[[1]][[readname]][[2]]

  starts <- seq(1, length(signal), by = segment - overlap)
  ends <- starts + segment - 1
  most_freq_vals <- as.numeric(names(sort(table(signal), decreasing = TRUE)[
    1:10
  ]))

  split_signal <- lapply(1:length(starts), function(i) {
    signal[starts[i]:ends[i]]
  })

  # replace NAs with 3 most frequent values (randomly sampled)
  #this is to avoid an error generated by cut() wrapped in color_matrix in gaf creator function
  #if all values would be equal, so the breaks would not be unique
  result_split <- lapply(split_signal, function(n) {
    replace(n, is.na(n), sample(most_freq_vals, sum(is.na(n)), T))
  })

  return(result_split)
}

#' Extracts features of poly(A) tails containing only A nucleotides for
#' training-set preparation.
#'
#' Training-set variant of the feature extraction wrapper designed
#' exclusively for A-only signals. It processes all reads in parallel via
#' \code{\link{extract_tail_data_trainingset}}, then applies an inverse
#' filtering criterion: only reads whose pseudomove vectors do \emph{not}
#' contain consecutive non-zero runs of length >= 4 are retained. This
#' ensures the resulting dataset represents pure homopolymer A tails
#' without modification artefacts.
#'
#' @details
#' The inverse filtering uses a sliding-window approach
#' (\code{stats::embed}) to detect runs of >= 4 consecutive non-zero
#' pseudomoves. Reads that pass this filter (i.e. have no such runs)
#' are collected as the A-only reference set. Unlike
#' \code{\link{create_tail_feature_list_trainingset}}, the returned list
#' does not include \code{zeromoved_readnames} or
#' \code{nonpseudomoved_readnames} categories.
#'
#' @param nanopolish Character string. Full path of the \code{.tsv} file
#'   produced by \code{nanopolish polya}.
#'
#' @param sequencing_summary Character string. Full path of the \code{.txt}
#'   file with the sequencing summary.
#'
#' @param workspace Character string. Full path of the directory containing
#'   basecalled multi-Fast5 files.
#'
#' @param num_cores Numeric \code{[1]}. Number of physical cores to use.
#'   Do not exceed 1 less than the number of cores at your disposal.
#'
#' @param basecall_group Character string \code{["Basecall_1D_000"]}. Name of
#'   the level in the Fast5 file hierarchy from which the data should be
#'   extracted.
#'
#' @param pass_only Logical \code{[TRUE]}. If \code{TRUE}, only reads tagged
#'   by nanopolish as \code{"PASS"} are retained. Otherwise, reads tagged as
#'   \code{"PASS"} or \code{"SUFFCLIP"} are included.
#'
#' @return A named list with one element:
#' \describe{
#'   \item{tail_feature_list}{Named list of per-read feature lists (as
#'     returned by \code{\link{extract_tail_data_trainingset}}) containing
#'     only reads with pure A tails.}
#' }
#' Always assign this returned list to a variable; printing the full list
#' to the console may crash the R session.
#'
#' @seealso \code{\link{create_tail_feature_list_trainingset}} for the
#'   non-A variant,
#'   \code{\link{create_tail_chunk_list_A}} for the next pipeline step,
#'   \code{\link{extract_tail_data_trainingset}} for per-read extraction,
#'   \code{\link{prepare_trainingset}} for the top-level wrapper.
#'
#' @importFrom foreach %dopar%
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' create_tail_feature_list_A(
#'   nanopolish = '/path/to/file',
#'   sequencing_summary = '/path/to/file',
#'   workspace = '/path/to/guppy/workspace',
#'   num_cores = 10,
#'   basecall_group = 'Basecall_1D_000')
#'
#' }
create_tail_feature_list_A <- function(nanopolish,
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
        ninetails::extract_tail_data_trainingset(
          x,
          polya_summary,
          workspace,
          basecall_group
        )
      })
    }


  #label each signal according to corresponding read name to avoid confusion
  squiggle_names <- polya_summary$readname
  names(tail_features_list) <- squiggle_names

  # prevent from running on reads which fulfill the pseudomove condition
  # this allows to create a set of reads with singular pseudomoves/pseudomoves shorter than 4 subsequent occurrences.
  # under this condition, reads without signal distortions would be directed to the list of interest
  # this condition is most versatile, as it takes into consideration signals without moves, with misdiagnosed moves
  # and signals with singular artifacts from z-score algo; so the A-containing dataset would be robust
  find_nonzeros <- function(x, y) rowSums(stats::embed(x != 0, y))
  tail_features_list <- Filter(
    function(x) all(find_nonzeros(x$tail_pseudomoves, 4) < 4),
    tail_features_list
  )

  #create final output
  tail_feature_list <- list()

  tail_feature_list[["tail_feature_list"]] <- tail_features_list

  # Done comm
  cat(paste0('[', as.character(Sys.time()), '] ', 'Done!', '\n', sep = ''))

  return(tail_feature_list)
}


#' Creates list of tail chunks containing only A nucleotides.
#'
#' Extracts overlapping fragments of poly(A) tail signals that do
#' \emph{not} contain non-A nucleotides and organises them in a nested
#' list keyed by read ID.
#'
#' @details
#' This training-set variant is designed for A-only reference data
#' preparation. Unlike \code{\link{create_tail_chunk_list_trainingset}}
#' (which centres chunks on modifications), this function uses
#' \code{\link{split_with_overlaps}} with \code{segment = 100} and
#' \code{overlap = 50}, effectively performing data augmentation by
#' producing overlapping windows across the entire tail signal.
#'
#' @param tail_feature_list List object produced by
#'   \code{\link{create_tail_feature_list_A}}.
#'
#' @param num_cores Numeric \code{[1]}. Number of physical cores to use.
#'   Do not exceed 1 less than the number of cores at your disposal.
#'
#' @return A named nested list organised by read IDs, where each element
#'   is a list of numeric vectors (each of length 100) representing the
#'   overlapping signal fragments.
#'
#' @seealso \code{\link{create_tail_feature_list_A}} for the preceding
#'   pipeline step,
#'   \code{\link{split_with_overlaps}} for the per-read segmentation logic,
#'   \code{\link{create_gaf_list_A}} for the next pipeline step,
#'   \code{\link{create_tail_chunk_list_trainingset}} for the non-A variant.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' create_tail_chunk_list_A(
#'   tail_feature_list = tail_feature_list,
#'   num_cores = 2)
#'
#' }
#'
create_tail_chunk_list_A <- function(tail_feature_list, num_cores) {

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

  #create empty list for extracted data
  tail_chunk_list = list()

  # creating cluster for parallel computing
  my_cluster <- parallel::makeCluster(num_cores)
  on.exit(parallel::stopCluster(my_cluster))
  doSNOW::registerDoSNOW(my_cluster)
  `%dopar%` <- foreach::`%dopar%`
  `%do%` <- foreach::`%do%`
  mc_options <- list(preschedule = TRUE, set.seed = FALSE, cleanup = FALSE)

  # header for progress bar
  cat(paste('Creating A-exclusive tail segmentation data...', '\n', sep = ''))

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
        ninetails::split_with_overlaps(
          x,
          tail_feature_list,
          segment = 100,
          overlap = 50
        )
      })
    }


  #rename first level of nested list accordingly
  names(tail_chunk_list) <- names(tail_feature_list[[1]])

  return(tail_chunk_list)
}

#' Produces list of GAFs containing exclusively A-nucleotides for neural
#' network training.
#'
#' Converts A-only signal chunks into Gramian Angular Field (GAF)
#' matrices in parallel. Each chunk is transformed via
#' \code{\link{combine_gafs}} and the resulting matrices are collected in
#' a flat named list suitable for direct input to the CNN training
#' pipeline.
#'
#' @details
#' This training-set variant is designed to work on signal fragments
#' that contain only A residues (produced by
#' \code{\link{create_tail_chunk_list_A}}). The naming convention for
#' output elements is \code{<readname>_<chunk_index>}.
#'
#' @param tail_chunk_list List object produced by
#'   \code{\link{create_tail_chunk_list_A}}.
#'
#' @param num_cores Numeric \code{[1]}. Number of physical cores to use.
#'   Do not exceed 1 less than the number of cores at your disposal.
#'
#' @return A named list of GAF matrices, where each element corresponds
#'   to one signal chunk and is named \code{<readname>_<index>}. Always
#'   assign this returned list to a variable; printing the full list to
#'   the console may crash the R session.
#'
#' @seealso \code{\link{create_tail_chunk_list_A}} for the preceding
#'   pipeline step,
#'   \code{\link{combine_gafs}} for the GAF transformation,
#'   \code{\link{create_gaf_list}} for the production counterpart,
#'   \code{\link{prepare_trainingset}} for the top-level wrapper.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' create_gaf_list_A(
#'   tail_chunk_list = tail_chunk_list,
#'   num_cores = 10)
#'
#' }
create_gaf_list_A <- function(tail_chunk_list, num_cores) {

  # Assertions
  if (missing(num_cores)) {
    stop(
      "Number of declared cores is missing. Please provide a valid num_cores argument.",
      call. = FALSE
    )
  }

  if (missing(tail_chunk_list)) {
    stop(
      "List of tail chunks is missing. Please provide a valid chunk_list argument.",
      call. = FALSE
    )
  }

  assert_condition(
    is.numeric(num_cores),
    "Declared core number must be numeric. Please provide a valid argument."
  )

  #create empty list for extracted data
  gaf_list = list()

  # creating cluster for parallel computing
  my_cluster <- parallel::makeCluster(num_cores)
  on.exit(parallel::stopCluster(my_cluster))
  doSNOW::registerDoSNOW(my_cluster)
  `%dopar%` <- foreach::`%dopar%`
  `%do%` <- foreach::`%do%`
  mc_options <- list(preschedule = TRUE, set.seed = FALSE, cleanup = FALSE)

  # header for progress bar
  cat(paste('Computing gramian angular fields...', '\n', sep = ''))

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
      lapply(tail_chunk_list[[i]], function(x) ninetails::combine_gafs(x))
    }


  #naming chunks based on readnames & indices
  for (chunk in seq_along(tail_chunk_list)) {
    names(tail_chunk_list[[chunk]]) <- paste(
      names(tail_chunk_list)[chunk],
      seq_along(tail_chunk_list[[chunk]]),
      sep = "_"
    )
  }

  chunk_names <- gsub(".*\\.", "", names(rapply(tail_chunk_list, class)))
  names(gaf_list) <- chunk_names

  return(gaf_list)
}

#' Filters read chunks containing non-adenosine nucleotides of interest for
#' neural network training-set preparation.
#'
#' Designed for use with synthetic spike-in data containing a single type
#' of non-A residue (G, C, or U in the context of 3'-homopolymer A). The
#' function retains only chunks whose pseudomove vectors contain a
#' consecutive run of the specified \code{value} (peak or valley) of
#' length >= 4 and whose start position is non-negative.
#'
#' @details
#' The filtering takes advantage of the characteristic pseudomove
#' patterns produced by the z-score signal filter:
#' \itemize{
#'   \item G nucleotides produce a \strong{peak} (pseudomove = +1).
#'   \item C and U nucleotides produce a \strong{valley} (pseudomove = -1).
#' }
#'
#' The function does \emph{not} distinguish between C and U (both produce
#' valleys); classification is handled downstream by the CNN. Therefore
#' it is essential that training datasets for C and U differ in the
#' transcript bodies they map to and/or are not sequenced in a single run.
#'
#' \strong{Important:} this function is \emph{not} suitable for preparing
#' the A-only reference set. Use
#' \code{\link{create_tail_feature_list_A}} for that purpose.
#'
#' Before proceeding, visual inspection of at least some filtered signals
#' is advisable. It may be necessary to manually adjust hardcoded
#' parameters (contact the developer/maintainer for details).
#'
#' @param tail_chunk_list List object produced by
#'   \code{\link{create_tail_chunk_list_trainingset}}.
#'
#' @param value Numeric \code{[1]}. Controls whether valleys (C/U) or peaks
#'   (G) are retained:
#' \describe{
#'   \item{-1}{Retain chunks containing valleys (C or U nucleotides).}
#'   \item{1}{Retain chunks containing peaks (G nucleotide).}
#' }
#'
#' @param num_cores Numeric \code{[1]}. Number of physical cores to use.
#'   Do not exceed 1 less than the number of cores at your disposal.
#'
#' @return A named list of filtered chunk sublists, where names correspond
#'   to read IDs. Empty reads (no chunks passing the filter) are removed.
#'
#' @seealso \code{\link{create_tail_chunk_list_trainingset}} for the
#'   preceding pipeline step,
#'   \code{\link{create_gaf_list}} for GAF conversion downstream,
#'   \code{\link{prepare_trainingset}} for the top-level wrapper.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # filtering G residue:
#' filter_nonA_chunks_trainingset(
#'   tail_chunk_list = list_object_with_tail_chunks,
#'   value = 1, num_cores = 2)
#'
#' # filtering C/U (user must know which type of residue is dealt with):
#' filter_nonA_chunks_trainingset(
#'   tail_chunk_list = list_object_with_tail_chunks,
#'   value = -1, num_cores = 2)
#'
#' }
filter_nonA_chunks_trainingset <- function(tail_chunk_list, value, num_cores) {

  #assertions
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
  assert_condition(
    is.numeric(value),
    "Provided value must be numeric. Please provide a valid argument."
  )
  assert_condition(
    value %in% c(-1, 1),
    "Provided value must be either 1 or -1. Please provide a valid argument."
  )

  #empty list for output
  filtered_input <- list()

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
    'Filtering training nonA chunks...',
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

  #filter the potential nonA containing chunks for training set preparation
  filtered_input <- foreach::foreach(
    i = seq_along(tail_chunk_list),
    .inorder = TRUE,
    .errorhandling = 'pass',
    .options.snow = opts,
    .options.multicore = mc_options
  ) %dopar%
    {
      filtered_output <- Filter(
        function(x) {
          any(with(rle(x$pseudomoves), lengths[values == value] >= 4)) &
            x$chunk_start_pos >= 0
        },
        tail_chunk_list[[i]]
      )
      lapply(filtered_output, function(x) x)
    }


  #remove empty sublists
  filtered_input <- Filter(function(x) length(x) > 0, filtered_input)

  #restore names in the list
  chunknames <- purrr::map_depth(filtered_input, 1, names) %>%
    unlist(use.names = F)
  chunknames <- unique(gsub("\\_[0-9]*$", "", chunknames))
  names(filtered_input) <- chunknames

  # Done comm
  cat(paste0('[', as.character(Sys.time()), '] ', 'Done!', '\n', sep = ''))

  return(filtered_input)
}


#' Filters out signals of a given nucleotide type for neural network
#' training-set preparation.
#'
#' Top-level convenience wrapper that orchestrates the complete
#' training-set production pipeline for a single nucleotide category.
#' Depending on the selected \code{nucleotide}, it chains the appropriate
#' feature extraction, chunk splitting, filtering, and GAF creation
#' functions into a single call.
#'
#' @details
#' The internal pipeline differs by nucleotide:
#' \describe{
#'   \item{\code{"A"}}{Uses the A-only branch:
#'     \code{\link{create_tail_feature_list_A}} \eqn{\rightarrow}
#'     \code{\link{create_tail_chunk_list_A}} \eqn{\rightarrow}
#'     \code{\link{create_gaf_list_A}}.}
#'   \item{\code{"C"}}{Uses the non-A branch with \code{value = -1}
#'     (valley filtering):
#'     \code{\link{create_tail_feature_list_trainingset}} \eqn{\rightarrow}
#'     \code{\link{create_tail_chunk_list_trainingset}} \eqn{\rightarrow}
#'     \code{\link{filter_nonA_chunks_trainingset}} \eqn{\rightarrow}
#'     \code{\link{create_gaf_list}}.}
#'   \item{\code{"G"}}{Uses the non-A branch with \code{value = 1}
#'     (peak filtering).}
#'   \item{\code{"U"}}{Uses the non-A branch with \code{value = -1}
#'     (valley filtering), same filter direction as C.}
#' }
#'
#' @param nucleotide Character. One of \code{"A"}, \code{"C"}, \code{"G"},
#'   or \code{"U"}. Defines the type of signal filtering applied to
#'   produce training data for the desired nucleotide context.
#'
#' @param nanopolish Character string. Full path of the \code{.tsv} file
#'   produced by \code{nanopolish polya}.
#'
#' @param sequencing_summary Character string. Full path of the \code{.txt}
#'   file with the sequencing summary.
#'
#' @param workspace Character string. Full path of the directory containing
#'   basecalled multi-Fast5 files.
#'
#' @param num_cores Numeric \code{[1]}. Number of physical cores to use.
#'   Do not exceed 1 less than the number of cores at your disposal.
#'
#' @param basecall_group Character string \code{["Basecall_1D_000"]}. Name of
#'   the level in the Fast5 file hierarchy from which the data should be
#'   extracted.
#'
#' @param pass_only Logical \code{[TRUE]}. If \code{TRUE}, only reads tagged
#'   by nanopolish as \code{"PASS"} are retained. Otherwise, reads tagged as
#'   \code{"PASS"} or \code{"SUFFCLIP"} are included.
#'
#' @return A named list of GAF matrices organised by
#'   \code{<read_ID>_<index>}. Always assign this returned list to a
#'   variable; printing the full list to the console may crash the R session.
#'
#' @seealso \code{\link{create_tail_feature_list_trainingset}},
#'   \code{\link{create_tail_feature_list_A}},
#'   \code{\link{create_tail_chunk_list_trainingset}},
#'   \code{\link{create_tail_chunk_list_A}},
#'   \code{\link{filter_nonA_chunks_trainingset}},
#'   \code{\link{create_gaf_list}},
#'   \code{\link{create_gaf_list_A}}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' prepare_trainingset(
#'   nucleotide = "A",
#'   nanopolish = '/path/to/file',
#'   sequencing_summary = '/path/to/file',
#'   workspace = '/path/to/guppy/workspace',
#'   num_cores = 10,
#'   basecall_group = 'Basecall_1D_000',
#'   pass_only = TRUE)
#'
#' }
#'
prepare_trainingset <- function(nucleotide,
                                nanopolish,
                                sequencing_summary,
                                workspace,
                                num_cores = 1,
                                basecall_group = "Basecall_1D_000",
                                pass_only = TRUE) {

  # nucleotide selection
  if (nucleotide == "A") {
    # process the A-containing data; split tails with overlaps
    tail_feature_list <- ninetails::create_tail_feature_list_A(
      nanopolish,
      sequencing_summary,
      workspace,
      num_cores,
      basecall_group,
      pass_only
    )
    tail_chunk_list <- ninetails::create_tail_chunk_list_A(
      tail_feature_list,
      num_cores
    )
    gafs_list <- ninetails::create_gaf_list_A(tail_chunk_list, num_cores)
  } else if (nucleotide == "C") {
    # process the C-containing data to filter only internal, C-containing chunks
    tail_feature_list <- ninetails::create_tail_feature_list_trainingset(
      nanopolish,
      sequencing_summary,
      workspace,
      num_cores,
      basecall_group,
      pass_only
    )
    tail_chunk_list <- ninetails::create_tail_chunk_list_trainingset(
      tail_feature_list,
      num_cores
    )
    filtered_chunk_list <- ninetails::filter_nonA_chunks_trainingset(
      tail_chunk_list,
      value = -1,
      num_cores
    )
    gafs_list <- ninetails::create_gaf_list(filtered_chunk_list, num_cores)
  } else if (nucleotide == "G") {
    # process the G-containing data to filter only internal, G-containing chunks
    tail_feature_list <- ninetails::create_tail_feature_list_trainingset(
      nanopolish,
      sequencing_summary,
      workspace,
      num_cores,
      basecall_group,
      pass_only
    )
    tail_chunk_list <- ninetails::create_tail_chunk_list_trainingset(
      tail_feature_list,
      num_cores
    )
    filtered_chunk_list <- ninetails::filter_nonA_chunks_trainingset(
      tail_chunk_list,
      value = 1,
      num_cores
    )
    gafs_list <- ninetails::create_gaf_list(filtered_chunk_list, num_cores)
  } else if (nucleotide == "U") {
    # process the U-containing data to filter only internal, U-containing chunks
    tail_feature_list <- ninetails::create_tail_feature_list_trainingset(
      nanopolish,
      sequencing_summary,
      workspace,
      num_cores,
      basecall_group,
      pass_only
    )
    tail_chunk_list <- ninetails::create_tail_chunk_list_trainingset(
      tail_feature_list,
      num_cores
    )
    filtered_chunk_list <- ninetails::filter_nonA_chunks_trainingset(
      tail_chunk_list,
      value = -1,
      num_cores
    )
    gafs_list <- ninetails::create_gaf_list(filtered_chunk_list, num_cores)
  } else {
    stop(
      "Wrong nucleotide selected. Please provide valid nucleotide argument.",
      call. = FALSE
    )
  }

  return(gafs_list)
}
