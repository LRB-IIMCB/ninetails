#' Extracts tail features of single RNA read from respective basecalled
#' multi-fast5 file.
#'
#' This is the version of the function useful to produce the training set
#' for CNN. Slightly different from the original one.
#'
#' @param readname character string. Name of the given read within the
#' analyzed dataset.
#'
#' @param polya_summary character string. The table containing data extracted
#' from nanopolish & sequencing summary (using extract_polya_data() function.
#'
#' @param workspace character string. Full path of the directory to search
#' the basecalled fast5 files in. The Fast5 files have to be multi-fast5 file.
#'
#' @param basecall_group character string ["Basecall_1D_000"]. Name of the
#' level in the fast5 file hierarchy from which the data should be extracted.
#'
#' @return A list containing read information organized by the read ID
#' is returned. Always assign this returned list to a variable, otherwise
#' the long list will be printed to the console, which may crash your R session.
#'
#' @export
#'
#' @examples
#'\dontrun{
#'
#' extract_tail_data_trainingset(readname = 'abc123de-fg45-6789-0987-6543hijk2109',
#'                               polya_summary = polya_summary_table,
#'                               workspace = '/path/to/folder/containing/multifast5s',
#'                               basecall_group = 'Basecall_1D_000')
#'
#'}
#'
#'

extract_tail_data_trainingset <- function(readname,
                                          polya_summary,
                                          workspace,
                                          basecall_group){

  #Assertions
  if (missing(readname)) {
    stop("Readname is missing. Please provide a valid readname argument.", call. =FALSE)
  }

  if (missing(workspace)) {
    stop("Directory with basecalled fast5s is missing. Please provide a valid workspace argument.", call. =FALSE)
  }

  if (missing(basecall_group)) {
    stop("Basecall group is missing. Please provide a valid basecall_group argument.", call. =FALSE)
  }

  if (missing(polya_summary)) {
    stop("Polya_summary is missing. Please provide a valid polya_summary argument.", call. =FALSE)
  }

  if (!is.data.frame(polya_summary) || nrow(polya_summary) == 0) {
    stop("Empty data frame provided as an input (polya_summary). Please provide valid input")
  }

  assertthat::assert_that(!is.na(workspace) && nchar(workspace) > 0,
                          msg = "Empty string provided as an input. Please provide a valid path to basecalled fast5 files.")

  assertthat::assert_that(is.character(workspace),
                          msg = paste0("Path to basecalled fast5 files is not a character string. Please provide a valid path to basecalled fast5 files."))
  assertthat::assert_that(is.character(readname),
                          msg = paste0("Given readname is not a character string. Please provide a valid readname argument."))

  # Extract data from fast5 file
  fast5_filenames <- polya_summary$filename
  fast5_readname <- paste0("read_",readname) # suffix "read_"
  fast5_file <- fast5_filenames[readname] # particular file name (eg. FAO12345.fast5)
  fast5_file_path <-file.path(workspace, fast5_file) #path to browsed fast5 file

  signal <- rhdf5::h5read(file.path(fast5_file_path),paste0(fast5_readname,"/Raw/Signal"))
  signal <- as.vector(signal) #vectorized sig (natively: array)

  # retrieving moves
  BaseCalled_template <- rhdf5::h5read(fast5_file_path,paste0(fast5_readname,"/Analyses/", basecall_group, "/BaseCalled_template"))
  move <- BaseCalled_template$Move #how the basecaller "moves" through the called sequence, and allows for a mapping from basecall to raw data
  move <- as.numeric(move) # change data type as moves are originally stored as raw
  #move <- factor(move, levels=c(0,1))

  #read parameter stored in channel_id group
  channel_id <- rhdf5::h5readAttributes(fast5_file_path,paste0(fast5_readname,"/channel_id")) # parent dir for attributes (within fast5 file)
  sampling_rate <- channel_id$sampling_rate # number of data points collected per second

  #read parameters (attrs) stored in basecall_1d_template
  basecall_1d_template <- rhdf5::h5readAttributes(fast5_file_path,paste0(fast5_readname,"/Analyses/Basecall_1D_000/Summary/basecall_1d_template")) # parent dir for attributes (within fast5)
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
  polya_end_position <- transcript_start_position -1

  # handle move data
  #this is to expand moves along raw signal
  moves_sample_wise_vector <- c(rep(move, each=stride), rep(NA, length(signal) - number_of_events))
  moves_tail_range <- moves_sample_wise_vector[polya_start_position:polya_end_position]

  # extract polya tail region of the signal
  # signal is winsorized here!
  signal <- ninetails::winsorize_signal(signal[polya_start_position:polya_end_position])

  # downsample (interpolate signal and moves) so the computation would be faster/easier to handle.
  signal <- round(stats::approx(signal, method = "linear", n=ceiling(0.2 * length(signal)))[[2]], digits=0)
  moves_tail_range <- round(stats::approx(moves_tail_range, method = "linear", n=ceiling(0.2 * length(moves_tail_range)))[[2]], digits=0)

  # filter signal to find local minima & maxima corresponding to potential C, G, U modifications
  pseudomoves <- filter_signal_by_threshold_trainingset(signal)

  extracted_data_single_list = list() # creating empty list for the extracted fast5 data

  extracted_data_single_list[["fast5_filename"]] <- polya_summary$filename[read_idx]
  extracted_data_single_list[["tail_signal"]] <- signal
  extracted_data_single_list[["tail_moves"]] <- moves_tail_range
  extracted_data_single_list[["tail_pseudomoves"]] <- pseudomoves

  return(extracted_data_single_list)
}


#' Extracts features of polyA tails of ONT RNA reads required for finding
#' non-A nucleotides within the given tails.
#'
#' This is the version of the function useful to produce the training set
#' for CNN. Slightly different from the original one.
#'
#' @param nanopolish character string. Full path of the .tsv file produced
#' by nanopolish polya function.
#'
#' @param sequencing_summary character string. Full path of the .txt file
#' with sequencing summary.
#'
#' @param workspace character string. Full path of the directory to search the
#' basecalled fast5 files in. The Fast5 files have to be multi-fast5 file.
#'
#' @param num_cores numeric [1]. Number of physical cores to use in processing
#' the data. Do not exceed 1 less than the number of cores at your disposal.
#'
#' @param basecall_group character string ["Basecall_1D_000"]. Name of the
#' level in the Fast5 file hierarchy from which the data should be extracted.
#'
#' @param pass_only logical [TRUE/FALSE]. If TRUE, only reads tagged by
#' nanopolish as "PASS" would be taken into consideration. Otherwise, reads
#' tagged as "PASS" & "SUFFCLIP" will be taken into account in analysis.
#' As a default, "TRUE" value is set.
#'
#' @return A list containing read information organized by the read ID
#' is returned. Always assign this returned list to a variable, otherwise
#' the long list will be printed to the console, which may crash your R session.
#'
#' @importFrom foreach %dopar%
#'
#' @export
#'
#' @examples
#'\dontrun{
#'
#' create_tail_feature_list_trainingset(nanopolish = '/path/to/file',
#'                                      sequencing_summary = '/path/to/file',
#'                                      workspace = '/path/to/guppy/workspace',
#'                                      num_cores = 10,
#'                                      basecall_group = 'Basecall_1D_000')
#'
#'}
#'
create_tail_feature_list_trainingset <- function(nanopolish,
                                                 sequencing_summary,
                                                 workspace,
                                                 num_cores,
                                                 basecall_group,
                                                 pass_only=TRUE){

  # variable binding (suppressing R CMD check from throwing an error)
  i <- NULL

  # Assertions
  if (missing(num_cores)) {
    stop("Number of declared cores is missing. Please provide a valid num_cores argument.", call. =FALSE)
  }

  if (missing(basecall_group)) {
    stop("Basecall group is missing. Please provide a valid basecall_group argument.", call. =FALSE)
  }

  if (missing(workspace)) {
    stop("Directory with basecalled fast5s (guppy workspace) is missing. Please provide a valid workspace argument.", call. =FALSE)
  }

  assertthat::assert_that(is.numeric(num_cores),
                          msg=paste0("Declared core number must be numeric. Please provide a valid argument."))

  # Extracting and processing polya & sequencing summary data
  polya_summary <- ninetails::extract_polya_data(nanopolish, sequencing_summary, pass_only)

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
  cat(paste0('[', as.character(Sys.time()), '] ','Extracting features of provided reads...', '\n', sep=''))

  # progress bar
  pb <- utils::txtProgressBar(min = 0,
                              max = length(polya_summary$filename),
                              style = 3,
                              width = 50,
                              char = "=",
                              file = stderr())
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)


  # parallel extraction
  tail_features_list <- foreach::foreach(i = seq_along(polya_summary$readname),
                                         .combine = c,
                                         .inorder = TRUE,
                                         .errorhandling = 'pass',
                                         .options.snow = opts,
                                         .options.multicore = mc_options) %dopar% {lapply(polya_summary$readname[i], function(x) ninetails::extract_tail_data_trainingset(x,polya_summary,workspace,basecall_group))}

  #close(pb)

  #label each signal according to corresponding read name to avoid confusion
  squiggle_names <- polya_summary$readname
  names(tail_features_list) <- squiggle_names

  # remove reads with only zero moved tails
  tail_features_list <- Filter(function(x) sum(x$tail_moves) !=0, tail_features_list)
  zeromoved_readnames <- squiggle_names[!(squiggle_names %in% names(tail_features_list))]

  # prevent from running on reads which do not fulfill the pseudomove condition
  tail_features_list <- Filter(function(x)any(with(rle(x$tail_pseudomoves), lengths[values!=0]>=4)), tail_features_list)

  # reads discarded because of unmet pseudomove condition
  #in this reads reported pseudomove chain is too short to be considered as potential modification
  nonpseudomoved_readnames <- squiggle_names[!(squiggle_names %in% c(zeromoved_readnames, names(tail_features_list)))]

  #create final output
  tail_feature_list <- list()

  tail_feature_list[["tail_feature_list"]] <- tail_features_list
  tail_feature_list[["zeromoved_readnames"]] <- zeromoved_readnames
  tail_feature_list[["nonpseudomoved_readnames"]] <- nonpseudomoved_readnames

  # Done comm
  cat(paste0('[', as.character(Sys.time()), '] ','Done!', '\n', sep=''))

  return(tail_feature_list)
}


#' Detection of outliers (peaks & valleys) in ONT signal using z-scores.
#'
#' This is the version of the function useful to produce the training set
#' for CNN. Slightly different from the original one.
#'
#' The function created based on the following source:
#' Brakel, J.P.G. van (2014). "Robust peak detection algorithm using z-scores".
#' Stack Overflow. Available at:
#' https://stackoverflow.com/questions/22583391/peak-signal-detection-in-realtime-timeseries-data/22640362#22640362
#' (version: 2020-11-08).
#'
#'
#' @param signal numeric vector. An ONT read fragment corresponding to the tail
#' region of the read of interest as delimited by nanopolish polya function (the
#' fragments are stored in tail_feature_list[[1]] produced by the
#' create_tail_feature_list() function.
#'
#' @return numeric vector of "pseudomoves" corresponding to the analyzed tail
#' region; containing values ranging from -1 to 1.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' filter_signal_by_threshold(signal=tail_feature_list[[1]][["readname"]][[4]])
#'
#'}
filter_signal_by_threshold_trainingset <- function(signal) {

  # variable binding
  x <- sd <- baseline <- std_cutoff <- NULL

  #assertions
  if (missing(signal)) {
    stop("Signal is missing. Please provide a valid signal argument [numeric vec].", call. =FALSE)
  }
  # reproducibility
  set.seed(123)

  # calibrate the algo on the sampled vals
  start_vals <- signal[1:10]
  most_freq_vals <- as.numeric(names(sort(table(signal),decreasing=TRUE)[1:20]))
  adjusted_signal <- c(sample(c(most_freq_vals, start_vals), 100, replace=TRUE), signal)

  # Empyrical parameters:
  adaptive_sampling_window <- 100 # datapoints window for adjusting algo
  SD_threshold <- 3.5 # how many SD tresholds from avg signal pseudomove should be reported

  pseudomoves <- rep(0,length(adjusted_signal))
  filtered_signal <- adjusted_signal

  baseline[adaptive_sampling_window] <- mean(adjusted_signal[1:adaptive_sampling_window], na.rm=TRUE)
  std_cutoff[adaptive_sampling_window] <- stats::sd(adjusted_signal[1:adaptive_sampling_window], na.rm=TRUE)
  for (i in (adaptive_sampling_window+1):length(adjusted_signal)){
    if (abs(adjusted_signal[i] - baseline[i-1]) > SD_threshold*std_cutoff[i-1]) {
      if (adjusted_signal[i] > baseline[i-1]) {
        pseudomoves[i] <- 1 #if they go up
      } else {
        pseudomoves[i] <- -1 # if they go down
      }
      filtered_signal[i] <- filtered_signal[i-1] #update
    } else {
      pseudomoves[i] <- 0 # uniform distr
      filtered_signal[i] <- adjusted_signal[i] #update
    }
    baseline[i] <- mean(filtered_signal[(i-adaptive_sampling_window):i], na.rm=TRUE)
    std_cutoff[i] <- stats::sd(filtered_signal[(i-adaptive_sampling_window):i], na.rm=TRUE)
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
  adjusted_pseudomoves[(length(adjusted_pseudomoves)-5):length(adjusted_pseudomoves)] <- 0 # end

  return(adjusted_pseudomoves)
}


#' Extracts fragments of poly(A) tail signal containing potential modifications
#' along with its delimitation (positional indices; coordinates) within the tail.
#'
#' This version of the function contains an additional category (pseudomoves)
#' within the each chunk sublist, therefore is intended for preparation the
#' training/validation data.
#'
#' @param readname character string. Name of the given read within the
#' analyzed dataset.
#'
#' @param tail_feature_list list object produced by create_tail_feature_list
#' function.
#'
#' @return a list object (nested) containing 4 categories: resulting fragments
#' (chunk_sequence), start coordinate of given fragment (chunk_start_pos)
#' its' end coordinate (chunk_end_pos) and "pseudomoves" (pseudomoves)
#' corresponding to the given signal chunk arranged by the signal ID and
#' positional indices (WARNING! from 3' end!).
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' split_tail_centered_trainingset(readname= "1234-anexample-r3adn4m3",
#'                                 tail_feature_list = tail_feature_list)
#'
#'}
split_tail_centered_trainingset <- function(readname,
                                            tail_feature_list) {
  #variable binding
  tail <- . <- NULL

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

  #extract required data
  signal <- tail_feature_list[[1]][[readname]][[2]]
  pseudomoves <- tail_feature_list[[1]][[readname]][[4]]

  # mod-centered chunk extraction
  # recompute rle
  mod_rle <- rle(pseudomoves)
  # pseudomoves filtered by condition (potentially modified - empyrical!)
  condition <- mod_rle$lengths >= 4 & mod_rle$values
  #beginning positions of filtered pseudomoves which satisfy conditions
  first_filtered_positions <- cumsum(c(1, utils::head(mod_rle$lengths,-1)))[condition]
  #length of pseudomoves satisfying condition
  filtered_length <- mod_rle$lengths[condition]
  #extracted coordinates (indices)
  start_positions <-  first_filtered_positions + floor(filtered_length/2) - 50
  end_positions <- start_positions + 99

  # extract signal chunks centered on potential modification
  list_1 <- lapply(1:length(start_positions), function(i) if(start_positions[i] >0) signal[start_positions[i]:end_positions[i]]
                   else c(rep(NA, abs(start_positions[i]-1)), signal[1:end_positions[i]]))

  # replace NAs with 3 most freq values (impute the data for gafs)
  most_freq_vals <- as.numeric(names(sort(table(signal),decreasing=TRUE)[1:3]))
  list_1 <- lapply(list_1, function(n) replace(n, is.na(n), sample(most_freq_vals,sum(is.na(n)),TRUE)))

  #add indices
  chunks_indices <- c(1:length(start_positions))

  # naming chunks based on names & indices
  chunk_names <- paste0(rep(readname, length(list_1)), '_', unlist(chunks_indices))
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
    purrr::map(~purrr::transpose(purrr::set_names(.,c('chunk_sequence', 'chunk_start_pos', 'chunk_end_pos', 'pseudomoves')))) %>%
    purrr::flatten(.)

  return(out)
}

#' Extracts fragments of poly(A) tails of ONT RNA reads containing non-A
#' nucleotides along their coordinates & appends the data to the nested
#' list organized by read IDs.
#'
#' This version of the function is intended to be used to produce training and
#' validation datasets.
#'
#' @param tail_feature_list list object produced by create_tail_feature_list
#' function.
#'
#' @param num_cores numeric [1]. Number of physical cores to use in processing
#' the data. Do not exceed 1 less than the number of cores at your disposal.
#'
#' @return a list object 9nested) containing the segmented tail data  (chunks,
#' coordinates) organized by the read IDs.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' create_tail_chunk_list_trainingset(tail_feature_list = tail_feature_list,
#'                                    num_cores = 2)
#'}
create_tail_chunk_list_trainingset <- function(tail_feature_list,
                                               num_cores){

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
  `%do%` <- foreach::`%do%`
  mc_options <- list(preschedule = TRUE, set.seed = FALSE, cleanup = FALSE)

  # header for progress bar
  cat(paste0('[', as.character(Sys.time()), '] ','Creating tail segmentation data...', '\n', sep=''))

  # progress bar
  pb <- utils::txtProgressBar(min = 0,
                              max = length(tail_feature_list[[1]]),
                              style = 3,
                              width = 50,
                              char = "=",
                              file = stderr())
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  #create empty list for extracted data
  tail_chunk_list = list()

  #parallel extraction
  tail_chunk_list <- foreach::foreach(i = seq_along(tail_feature_list[[1]]),
                                      .combine = c,
                                      .inorder = TRUE,
                                      .errorhandling = 'pass',
                                      .options.snow = opts,
                                      .options.multicore = mc_options) %dopar% {
                                        lapply(names(tail_feature_list[[1]][i]), function(x) ninetails::split_tail_centered_trainingset(x,tail_feature_list))
                                      }

  #close(pb)

  #rename first level of nested list accordingly
  names(tail_chunk_list) <- names(tail_feature_list[[1]])

  # Done comm
  cat(paste0('[', as.character(Sys.time()), '] ','Done!', '\n', sep=''))

  return(tail_chunk_list)
}


#' Splits signal to overlapping fragments of equal length.
#'
#' In case if the signal is not completely divisible by given segment length,
#' the function fills missing data (NAs) with 10 most frequent values from the
#' entire signal (randomly sampled).
#'
#' @param readname character string. Name of the given ONT signal.
#'
#' @param tail_feature_list list object produced by create_tail_feature_list
#' function.
#'
#' @param segment numeric [1]. Length of the chunk(s) to be created.
#'
#' @param overlap numeric [1]. Length of the overlap between the chunks.
#'
#' @return a list object containing split chunks of signal.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' split_with_overlaps <- function(signal='12fdcb3-ewfd543-34552-1234ddta345',
#'                                 segment = 300, overlap 100)
#'
#' }
split_with_overlaps <- function(readname,
                                tail_feature_list,
                                segment,
                                overlap){

  #extract required data
  signal <- tail_feature_list[[1]][[readname]][[2]]

  starts <- seq(1, length(signal), by=segment-overlap)
  ends   <- starts + segment - 1
  most_freq_vals <- as.numeric(names(sort(table(signal),decreasing=TRUE)[1:10]))

  split_signal <- lapply(1:length(starts), function(i) signal[starts[i]:ends[i]])

  # replace NAs with 3 most frequent values (randomly sampled)
  #this is to avoid an error generated by cut() wrapped in color_matrix in gaf creator function
  #if all values would be equal, so the breaks would not be unique
  result_split <- lapply(split_signal, function(n) replace(n, is.na(n), sample(most_freq_vals,sum(is.na(n)),T)))

  return(result_split)
}

#' Extracts features of poly(A) tails of ONT RNA reads required for finding
#' non-A nucleotides within the given tails.
#'
#' This is the version of the function useful to produce the training set
#' for CNN. Slightly different from the original one. This version allows
#' to prepare the dataset of only A containing signals.
#'
#'
#' @param nanopolish character string. Full path of the .tsv file produced
#' by nanopolish polya function.
#'
#' @param sequencing_summary character string. Full path of the .txt file
#' with sequencing summary.
#'
#' @param workspace character string. Full path of the directory to search the
#' basecalled fast5 files in. The Fast5 files have to be multi-fast5 file.
#'
#' @param num_cores numeric [1]. Number of physical cores to use in processing
#' the data. Do not exceed 1 less than the number of cores at your disposal.
#'
#' @param basecall_group character string ["Basecall_1D_000"]. Name of the
#' level in the Fast5 file hierarchy from which the data should be extracted.
#'
#' @param pass_only logical [TRUE/FALSE]. If TRUE, only reads tagged by
#' nanopolish as "PASS" would be taken into consideration. Otherwise, reads
#' tagged as "PASS" & "SUFFCLIP" will be taken into account in analysis.
#' As a default, "TRUE" value is set.
#'
#' @return A list containing read information organized by the read ID
#' is returned. Always assign this returned list to a variable, otherwise
#' the long list will be printed to the console, which may crash your R session.
#'
#' @importFrom foreach %dopar%
#'
#' @export
#'
#' @examples
#'\dontrun{
#'
#' create_tail_feature_list_A(nanopolish = '/path/to/file',
#'                                      sequencing_summary = '/path/to/file',
#'                                      workspace = '/path/to/guppy/workspace',
#'                                      num_cores = 10,
#'                                      basecall_group = 'Basecall_1D_000')
#'}
create_tail_feature_list_A <- function(nanopolish,
                                       sequencing_summary,
                                       workspace,
                                       num_cores,
                                       basecall_group,
                                       pass_only=TRUE){

  # variable binding (suppressing R CMD check from throwing an error)
  i <- NULL

  # Assertions
  if (missing(num_cores)) {
    stop("Number of declared cores is missing. Please provide a valid num_cores argument.", call. =FALSE)
  }

  if (missing(basecall_group)) {
    stop("Basecall group is missing. Please provide a valid basecall_group argument.", call. =FALSE)
  }

  if (missing(workspace)) {
    stop("Directory with basecalled fast5s (guppy workspace) is missing. Please provide a valid workspace argument.", call. =FALSE)
  }

  assertthat::assert_that(is.numeric(num_cores),
                          msg=paste0("Declared core number must be numeric. Please provide a valid argument."))

  # Extracting and processing polya & sequencing summary data
  polya_summary <- ninetails::extract_polya_data(nanopolish, sequencing_summary, pass_only)

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
  cat(paste0('[', as.character(Sys.time()), '] ','Extracting features of provided reads...', '\n', sep=''))

  # progress bar
  pb <- utils::txtProgressBar(min = 0,
                              max = length(polya_summary$filename),
                              style = 3,
                              width = 50,
                              char = "=",
                              file = stderr())
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)


  # parallel extraction
  tail_features_list <- foreach::foreach(i = seq_along(polya_summary$readname),
                                         .combine = c,
                                         .inorder = TRUE,
                                         .errorhandling = 'pass',
                                         .options.snow = opts,
                                         .options.multicore = mc_options) %dopar% {lapply(polya_summary$readname[i], function(x) ninetails::extract_tail_data_trainingset(x,polya_summary,workspace,basecall_group))}

  #close(pb)

  #label each signal according to corresponding read name to avoid confusion
  squiggle_names <- polya_summary$readname
  names(tail_features_list) <- squiggle_names

  # prevent from running on reads which fulfill the pseudomove condition
  # this allows to create a set of reads with singular pseudomoves/pseudomoves shorter than 4 subsequent occurrences.
  # under this condition, reads without signal distortions would be directed to the list of interest
  # this condition is most versatile, as it takes into consideration signals without moves, with misdiagnosed moves
  # and signals with singular artifacts from z-score algo; so the A-containing dataset would be robust
  find_nonzeros <- function (x, y) rowSums(stats::embed(x != 0, y))
  tail_features_list <- Filter(function(x) all(find_nonzeros(x$tail_pseudomoves, 4) < 4), tail_features_list)

  #create final output
  tail_feature_list <- list()

  tail_feature_list[["tail_feature_list"]] <- tail_features_list

  # Done comm
  cat(paste0('[', as.character(Sys.time()), '] ','Done!', '\n', sep=''))

  return(tail_feature_list)
}


#' Creates list of tail chunks containing As exclusively.
#'
#' Extracts fragments of polyA tails of ONT RNA reads containing only A
#' nucleotides along their coordinates & appends the data to the nested
#' list organized by read IDs.
#'
#' This version of the function is intended to be used to produce training and
#' validation datasets. It is designed to filter signal fragments that DO NOT
#' contain potential non-A nucleotides.
#'
#' This function segments the signals with overlaps, making the data multiplied (data augmentation).
#'
#' @param tail_feature_list list object produced by create_tail_feature_list
#' function.
#'
#' @param num_cores numeric [1]. Number of physical cores to use in processing
#' the data. Do not exceed 1 less than the number of cores at your disposal.
#'
#' @return a list object (nested) containing the segmented tail data  (chunks,
#' coordinates) organized by the read IDs.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' create_tail_chunk_list_A(tail_feature_list = tail_feature_list,
#'                          num_cores = 2)
#'
#'}
#'
create_tail_chunk_list_A <- function(tail_feature_list,
                                     num_cores){
  #variable binding
  i <- NULL

  # initial assertions
  if (missing(num_cores)) {
    stop("Number of declared cores is missing. Please provide a valid num_cores argument.", call. =FALSE)
  }

  if (missing(tail_feature_list)) {
    stop("List of features is missing. Please provide a valid tail_feature_list argument.", call. =FALSE)
  }

  assertthat::assert_that(is.numeric(num_cores),
                          msg = paste("Declared core number must be numeric. Please provide a valid argument."))
  assertthat::assert_that(is.list(tail_feature_list),
                          msg = paste("Given tail_feature_list is not a list (class). Please provide valid file format."))

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
  cat(paste('Creating A-exclusive tail segmentation data...', '\n', sep=''))

  # progress bar
  pb <- utils::txtProgressBar(min = 0,
                              max = length(tail_feature_list[[1]]),
                              style = 3,
                              width = 50,
                              char = "=",
                              file = stderr())
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)


  #parallel extraction
  tail_chunk_list <- foreach::foreach(i = seq_along(tail_feature_list[[1]]),
                                      .combine = c,
                                      .inorder = TRUE,
                                      .errorhandling = 'pass',
                                      .options.snow = opts,
                                      .options.multicore = mc_options) %dopar% {
                                        lapply(names(tail_feature_list[[1]][i]), function(x) ninetails::split_with_overlaps(x,
                                                                                                                            tail_feature_list,
                                                                                                                            segment = 100,
                                                                                                                            overlap = 50))
                                      }

  #close(pb)

  #rename first level of nested list accordingly
  names(tail_chunk_list) <- names(tail_feature_list[[1]])

  return(tail_chunk_list)
}

#' Produces list of GAFs containing exclusively A-nucleotides for neural net
#' training.
#'
#' This version of the function is intended to be used to produce training and
#' validation datasets. It is designed to work on signal fragments that
#' contain only A residues.
#'
#' @param tail_chunk_list character string. The list object produced
#' by create_chunk_list function.
#'
#' @param num_cores numeric [1]. Number of physical cores to use in processing
#' the data. Do not exceed 1 less than the number of cores at your disposal.
#'
#' @return A list of gaf matrices organized by the read ID_index
#' is returned. Always assign this returned list to a variable, otherwise
#' the long list will be printed to the console, which may crash your R session.
#'
#' @export
#'
#' @examples
#'\dontrun{
#'
#' create_gaf_list_A(feature_list = tail_chunk_list, num_cores = 10)
#'}
create_gaf_list_A <- function(tail_chunk_list, num_cores){

  #variable binding
  i <- NULL

  # Assertions
  if (missing(num_cores)) {
    stop("Number of declared cores is missing. Please provide a valid num_cores argument.", call. =FALSE)
  }

  if (missing(tail_chunk_list)) {
    stop("List of tail chunks is missing. Please provide a valid chunk_list argument.", call. =FALSE)
  }

  assertthat::assert_that(is.numeric(num_cores),
                          msg=paste("Declared core number must be numeric. Please provide a valid argument."))

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
  cat(paste('Computing gramian angular fields...', '\n', sep=''))

  # progress bar
  pb <- utils::txtProgressBar(min = 0,
                              max = length(tail_chunk_list),
                              style = 3,
                              width = 50,
                              char = "=",
                              file = stderr())
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)


  gaf_list <- foreach::foreach(i = seq_along(tail_chunk_list),
                                .combine = c,
                                .inorder = TRUE,
                                .errorhandling = 'pass',
                                .options.snow = opts,
                               .options.multicore = mc_options) %dopar% {lapply(tail_chunk_list[[i]], function(x) ninetails::combine_gafs(x))}


  #close(pb)

  #naming chunks based on readnames & indices
  for (chunk in seq_along(tail_chunk_list)) {
    names(tail_chunk_list[[chunk]]) <- paste(names(tail_chunk_list)[chunk], seq_along(tail_chunk_list[[chunk]]), sep = "_")
  }

  chunk_names <- gsub(".*\\.", "", names(rapply(tail_chunk_list, class)))
  names(gaf_list) <- chunk_names

  return(gaf_list)
}

#' Filters read chunks containing nonadenosine nucleotides of interest for
#' neural net training set preparation.
#'
#' The function is designed to be used on a generated set of synthetic
#' spike-ins (laboratory produced) containing a particular type of residue
#' (G, C or U in the context of 3'-homopolymer A, respectively).
#'
#' IMPORTANT NOTE!
#' The function is not suitable for preparing a set containing
#' pure (containing only A) fragments of polyA tails.
#'
#' This function currently allows to produce the filtered list of chunks
#' containing the signal deviations corresponding to the given non-A nucleotide
#' of interest. In its current form, the function allows filtering signal
#' fragments having C, G or U nucleotides (one category at a time).
#'
#' The function uses as input a list of fragments (chunks) centered on a
#' deviation in the signal produced by create_tail_chunk_list() function.
#'
#' It takes advantage of the fact that when filtering a signal based on moving
#' average & standard deviation filters, a characteristic pattern of "pseudomoves"
#' is produced. Which corresponds to changes in the signal visible even
#' to the naked eye. Empirically, it has been established that the presence
#' of nucleotide G results in a peak, while the presence of nucleotides C
#' and U results in a valley. These deviations mostly have specific parameters
#' (depth and width) based on which filtering criteria can be established.
#'
#' Before performing filtering, the user needs to know what type of nucleotide
#' (G, C or U) they are dealing with in a given data set. Otherwise, it may be
#' difficult to interpret the results, and the network training itself may
#' not lead to satisfactory results.
#'
#' The filtering procedure is controlled by the value parameter. It determines
#' whether chunks containing a peak (value then takes the value 1) or a valley
#' (value takes the value -1) are kept in the output.
#'
#' The function does not distinguish between fragments containing C and U
#' (both produce valley; the classification is handled by CNN). Therefore,
#' it is important that the training datasets differ in the transcript bodies
#' to which they can be mapped and/or are not run together in a single
#' sequencing run.
#'
#' Before proceeding further, it is advisable to perform a visual inspection
#' of at least some of the filtered signals. It may be necessary to manually
#' adjust the hardcoded parameters of other functions (contact dev/maintainer
#' for further details).
#'
#' @param tail_chunk_list character string. The list object produced
#' by create_chunk_list function.
#'
#' @param value numeric [1]. A parameter that controls whether valleys (C, U)
#' or peaks (G) are filtered. For C, U nucleotides it takes the value of -1,
#' for G nucleotide it takes the value of 1.
#'
#' @param num_cores numeric [1]. Number of physical cores to use in processing
#' the data. Do not exceed 1 less than the number of cores at your disposal.
#'
#' @return a list of signal chunks filtered based on the user-defined value
#' parameter.
#'
#' @export
#'
#' @examples
#'\dontrun{
#'
#' # filtering G residue:
#' filter_nonA_chunks_trainingset(tail_chunk_list = list_object_with_tail_chunks,
#'                                value = 1, num_cores = 2)
#'
#' # filtering C/U (user must know which type of residue is dealing with):
#' filter_nonA_chunks_trainingset(tail_chunk_list = list_object_with_tail_chunks,
#'                                value = -1, num_cores = 2)
#' }
filter_nonA_chunks_trainingset <- function(tail_chunk_list,
                                           value,
                                           num_cores){
  #variable binding
  i <- NULL

  #assertions
  if (missing(num_cores)) {
    stop("Number of declared cores is missing. Please provide a valid num_cores argument.", call. =FALSE)
  }

  if (missing(tail_chunk_list)) {
    stop("List of tail chunks is missing. Please provide a valid tail_chunk_list argument.", call. =FALSE)
  }

  assertthat::assert_that(is.list(tail_chunk_list),
                          msg = paste0("Given tail_chunk_list is not a list (class). Please provide valid file format."))
  assertthat::assert_that(is.numeric(num_cores),
                          msg=paste0("Declared core number must be numeric. Please provide a valid argument."))
  assertthat::assert_that(is.numeric(value),
                          msg=paste0("Provided value must be numeric. Please provide a valid argument."))
  assertthat::assert_that(value %in% c(-1,1),
                          msg=paste0("Provided value must be either 1 or -1. Please provide a valid argument."))

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
  cat(paste0('[', as.character(Sys.time()), '] ','Filtering training nonA chunks...', '\n', sep=''))

  # progress bar
  pb <- utils::txtProgressBar(min = 0,
                              max = length(tail_chunk_list),
                              style = 3,
                              width = 50,
                              char = "=",
                              file = stderr())
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  #filter the potential nonA containing chunks for training set preparation
  filtered_input <- foreach::foreach(i = seq_along(tail_chunk_list),
                                     .inorder = TRUE,
                                     .errorhandling = 'pass',
                                     .options.snow = opts,
                                     .options.multicore = mc_options) %dopar% {
    filtered_output <-  Filter(function(x) any(with(rle(x$pseudomoves), lengths[values==value]>=4)) & x$chunk_start_pos>=0, tail_chunk_list[[i]])
    lapply(filtered_output, function(x) x)
  }

  #close(pb)

  #remove empty sublists
  filtered_input <- Filter(function(x) length(x) > 0, filtered_input)

  #restore names in the list
  chunknames <- purrr::map_depth(filtered_input, 1, names) %>% unlist(use.names = F)
  chunknames <- unique(gsub("\\_[0-9]*$","",chunknames))
  names(filtered_input) <- chunknames

  # Done comm
  cat(paste0('[', as.character(Sys.time()), '] ','Done!', '\n', sep=''))

  return(filtered_input)
}


#' Filters out signals of given type of interest.
#'
#' Filters out signals corresponding to a given category of reads
#' (containing A-nucleotides alone or particular types of non-adenosine
#' nucleotides: G, C or U) in order to prepare the set for neural network
#' training.
#'
#' @param nucleotide character. One of the following ["A"/"C"/"G","U"].
#' This parameter defines the type of filtering approach applied to the input
#' dataset in order to obtain the signals potentially carrying the
#' desired nucleotide context.
#'
#' @param nanopolish character string. Full path of the .tsv file produced
#' by nanopolish polya function.
#'
#' @param sequencing_summary character string. Full path of the .txt file
#' with sequencing summary.
#'
#' @param workspace character string. Full path of the directory to search the
#' basecalled fast5 files in. The Fast5 files have to be multi-fast5 file.
#'
#' @param num_cores numeric [1]. Number of physical cores to use in processing
#' the data. Do not exceed 1 less than the number of cores at your disposal.
#'
#' @param basecall_group character string ["Basecall_1D_000"]. Name of the
#' level in the fast5 file hierarchy from which the data should be extracted.
#'
#' @param pass_only logical [TRUE/FALSE]. If TRUE, only reads tagged by
#' nanopolish as "PASS" would be taken into consideration. Otherwise, reads
#' tagged as "PASS" & "SUFFCLIP" will be taken into account in analysis.
#' As a default, "TRUE" value is set.
#'
#' @return A list of filtered GAF matrices organized by the read ID_index
#' is returned. Always assign this returned list to a variable, otherwise
#' the long list will be printed to the console, which may crash your R session.
#'
#' @export
#'
#' @examples
#'\dontrun{
#'
#' prepare_trainingset(nucleotide="A"
#'                     nanopolish = '/path/to/file',
#'                     sequencing_summary = '/path/to/file',
#'                     workspace = '/path/to/guppy/workspace',
#'                     num_cores = 10,
#'                     basecall_group = 'Basecall_1D_000',
#'                     pass_only=TRUE)
#'}
#'
prepare_trainingset <- function(nucleotide,
                                nanopolish,
                                sequencing_summary,
                                workspace,
                                num_cores=1,
                                basecall_group="Basecall_1D_000",
                                pass_only=TRUE){

  # nucleotide selection
  if (nucleotide=="A"){
    # process the A-containing data; split tails with overlaps
    tail_feature_list <- ninetails::create_tail_feature_list_A(nanopolish,sequencing_summary,workspace,num_cores,basecall_group,pass_only)
    tail_chunk_list <- ninetails::create_tail_chunk_list_A(tail_feature_list, num_cores)
    gafs_list <- ninetails::create_gaf_list_A(tail_chunk_list, num_cores)

  } else if (nucleotide=="C") {
    # process the C-containing data to filter only internal, C-containing chunks
    tail_feature_list <- ninetails::create_tail_feature_list_trainingset(nanopolish,sequencing_summary,workspace,num_cores,basecall_group,pass_only)
    tail_chunk_list <- ninetails::create_tail_chunk_list_trainingset(tail_feature_list, num_cores)
    filtered_chunk_list <- ninetails::filter_nonA_chunks_trainingset(tail_chunk_list, value=-1, num_cores)
    gafs_list <- ninetails::create_gaf_list(filtered_chunk_list, num_cores)

  } else if (nucleotide=="G") {
    # process the G-containing data to filter only internal, G-containing chunks
    tail_feature_list <- ninetails::create_tail_feature_list_trainingset(nanopolish,sequencing_summary,workspace,num_cores,basecall_group,pass_only)
    tail_chunk_list <- ninetails::create_tail_chunk_list_trainingset(tail_feature_list, num_cores)
    filtered_chunk_list <- ninetails::filter_nonA_chunks_trainingset(tail_chunk_list, value=1, num_cores)
    gafs_list <- ninetails::create_gaf_list(filtered_chunk_list, num_cores)

  } else if (nucleotide=="U") {
    # process the U-containing data to filter only internal, U-containing chunks
    tail_feature_list <- ninetails::create_tail_feature_list_trainingset(nanopolish,sequencing_summary,workspace,num_cores,basecall_group,pass_only)
    tail_chunk_list <- ninetails::create_tail_chunk_list_trainingset(tail_feature_list, num_cores)
    filtered_chunk_list <- ninetails::filter_nonA_chunks_trainingset(tail_chunk_list, value=-1, num_cores)
    gafs_list <- ninetails::create_gaf_list(filtered_chunk_list, num_cores)

  } else {
    stop("Wrong nucleotide selected. Please provide valid nucleotide argument.", call. =FALSE)
  }

  return(gafs_list)
}
