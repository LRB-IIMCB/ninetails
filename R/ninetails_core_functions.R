#' Extracts features of multiple RNA reads from nanopolish output
#' & sequencing summary file.
#'
#' This function extracts features of poly(A) tails of selected RNA reads
#' from the output table provided by nanopolish polya function and sequencing
#' summary provided by the sequencer. Filenames are taken from the sequencing
#' summary file. Only reads with tail lengths estimated as >= 10nt by nanopolish
#' polya function are taken into account.
#'
#' @param nanopolish character string. Either full path of the .tsv file produced
#' by nanopolish polya function or an environment object containing nanopolish
#' data.
#'
#' @param sequencing_summary character string. Either full path of the .txt file
#' with sequencing summary or an environment object containing sequencing
#' summary data.
#'
#' @param pass_only logical [TRUE/FALSE]. If TRUE, only reads tagged by
#' nanopolish as "PASS" would be taken into consideration. Otherwise, reads
#' tagged as "PASS" & "SUFFCLIP" will be taken into account in analysis.
#' As a default, "TRUE" value is set.
#'
#' @return A tibble/df containing read information organized by the read ID
#' is returned. Always assign this returned tibble to a variable. Otherwise
#' the long tibble will be printed to the console, which may crash your
#' R session.
#'
#' @export
#'
#' @examples
#'\dontrun{
#'
#' ninetails::extract_polya_data(nanopolish = '/path/to/nanopolish/polya/output.tsv',
#'                    sequencing_summary = '/path/to/sequencing_summary.txt',
#'                    pass_only = TRUE)
#'
#'}
extract_polya_data <- function(nanopolish,
                               sequencing_summary,
                               pass_only = TRUE){

  # variable binding (suppressing R CMD check from throwing an error)
  readname <- polya_start <- transcript_start <- polya_length <- qc_tag  <- filename <- read_id <- NULL

  if (missing(nanopolish)) {
    stop("Nanopolish polya output is missing. Please provide a valid nanopolish argument.", .call = FALSE)
  }

  if (missing(sequencing_summary)) {
    stop("Sequencing summary file is missing. Please provide a valid sequencing_summary argument.", .call = FALSE)
  }

  assertthat::assert_that(is.logical(pass_only),
                          msg="Please provide TRUE/FALSE values for pass_only parameter")

  # Accept either path to file or in-memory file - PK's GH issue
  if (checkmate::test_string(nanopolish)) {
    # if string provided as an argument, read from file
    checkmate::assert_file_exists(nanopolish)
    nanopolish_polya_table <- vroom::vroom(nanopolish,
                                           col_select=c(readname, polya_start, transcript_start, polya_length, qc_tag),
                                           show_col_types = FALSE)
  } else {
    # make sure that nanopolish is an object with rows
    if (!is.data.frame(nanopolish) || nrow(nanopolish) == 0) {
      stop("Empty data frame provided as an input (nanopolish). Please provide valid input")
    }

    nanopolish_polya_table <- nanopolish[,c("readname","polya_start","transcript_start","polya_length","qc_tag")]
  }

  # Accept either path to file or in-memory file - PK's GH issue
  if (checkmate::test_string(sequencing_summary)) {
    # if string provided as an argument, read from file
    checkmate::assert_file_exists(sequencing_summary)
    sequencing_summary_table <- vroom::vroom(sequencing_summary,
                                             col_select = c(filename, read_id),
                                             show_col_types = FALSE)
  } else {
    if (!is.data.frame(sequencing_summary) || nrow(sequencing_summary) == 0) {
      stop("Empty data frame provided as an input (sequencing_summary). Please provide valid input")
    }

    sequencing_summary_table <- sequencing_summary
  }

  #rename read id column
  names(sequencing_summary_table)[names(sequencing_summary_table) == "read_id"] <- "readname"

  # Add filtering criterion: select only pass or pass $ suffclip
  if(pass_only == TRUE){
    polya_summary <- dplyr::left_join(nanopolish_polya_table[which(nanopolish_polya_table$qc_tag=="PASS"),],
                                      sequencing_summary_table,
                                      by="readname")
  } else {
    polya_summary <- dplyr::left_join(nanopolish_polya_table[which(nanopolish_polya_table$qc_tag %in% c("PASS", "SUFFCLIP")),],
                                      sequencing_summary_table,
                                      by="readname")
  }

  # Add filtering criterion: tail length >= 10 nt
  polya_summary <- dplyr::filter(polya_summary, polya_length >=10)

  # Prevent bugs from custom models used for mapping (secondary alignments)
  polya_summary <- unique(polya_summary)

  names(polya_summary$filename) <- polya_summary$readname #named vec of filenames (names = readnames)
  attr(polya_summary, 'spec') <- NULL #drop attributes left by vroom

  return(polya_summary)
}


#' Extracts tail features of single RNA read from respective basecalled
#' multi-fast5 file.
#'
#' This function extracts metadata of RNA read from multi-fast5 file basecalled
#' by guppy. Filenames are taken from the sequencing
#' summary file. Tail signal (as delimited by nanopolish polya function) is
#' extracted (before the extraction step, the entire signal is winsorized
#' so the cliffs are removed), and downsampled (to facilitate further analysis).
#'
#' @param readname character string. Name of the given read within the
#' analyzed dataset.
#'
#' @param polya_summary character string. The table containing data extracted
#' from nanopolish & sequencing summary (using \code{\link{extract_polya_data}} function.
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
#' ninetails::extract_tail_data(readname = 'abc123de-fg45-6789-0987-6543hijk2109',
#'                              polya_summary = polya_summary_table,
#'                              workspace = '/path/to/folder/containing/multifast5s',
#'                              basecall_group = 'Basecall_1D_000')
#'
#'}
extract_tail_data <- function(readname,
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

  if (!is.data.frame(polya_summary) || nrow(polya_summary) == 0){
    stop("Empty data frame provided as an input (polya_summary). Please provide a valid input table.",call. =FALSE )
  }


  assertthat::assert_that(is.character(workspace),
                          msg = paste0("Path to basecalled fast5 files is not a character string. Please provide a valid path to basecalled fast5 files."))
  assertthat::assert_that(checkmate::test_string(workspace, null.ok=F, min.chars=1),
                          msg="Empty string provided as an input. Please provide a valid path to basecalled fast5 files.")
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
  basecall_1d_template <- rhdf5::h5readAttributes(fast5_file_path,paste0(fast5_readname,"/Analyses/", basecall_group,"/Summary/basecall_1d_template")) # parent dir for attributes (within fast5); fixed hardcoding bug
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
  pseudomoves <- ninetails::filter_signal_by_threshold(signal)

  extracted_data_single_list = list() # creating empty list for the extracted fast5 data

  extracted_data_single_list[["fast5_filename"]] <- polya_summary$filename[read_idx]
  extracted_data_single_list[["tail_signal"]] <- signal
  extracted_data_single_list[["tail_moves"]] <- moves_tail_range
  #added for robustness
  extracted_data_single_list[["tail_pseudomoves"]] <- pseudomoves
  # excluded to reduce the output
  #extracted_data_single_list[["polya_start_pos"]] <- polya_start_position
  #extracted_data_single_list[["transcript_start_pos"]] <- transcript_start_position
  #extracted_data_single_list[["tail_length"]] <- polya_summary$polya_length[read_idx]

  return(extracted_data_single_list)
}


#' Extracts features of polyA tails of ONT RNA reads required for finding
#' non-A nucleotides within the given tails.
#'
#' This function extracts tail features of RNA reads from multi-fast5 files
#' basecalled by guppy and polyA tail characteristics (coordinates) created
#' by nanopolish polya function. Filenames are taken from the sequencing
#' summary file.
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
#'tfl <- ninetails::create_tail_feature_list(
#'  nanopolish = system.file('extdata',
#'                           'test_data',
#'                           'nanopolish_output.tsv',
#'                           package = 'ninetails'),
#'  sequencing_summary = system.file('extdata',
#'                                   'test_data',
#'                                   'sequencing_summary.txt',
#'                                   package = 'ninetails'),
#'  workspace = system.file('extdata',
#'                          'test_data',
#'                          'basecalled_fast5',
#'                          package = 'ninetails'),
#'  num_cores = 2,
#'  basecall_group = 'Basecall_1D_000',
#'  pass_only=TRUE)
#'
#'}
#'
create_tail_feature_list <- function(nanopolish,
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
                                         .options.multicore = mc_options) %dopar% {lapply(polya_summary$readname[i], function(x) ninetails::extract_tail_data(x,polya_summary,workspace,basecall_group))}

  #close(pb)

  #label each signal according to corresponding read name to avoid confusion
  # this deals with issue#5 (nanopolish records do not exactly match seqsummary file)
  squiggle_names <- as.vector(sapply(tail_features_list, function(x) attributes(x[[1]])$names))
  tail_features_list <- stats::setNames(tail_features_list, squiggle_names)

  #former solution:
  #squiggle_names <- polya_summary$readname
  #names(tail_features_list) <- squiggle_names

  # remove reads with only zero moved tails
  tail_features_list <- Filter(function(x) sum(x$tail_moves) !=0, tail_features_list)
  zeromoved_readnames <- squiggle_names[!(squiggle_names %in% names(tail_features_list))]

  # prevent from running on reads which do not fulfill the pseudomove condition
  tail_features_list <- Filter(function(x)any(with(rle(x$tail_pseudomoves), lengths[values!=0]>=5)), tail_features_list)

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
#' The function enables detection in the signal corresponding to the poly(A)
#' tail of areas in which the signal value significantly deviates from the
#' values typical of an adenosine homopolymer.
#'
#' The function results in a vector of "pseudomoves" with values in the range
#' -1:1, where -1 corresponds to signals that are significantly less than the
#' mean and standard deviation of the typical A-homopolymer, 1 corresponds to
#' signals that are significantly higher than the mean and standard deviation
#' of the typical A-homopolymer, and 0 corresponds to the typical values
#' for the A-homopolymer.
#'
#' The "pseudomoves" vector allows more accurate calibration of the nucleotide
#' positions of potential non-adenosine residues than the moves produced
#' by the guppy basecaller.
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
#' \code{\link{create_tail_feature_list}} function.
#'
#' @return numeric vector of "pseudomoves" corresponding to the analyzed tail
#' region; containing values ranging from -1 to 1.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' ninetails::filter_signal_by_threshold(signal=tail_feature_list[[1]][["readname"]][[4]])
#'
#'}
filter_signal_by_threshold <- function(signal) {

  # variable binding
  x <- sd <- baseline <- std_cutoff <- NULL

  #assertions
  if (missing(signal)) {
    stop("Signal is missing. Please provide a valid signal argument [numeric vec].", call. =FALSE)
  }

  assertthat::assert_that(is.numeric(signal),
                          msg=paste0("Signal must be numeric. Please provide a valid argument."))

  # reproducibility
  set.seed(123)

  # calibrate the algo on the sampled vals
  start_vals <- signal[1:10]
  most_freq_vals <- as.numeric(names(sort(table(signal),decreasing=TRUE)[1:20]))
  adjusted_signal <- c(sample(c(most_freq_vals, start_vals), 100, replace=TRUE), signal)

  #interpolatedsignal <- round(stats::approx(signal, method = "linear", n=ceiling(0.05 * length(signal)))[[2]], digits=0)
  #adjusted_signal <- c(sample(interpolatedsignal, 100, replace=TRUE), signal)

  # Empyrical parameters:
  adaptive_sampling_window <- 100 # datapoints window for adjusting algo
  SD_threshold <- 3.5 # how many SD tresholds from avg signal pseudomove should be reported

  pseudomoves <- rep(0,length(adjusted_signal))
  #filtered_signal <- adjusted_signal[1:adaptive_sampling_window] #first sampling
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

#' Extracts fragments of poly(A) tail signal containing potential modifications
#' along with its delimitation (positional indices; coordinates) within the tail.
#'
#' This function finds areas in the poly(A) tail signal containing potential
#' non-A residues. Extracts strings spanning 100 data points (extrapolated)
#' where the potential modification is always at the center of a given
#' extracted fragment.
#'
#' It extracts the region of interest based on 2 assumptions (conditions):
#' the presence of significant raw signal distortion (above the given threshold)
#' recorded by the thresholding algo as "pseudomove" and the transition of state
#' (move==1) recorded by Guppy. If move==0 are present exclusively within the given
#' signal chunk, then this chunk is dropped from the analysis (as the distortion
#' is most likely caused by sequencing artifact, not the non-A residue itself
#' - this is introduced in order to minimize alse positives)
#'
#' If the data indicating the presence of modifications are in the regions
#' of the ends of the signal (3' or 5'), then the missing upstream or downstream
#' data are imputed based on the most frequent values in the entire signal.
#'
#' The function returns a list object (nested) which consists of the individual
#' fragments of the input signal (chunk_sequence) and the coordinates
#' of the given fragment in the input signal: start position (chunk_start_pos)
#' and end position (chunk_end_pos).
#'
#' @param readname character string. Name of the given read within the
#' analyzed dataset.
#'
#' @param tail_feature_list list object produced by \code{\link{create_tail_feature_list}}
#' function.
#'
#' @return a list object (nested) containing 3 categories: resulting fragments
#' (chunk_sequence), start coordinate of given fragment (chunk_start_pos)
#' its' end coordinate (chunk_end_pos) arranged by the signal ID and positional
#' indices (WARNING! from 3' end!).
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' ninetails::split_tail_centered(readname= "1234-anexample-r3adn4m3",
#'                     tail_feature_list = tail_feature_list)
#'
#'}
split_tail_centered <- function(readname,
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
  moves <- tail_feature_list[[1]][[readname]][[3]]
  pseudomoves <- tail_feature_list[[1]][[readname]][[4]]


  # mod-centered chunk extraction
  # recompute rle
  mod_rle <- rle(pseudomoves)
  # pseudomoves filtered by condition (potentially decorated - empyrical!)
  condition <- mod_rle$lengths >= 5 & mod_rle$values
  # beginning positions of filtered pseudomoves which satisfy conditions
  first_filtered_positions <- cumsum(c(1, utils::head(mod_rle$lengths,-1)))[condition]
  # length of pseudomoves satisfying condition
  filtered_length <- mod_rle$lengths[condition]
  # extracted coordinates (indices)
  start_positions <-  first_filtered_positions + floor(filtered_length/2) - 50
  end_positions <- start_positions + 99

  # extract signal chunks centered on potential modification
  list_1 <- lapply(1:length(start_positions), function(i) if(start_positions[i] >0) signal[start_positions[i]:end_positions[i]]
                   else c(rep(NA, abs(start_positions[i]-1)), signal[1:end_positions[i]]))

  # replace NAs with 5 most freq values (impute the data for gafs)
  most_freq_vals <- as.numeric(names(sort(table(signal),decreasing=TRUE)[1:5]))
  list_1 <- lapply(list_1, function(n) replace(n, is.na(n), sample(most_freq_vals,sum(is.na(n)),TRUE)))

  #add indices
  chunks_indices <- c(1:length(start_positions))

  # naming chunks based on names & indices
  chunk_names <- paste0(rep(readname, length(list_1)), '_', unlist(chunks_indices))
  names(list_1) <- chunk_names

  # retrieve coordinates as list_2 and list_3:
  list_2 <- as.list(start_positions)
  #names(list_2) <- chunk_names
  list_3 <- as.list(end_positions)
  #names(list_3) <- chunk_names

  # retrieve move vector fragments based on the same positional info:
  list_4 <- lapply(1:length(start_positions), function(i) if(start_positions[i] >0) moves[start_positions[i]:end_positions[i]]
                   else c(rep(NA, abs(start_positions[i]-1)), moves[1:end_positions[i]]))

  out <- mget(ls(pattern = '^list.*\\d$')) %>%
    split(sub("_\\d+$", '', names(.))) %>%
    purrr::map(~purrr::transpose(purrr::set_names(.,c('chunk_sequence', 'chunk_start_pos', 'chunk_end_pos', 'chunk_moves')))) %>%
    purrr::flatten(.)


  # remove chunks with discordance between moves & pseudomoves reported
  # if pseudomoves registered while moves not - remove chunk
  # this usually causes the loss of terminal positions, but this is to avoid
  # inheritance of nanopolish error as well as the warn labeling
  out <- Filter(function(x) sum(x$chunk_moves) !=0, out)

  return(out)
}


#' Creates list of polyA tail chunks centered on significant signal deviations.
#'
#' Extracts fragments of polyA tails of ONT RNA reads potentially containing
#' non-A nucleotides along their coordinates & appends the data to the nested
#' list organized by read IDs.
#'
#' @param tail_feature_list list object produced by \code{\link{create_tail_feature_list}}
#' function.
#' @param num_cores numeric [1]. Number of physical cores to use in processing
#' the data. Do not exceed 1 less than the number of cores at your disposal.
#'
#' @return a list object (nested) containing the segmented tail data (chunks,
#' coordinates) organized by the read IDs.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' tcl <- ninetails::create_tail_chunk_list(tail_feature_list = tfl,
#'                                          num_cores = 3)
#'
#'}
create_tail_chunk_list <- function(tail_feature_list,
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

  #output list
  tail_chunk_list <- list()

  #parallel extraction
  tail_chunk_list <- foreach::foreach(i = seq_along(tail_feature_list[[1]]),
                                      .combine = c,
                                      .inorder = TRUE,
                                      .errorhandling = 'pass',
                                      .options.snow = opts,
                                      .options.multicore = mc_options) %dopar% {
                                        lapply(names(tail_feature_list[[1]][i]), function(x) ninetails::split_tail_centered(x,tail_feature_list))
                                      }

  #close(pb)

  #rename first level of nested list accordingly
  names(tail_chunk_list) <- names(tail_feature_list[[1]])


  #drop moves variable - recursive fn for prunning moves vector from the output
  # this is to save the memory space and minimize the size of vars
  .prune_moves <- function(i)
    lapply(i, function(x)
      if (is.list(x)) {
        if(!is.null(names(x))) .prune_moves(x[names(x)!="chunk_moves"]) else .prune_moves(x)
      } else x
    )

  tail_chunk_list <- .prune_moves(tail_chunk_list)

  #remove empty elements from the list (clean the output)
  tail_chunk_list <-rrapply::rrapply(tail_chunk_list, condition = Negate(is.null), how = "prune")


  # Done comm
  cat(paste0('[', as.character(Sys.time()), '] ','Done!', '\n', sep=''))

  return(tail_chunk_list)
}


#' Converts ONT signal to Gramian Angular Field.
#'
#' This function represents time series data (ont squiggle) in a polar coordinate
#' system instead of the typical Cartesian coordinates. This is a Pythonic
#' pyts.image equivalent written in R language.
#'
#' Two methods of such transformation are available: Gramian angular summation
#' field (GASF) and Gramian angular difference field (GADF).
#'
#' @param tail_chunk numeric. A numeric vector representing signal chunk
#' within the analyzed dataset.
#'
#' @param method character string specifying the type of Gramian Angular Field:
#' "s" can be used to produce summation field (GASF) and "d" to produce
#' difference field (GADF). Defaults to summation ["s"].
#'
#' @return an array (100,100,1) with values (coordinates) representing ONT signal.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' ninetails::create_gaf(tail_chunk = tail_chunk,
#'                       method="s")
#'
#'}

create_gaf <- function(tail_chunk, method="s"){

  #assertions
  if (missing(tail_chunk)) {
    stop("Tail_chunk is missing. Please provide a valid tail_chunk argument.", call. =FALSE)
  }
  if (missing(method)) {
    stop("Transformation method is missing. Please provide a valid method argument.", call. =FALSE)
  }
  assertthat::assert_that(is.numeric(tail_chunk),
                          msg=paste0("Provided tail_chunk must be numeric. Please provide a valid argument."))
  assertthat::assert_that(is.character(method),
                          msg=paste0("Provided method must be character string. Please provide a valid argument."))
  assertthat::assert_that(length(tail_chunk) == 100,
                          msg=paste0("Provided chunks of wrong length. The chunk length should be equal to 100. Please provide a valid tail_chunk."))


  # rescale values so that all of them fall in the interval [-1, 1]:
  tail_chunk <- (tail_chunk-max(tail_chunk)+(tail_chunk-min(tail_chunk)))/(max(tail_chunk)-min(tail_chunk))

  if (method=="s") {
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
    tail_chunk <- array(t(tail_chunk), c(100,100,1))

  } else if (method=="d") {
    # computing arc sine instead of arc cosine
    # calculating sine instead of cosine
    tail_chunk <- asin(tail_chunk)
    tail_chunk <- cbind(replicate(length(tail_chunk), tail_chunk))
    tail_chunk <- tail_chunk + t(tail_chunk)
    tail_chunk <- sin(tail_chunk)
    tail_chunk <- array(t(tail_chunk), c(100,100,1))
  } else {
    stop("Wrong GAF method definition. Please provide 's' for summation or 'd' for difference GAF.")
  }

  # rescale values so that all of them fall in the interval [0, 1]:
  tail_chunk <- round((tail_chunk-min(tail_chunk))/(max(tail_chunk)-min(tail_chunk)), 4)

  return(tail_chunk)
}

#' Creates a two-dimensional array containing GASF and GADF resulting from
#' the transformation of a given ONT tail chunk.
#'
#' This function allows for the classification of signal fragments based
#' on angular coordinates generated by two methods (summation & difference)
#' simultaneously.
#'
#' Using this approach increases the sensitivity of the classification.
#' It overcomes the limitations of each method.
#'
#' @param tail_chunk numeric. A numeric vector representing signal chunk
#' within the analyzed dataset.
#'
#' @return an array (100,100,2) with values (coordinates) representing GASF
#' (first dimension) and GADF (second dimension) produced by the
#' \code{\link{create_gaf}} function applied to given fragment (tail chunk)
#' of analyzed ONT signal.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' ninetails::combine_gafs(tail_chunk = tail_chunk)
#'
#'}
combine_gafs <- function(tail_chunk){

  #assertions
  if (missing(tail_chunk)) {
    stop("Tail_chunk is missing. Please provide a valid tail_chunk argument.", call. =FALSE)
  }

  assertthat::assert_that(is.numeric(tail_chunk),
                          msg=paste0("Provided tail_chunk must be numeric. Please provide a valid argument."))

  #create gasf & gaf
  GASF <- ninetails::create_gaf(tail_chunk=tail_chunk, method="s")
  GADF <- ninetails::create_gaf(tail_chunk=tail_chunk, method="d")

  #create array of gafs
  combined_gafs <- array(c(GASF, GADF), dim = c(100, 100, 2))

  return(combined_gafs)
}



#' Creates list of gramian angular matrices produced based on
#' list of splitted tails (tail chunks).
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
#' @importFrom foreach %dopar%
#' @importFrom utils head
#'
#' @export
#'
#' @examples
#'\dontrun{
#'
#' gl <- ninetails::create_gaf_list(tail_chunk_list = tcl,
#'                                  num_cores = 2)
#'
#'}
create_gaf_list <- function(tail_chunk_list,
                            num_cores){

  #variable biding
  i <- NULL

  # Assertions
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
  cat(paste0('[', as.character(Sys.time()), '] ','Computing gramian angular fields...', '\n', sep=''))

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
                               .options.multicore = mc_options) %dopar% {
                                 lapply(tail_chunk_list[[i]], function(x_ij) ninetails::combine_gafs(x_ij[['chunk_sequence']]))
                               }

  #close(pb)

  # Done comm
  cat(paste0('[', as.character(Sys.time()), '] ','Done!', '\n', sep=''))

  return(gaf_list)
}

#' Performs classification of matrices produced from read chunks with CNN.
#'
#' This function allows to assign gafs corresponding to the given signals
#' into one of 4 categories (A, C, G, U, respectively). This function in its
#' current implementation allows to load the pretrained CNN model. It uses
#' tensorflow backend to produce the classification output.
#'
#' @param gaf_list [list] A list of gaf matrices organized by the read ID_index.
#'
#' @return a list of gafs predictions based on used model.
#'
#' @export
#'
#' @examples
#'\dontrun{
#'
#' pl <- ninetails::predict_gaf_classes(gl)
#'
#'}
predict_gaf_classes <- function(gaf_list){

  #assertions
  if (missing(gaf_list)) {
    stop("List of transformed signal chunks is missing. Please provide a valid gaf_list argument.", call. =FALSE)
  }

  chunknames <- names(gaf_list)
  names(gaf_list)  <- NULL

  gaf_list <- simplify2array(gaf_list)
  gaf_list <- aperm(gaf_list, c(4,1,2,3))

  # Output info
  cat(paste0('[', as.character(Sys.time()), '] ','Classifying gramian angular fields...', '\n', sep=''))

  keras_model <- ninetails::load_keras_model()

  #predict chunk class
  predicted_gaf_classes <- keras_model %>% stats::predict(gaf_list) %>% keras::k_argmax()
  predicted_gaf_classes <- as.numeric(predicted_gaf_classes)

  predicted_list = list() # creating empty list for the extracted  data

  predicted_list[["chunkname"]] <- chunknames
  predicted_list[["prediction"]] <- predicted_gaf_classes

  # Output info
  cat(paste0('[', as.character(Sys.time()), '] ','Done!', '\n', sep=''))

  return(predicted_list)
}


#' Creates the list object containing tabular outputs of ninetails pipeline.
#'
#' The read_classes dataframe contains results of read classification.
#' The sequencing reads are assigned to classes based on whether the initial
#' conditions are met or not (e.g. sufficient read quality, sufficient length
#' (>=10 nt), move transition presence, local signal anomaly detected etc.).
#' According to this, reads are assigned into 3 main categories: decorated,
#' blank, unclassified (class column).
#'
#' The more detailed info regarding read classification is stored within
#' 'comments' column. To make the output more compact, it contains codes
#' as follows:\itemize{
#' \item IRL - insufficient read length
#' \item QCF - nanopolish qc failed
#' \item MAU - move transition absent, nonA residue undetected
#' \item MPU - move transition present, nonA residue undetected
#' \item NIN - not included in the analysis (pass only = T)
#' \item YAY - move transition present, nonA residue detected
#' }
#'
#' The nonadenosine_residues contains detailed positional info regarding all
#' potential nonadenosine residues detected.
#'
#' @param tail_feature_list list object produced by \code{\link{create_tail_feature_list}}
#' function.
#'
#' @param tail_chunk_list list object produced by \code{\link{create_tail_chunk_list}}
#' function.
#'
#' @param nanopolish character string. Full path of the .tsv file produced
#' by nanopolish polya function or name of in-memory file (environment object).
#'
#' @param predicted_list a list object produced by \code{\link{predict_classes}} function.
#'
#' @param num_cores numeric [1]. Number of physical cores to use in processing
#' the data. Do not exceed 1 less than the number of cores at your disposal.
#'
#' @param pass_only logical [TRUE/FALSE]. If TRUE, only reads tagged by
#' nanopolish as "PASS" would be taken into consideration. Otherwise, reads
#' tagged as "PASS" & "SUFFCLIP" will be taken into account in analysis.
#' As a default, "TRUE" value is set.
#'
#' @param qc logical [TRUE/FALSE]. If TRUE, the quality control of the output
#' predictions would be performed. This means that the reads/non-A residue
#' positions in terminal nucleotides, which are most likely artifacts, are
#' labeled accordingly as "-WARN" (residues recognized as non-A due to
#' nanopolish segmentation error which is inherited from nanopolish,
#' as ninetails uses nanopolish segmentation). It is then up to user, whether
#' they would like to include or discard such reads from their pipeline. However,
#' it is advised to treat them with caution. By default, the qc option is enabled
#' (this parameter is set to TRUE).
#'
#' @return This function returns a list object containing two dataframes:
#' "read_classes" and "nonadenosine_residues" with the final output.
#' First dataframe contains initial indications, whether the given read was
#' classified or omitted (with reason) and if classified, whether read was
#' recognized as decorated (containing non-adenosine residue) or not.
#' The second dataframe contains detailed info on type and estimated positions
#' of non-adenosine residues detected.
#'
#' @export
#'
#' @examples
#'\dontrun{
#'
#' create_outputs(tail_feature_list = tail_feature_list,
#'                tail_chunk_list = tail_chunk_list,
#'                nanopolish = '/path/to/nanopolish_output.tsv',
#'                predicted_list = predicted_list,
#'                num_cores = 2,
#'                pass_only=TRUE,
#'                qc=TRUE)
#'}
#'
#'
create_outputs <- function(tail_feature_list,
                           tail_chunk_list,
                           nanopolish,
                           predicted_list,
                           num_cores,
                           pass_only=TRUE,
                           qc=TRUE){

  #variable binding
  readname <- polya_length <- qc_tag<- i <- chunkname <- contig <- est_nonA_pos <- NULL

  #assertions

  if (missing(tail_feature_list)) {
    stop("List of tail features is missing. Please provide a valid tail_feature_list argument.", call. =FALSE)
  }

  if (missing(tail_chunk_list)) {
    stop("List of tail chunks is missing. Please provide a valid tail_chunk_list argument.", call. =FALSE)
  }

  if (missing(nanopolish)) {
    stop("Nanopolish polya output is missing. Please provide a valid nanopolish argument.", .call = FALSE)
  }

  if (missing(predicted_list)) {
    stop("List of predictions is missing. Please provide a valid predicted_list argument.", call. =FALSE)
  }

  if (missing(num_cores)) {
    stop("Number of declared cores is missing. Please provide a valid num_cores argument.", call. =FALSE)
  }

  assertthat::assert_that(is.list(tail_feature_list),
                          msg = paste0("Given tail_feature_list is not a list (class). Please provide valid object."))
  assertthat::assert_that(is.list(tail_chunk_list),
                          msg = paste0("Given tail_chunk_list is not a list (class). Please provide valid object."))
  assertthat::assert_that(is.list(predicted_list),
                          msg = paste0("Given predicted_list is not a list (class). Please provide valid object."))

  if (checkmate::test_string(nanopolish)) {
    # if string provided as an argument, read from file
    checkmate::assert_file_exists(nanopolish)
    nanopolish_polya_table <- vroom::vroom(nanopolish,
                                           col_select=c(readname, contig, polya_length, qc_tag),
                                           show_col_types = FALSE)
  } else {
    # make sure that nanopolish is an object with rows
    if (!is.data.frame(nanopolish) || nrow(nanopolish) == 0) {
      stop("Empty data frame provided as an input (nanopolish). Please provide valid input")
    }

    nanopolish_polya_table <- nanopolish[,c("readname", "contig","polya_length","qc_tag")]
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
  cat(paste0('[', as.character(Sys.time()), '] ','Retrieving estimated length data...', '\n', sep=''))

  #set progressbar
  pb <- utils::txtProgressBar(min = 0,
                              max = length(names(tail_feature_list[[1]])),
                              style = 3,
                              width = 50,
                              char = "=",
                              file = stderr())
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  tail_length_list <- foreach::foreach(i = seq_along(tail_feature_list[[1]]),
                                       .combine = c, .inorder = TRUE,
                                       .errorhandling = 'pass',
                                       .options.snow = opts,
                                       .options.multicore = mc_options) %dopar% {
                                         lapply(names(tail_feature_list[[1]][i]), function(x) length(tail_feature_list[[1]][[x]][[2]]))
                                       }

  #close(pb)

  #coerce to df, add names
  tail_length_list <- do.call("rbind.data.frame", tail_length_list)
  tail_length_list$readname <- names(tail_feature_list[["tail_feature_list"]])
  colnames(tail_length_list) <- c("signal_length", "readname")

  #merge data from feature list with nanopolish estimations
  tails_tail_feature_list <- dplyr::left_join(tail_length_list, nanopolish_polya_table, by="readname")


  ### Chunk positional data
  #create empty list for extracted data
  non_a_position_list <- list()

  # header for progress bar
  cat(paste0('[', as.character(Sys.time()), '] ','Retrieving position calibrating data...', '\n', sep=''))

  #set progressbar
  # progress bar
  pb <- utils::txtProgressBar(min = 0,
                              max = length(tail_chunk_list),
                              style = 3,
                              width = 50,
                              char = "=",
                              file = stderr())
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)


  non_a_position_list <- foreach::foreach(i = seq_along(tail_chunk_list),
                                          .combine = c, .inorder = TRUE,
                                          .errorhandling = 'pass',
                                          .options.snow = opts,
                                          .options.multicore = mc_options) %dopar% {
                                            lapply(tail_chunk_list[[i]], function(x) x[['chunk_start_pos']]+50)
                                          }
  #close(pb)

  # coerce to df, add names
  non_a_position_list <- do.call("rbind.data.frame", non_a_position_list)
  chunknames <-  purrr::map_depth(tail_chunk_list, 1, names) %>% unlist(use.names = F)
  non_a_position_list$chunkname <- chunknames
  non_a_position_list<- non_a_position_list %>% dplyr::mutate(readname = gsub('_.*','',chunkname))
  colnames(non_a_position_list)[1] <- c("centr_signal_pos")

  #merge data from feature list with nanopolish estimations
  non_a_position_list <- dplyr::left_join(non_a_position_list, tails_tail_feature_list, by="readname")

  # HANDLE PREDICTIONS
  moved_chunks_table <- data.frame(t(Reduce(rbind, predicted_list)))
  # rename cols
  colnames(moved_chunks_table) <- c("chunkname","prediction")
  # extract readnames
  moved_chunks_table$readname <- sub('\\_.*', '', moved_chunks_table$chunkname)

  #substitute predictions with letter code
  # moved_chunks_table$prediction[moved_chunks_table$prediction ==0] <- "A"
  # moved_chunks_table$prediction[moved_chunks_table$prediction ==1] <- "C"
  # moved_chunks_table$prediction[moved_chunks_table$prediction ==2] <- "G"
  # moved_chunks_table$prediction[moved_chunks_table$prediction ==3] <- "U"

  # vectorized substitution of predictions with letter code
  prediction_dict <- c("0" = "A", "1" = "C", "2" = "G", "3" = "U")
  moved_chunks_table$prediction <- prediction_dict[as.character(moved_chunks_table$prediction)]

  #extract reads with move==1 and modification absent (not detected)
  moved_blank_readnames <- names(which(with(moved_chunks_table, tapply(prediction, readname, unique) == 'A')))

  # cleaned chunks_table
  moved_chunks_table <- subset(moved_chunks_table, !(readname %in% moved_blank_readnames))
  # delete A-containing rows
  moved_chunks_table <- moved_chunks_table[!(moved_chunks_table$prediction=="A"),]
  #merge data from feats & predictions
  moved_chunks_table <- dplyr::left_join(moved_chunks_table, non_a_position_list, by=c("readname", "chunkname"))

  #estimate non-A nucleotide position
  moved_chunks_table$est_nonA_pos <- round(moved_chunks_table$polya_length-((moved_chunks_table$polya_length*moved_chunks_table$centr_signal_pos)/moved_chunks_table$signal_length), digits=2)

  #clean up the output nonA table:
  moved_chunks_table <- moved_chunks_table[,c(3,6,2,9,7,8)]

  # Handle other (discarded) reads:
  discarded_reads <- nanopolish_polya_table[!nanopolish_polya_table$readname %in% moved_chunks_table$readname,]

  # Add filtering criterion: select only pass or pass $ suffclip
  if(pass_only == TRUE){
    discarded_reads <- discarded_reads %>%
      dplyr::filter(!readname %in% moved_chunks_table$readname) %>%
      dplyr::mutate(comments = dplyr::case_when(polya_length < 10 ~ "IRL",
                                                qc_tag == "SUFFCLIP" ~ "NIN",
                                                qc_tag == "ADAPTER" ~ "QCF",
                                                qc_tag == "NOREGION" ~ "QCF",
                                                qc_tag == "READ_FAILED_LOAD" ~ "QCF",
                                                readname %in% moved_blank_readnames ~ "MPU",
                                                TRUE ~ "MAU"),
                    class = dplyr::case_when(polya_length < 10 ~ "unclassified",
                                             readname %in% moved_blank_readnames ~ "blank",
                                             comments == "MAU" ~ "blank",
                                             TRUE ~ "unclassified"))
  } else {
    discarded_reads <- discarded_reads %>%
      dplyr::filter(!readname %in% moved_chunks_table$readname) %>%
      dplyr::mutate(comments = dplyr::case_when(polya_length < 10 ~ "IRL",
                                                qc_tag == "ADAPTER" ~ "QCF",
                                                qc_tag == "NOREGION" ~ "QCF",
                                                qc_tag == "READ_FAILED_LOAD" ~ "QCF",
                                                readname %in% moved_blank_readnames ~ "MPU",
                                                TRUE ~ "MAU"),
                    class = dplyr::case_when(polya_length < 10 ~ "unclassified",
                                             readname %in% moved_blank_readnames ~ "blank",
                                             comments == "MAU" ~ "blank",
                                             TRUE ~ "unclassified"))
  }


  decorated_reads <- nanopolish_polya_table[nanopolish_polya_table$readname %in% moved_chunks_table$readname,]


  # decorated_reads <- decorated_reads %>% dplyr::mutate(class = "decorated",
  #                                                      comments = "YAY")
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


  if(qc == TRUE){
    # filter the outermost positions (termini!) from position data:
    # empirically tested constraints! report without terminal data
    moved_chunks_table_trimmed <- moved_chunks_table %>% dplyr::filter(!(est_nonA_pos < 2) & !(est_nonA_pos > polya_length-2))
    moved_chunks_table_discarded <- subset(moved_chunks_table,!(readname %in% moved_chunks_table_trimmed$readname))


    # subset potential artifacts from read_classes
    potential_artifacts <- subset(nanopolish_polya_table, readname %in% moved_chunks_table_discarded$readname)

    decorated_reads_edited <- nanopolish_polya_table %>%
      dplyr::mutate(class = dplyr::case_when(
        readname %in% potential_artifacts$readname ~ paste0(class,"-WARN"),
        TRUE ~ paste0(class)))


    # label potential artifacts in nonadenosine residue dataframe
    moved_chunks_table_qc <- moved_chunks_table %>%
      dplyr::mutate(prediction=dplyr::case_when(est_nonA_pos < 2 ~ paste0(prediction, "-WARN"),
                                                est_nonA_pos > polya_length-2 ~ paste0(prediction, "-WARN"),
                                                TRUE~ paste0(prediction)))


    #CREATE FINAL OUTPUT
    #prevent potential bugs inherited from nanopolish multimapping
    decorated_reads_edited <- unique(decorated_reads_edited)
    moved_chunks_table_qc <- unique(moved_chunks_table_qc)

    ninetails_output[['read_classes']] <- decorated_reads_edited
    ninetails_output[['nonadenosine_residues']] <- moved_chunks_table_qc

  } else{
    #CREATE FINAL OUTPUT
    #prevent potential bugs inherited from nanopolish multimapping
    nanopolish_polya_table <- unique(nanopolish_polya_table)
    moved_chunks_table <- unique(moved_chunks_table)

    ninetails_output[['read_classes']] <- nanopolish_polya_table
    ninetails_output[['nonadenosine_residues']] <- moved_chunks_table
  }

  return(ninetails_output)

}



