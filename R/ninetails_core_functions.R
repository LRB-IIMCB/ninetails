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
#'
#'
#' @examples
#'\dontrun{
#'
#' extract_polya_data(nanopolish = '/path/to/nanopolish/polya/output.tsv',
#'                    sequencing_summary = '/path/to/sequencing_summary.txt',
#'                    pass_only = TRUE)
#'
#'}
#'
#'

extract_polya_data <- function(nanopolish, sequencing_summary, pass_only = TRUE){

  # variable binding (suppressing R CMD check from throwing an error)
  readname <- polya_start <- transcript_start <- polya_length <- qc_tag  <- filename <- read_id <- NULL

    if (missing(nanopolish)) {
    stop("Nanopolish polya output is missing. Please provide a valid nanopolish argument.", .call = FALSE)
  }

  if (missing(sequencing_summary)) {
    stop("Sequencing summary file is missing. Please provide a valid sequencing_summary argument.", .call = FALSE)
  }

  assertthat::assert_that(assertive::is_a_bool(pass_only),
                          msg="Please provide TRUE/FALSE values for pass_only parameter")

  if (assertive::is_character(nanopolish)) {

    assertthat::assert_that(assertive::is_existing_file(nanopolish),
                            msg=paste0("File ",nanopolish," does not exist",sep=""))

    nanopolish_polya_table <- vroom::vroom(nanopolish,
                                           col_select=c(readname, polya_start, transcript_start, polya_length, qc_tag),
                                           show_col_types = FALSE)
  }
  else if (assertive::has_rows(nanopolish)) {
    nanopolish_polya_table <- nanopolish[,c("readname","polya_start","transcript_start","polya_length","qc_tag")]
  }
  else {
    stop("Wrong nanopolish parameter. Please provide filepath or object.")
  }

  #in case if smth goes wrong with rows
  assertthat::assert_that(assertive::has_rows(nanopolish_polya_table),
                          msg = "Empty data frame provided as an input (nanopolish). Please provide valid input")


  if (assertive::is_character(sequencing_summary)) {

    assertthat::assert_that(assertive::is_existing_file(sequencing_summary),
                            msg=paste0("File ",sequencing_summary," does not exist",sep=""))

    sequencing_summary_table <- vroom::vroom(sequencing_summary,
                                             col_select = c(filename, read_id),
                                             show_col_types = FALSE)
  }
  else if (assertive::has_rows(sequencing_summary)) {
    sequencing_summary_table <- sequencing_summary
  }
  else {
    stop("Wrong sequencing_summary parameter. Please provide filepath or object.")
  }

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
#' extract_tail_data(readname = 'abc123de-fg45-6789-0987-6543hijk2109',
#'                  polya_summary = polya_summary_table,
#'                  workspace = '/path/to/folder/containing/multifast5s',
#'                  basecall_group = 'Basecall_1D_000')
#'
#'}
#'
#'

extract_tail_data <- function(readname, polya_summary, workspace, basecall_group){

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

  assertthat::assert_that(assertive::has_rows(polya_summary),
                          msg = paste0("Empty data frame provided as an input (polya_summary). Please provide a valid input table."))
  assertthat::assert_that(assertive::is_a_non_missing_nor_empty_string(workspace),
                          msg = paste0("Empty string provided as an input. Please provide a valid path to basecalled fast5 files."))
  assertthat::assert_that(assertive::is_character(workspace),
                          msg = paste0("Path to basecalled fast5 files is not a character string. Please provide a valid path to basecalled fast5 files."))
  assertthat::assert_that(assertive::is_character(readname),
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
  signal <- signal[polya_start_position:polya_end_position]

  # downsample (interpolate signal and moves) so the computation would be faster/easier to handle.
  signal <- round(stats::approx(signal, method = "linear", n=ceiling(0.2 * length(signal)))[[2]], digits=0)
  moves_tail_range <- round(stats::approx(moves_tail_range, method = "linear", n=ceiling(0.2 * length(moves_tail_range)))[[2]], digits=0)


  extracted_data_single_list = list() # creating empty list for the extracted fast5 data

  extracted_data_single_list[["fast5_filename"]] <- polya_summary$filename[read_idx]
  extracted_data_single_list[["polya_start_pos"]] <- polya_start_position
  extracted_data_single_list[["transcript_start_pos"]] <- transcript_start_position
  extracted_data_single_list[["tail_signal"]] <- signal
  extracted_data_single_list[["tail_moves"]] <- moves_tail_range
  extracted_data_single_list[["tail_length"]] <- polya_summary$polya_length[read_idx]


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
#' create_tail_feature_list(nanopolish = '/path/to/file',
#'                          sequencing_summary = '/path/to/file',
#'                          workspace = '/path/to/guppy/workspace',
#'                          num_cores = 10,
#'                          basecall_group = 'Basecall_1D_000')
#'
#'}
#'
create_tail_feature_list <- function(nanopolish, sequencing_summary, workspace, num_cores, basecall_group, pass_only=TRUE){

  # variable binding (suppressing R CMD check from throwing an error)
  nam <- NULL

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

  assertthat::assert_that(assertive::is_numeric(num_cores),
                          msg=paste0("Declared core number must be numeric. Please provide a valid argument."))

  # Extracting and processing polya & sequencing summary data
  polya_summary <- extract_polya_data(nanopolish, sequencing_summary, pass_only)

  # this is list of indexes required for parallel computing; the main list of reads is split for chunks (does not accept readname)
  index_list = split(1:length(names(polya_summary$filename)),
                     ceiling(1:length(names(polya_summary$filename))/100))

  #create empty list for extracted fast5 data
  tail_feature_list = list()

  # header for progress bar
  cat(paste0('[', as.character(Sys.time()), '] ','Extracting features of provided reads...', '\n', sep=''))

  # progress bar
  pb <- utils::txtProgressBar(min = 0,
                              max = length(index_list),
                              style = 3,
                              width = 50,
                              char = "=",
                              file= stderr())

  # use selected number of cores
  doParallel::registerDoParallel(cores = num_cores)

  # loop for parallel extraction
  for (indx in 1:length(index_list)){

    # work on subsets of reads in parallel
    tail_feature_list <- c(tail_feature_list, foreach::foreach(nam = names(polya_summary$filename)[index_list[[indx]]]) %dopar% extract_tail_data(nam,polya_summary,workspace,basecall_group))

    utils::setTxtProgressBar(pb, indx)

  }

  close(pb)

  #stop cluster
  doParallel::stopImplicitCluster()

  #label each signal according to corresponding read name to avoid confusion
  squiggle_names <- polya_summary$readname
  names(tail_feature_list) <- squiggle_names

  # remove reads with only zero moved tails
  tail_feature_list <- Filter(function(x) sum(x$tail_moves) !=0, tail_feature_list)

  # Done comm
  cat(paste0('[', as.character(Sys.time()), '] ','Done!', '\n', sep=''))

  return(tail_feature_list)

}



#' Creates vector of overlapping tail fragments extracted from
#' nanopolish output and fast5 files, in which move value = 1 occurs.
#'
#'
#' @param readname character string. Name of the given read within the
#' analyzed dataset.
#'
#' @param tail_feature_list list object produced by create_tail_feature_list
#' function.
#'
#' @param segment numeric [1]. Length of the chunk(s) to be created.
#'
#' @param overlap numeric [1]. Length of the overlap between the chunks.
#'
#' @return A vector of indices of moved chunks organized by their order in
#' an original list is returned.
#'
#' @export
#'
#' @examples
#'\dontrun{
#'
#' split_with_overlaps_moved(readname = "name_of_given_read_of_interest",
#'                           tail_feature_list = list_of_tail_features,
#'                           segment = 100
#'                           overlap = 10)
#'
#'}
#'
#'
split_with_overlaps_moved <- function(readname, tail_feature_list, segment, overlap) {
  # avoiding 'no visible binding for global variable' error
  x <- NULL

  #assertions
  if (missing(readname)) {
    stop("Readname is missing. Please provide a valid readname argument.", call. =FALSE)
  }

  if (missing(tail_feature_list)) {
    stop("List of tail features is missing. Please provide a valid tail_feature_list argument.", call. =FALSE)
  }

  assertthat::assert_that(assertive::is_character(readname),
                          msg = paste0("Given readname is not a character string. Please provide a valid readname."))
  assertthat::assert_that(assertive::is_list(tail_feature_list),
                          msg = paste0("Given tail_feature_list is not a list (class). Please provide valid file format."))

  checkmate::assert_integerish(segment,
                               tol = sqrt(.Machine$double.eps),
                               lower = 1,
                               upper = Inf,
                               any.missing = TRUE,
                               all.missing = TRUE,
                               len = NULL,
                               min.len = 1L,
                               max.len = 1,
                               unique = FALSE,
                               sorted = FALSE,
                               names = NULL,
                               typed.missing = FALSE,
                               null.ok = FALSE,
                               coerce = FALSE,
                               .var.name = checkmate::vname(x),
                               add = NULL)

  checkmate::assert_integerish(overlap,
                               tol = sqrt(.Machine$double.eps),
                               lower = 0,
                               upper = Inf,
                               any.missing = TRUE,
                               all.missing = TRUE,
                               len = NULL,
                               min.len = 1L,
                               max.len = 1,
                               unique = FALSE,
                               sorted = FALSE,
                               names = NULL,
                               typed.missing = FALSE,
                               null.ok = FALSE,
                               coerce = FALSE,
                               .var.name = checkmate::vname(x),
                               add = NULL)


  assertthat::assert_that(overlap<segment, msg="Segment value must be greater than overlap value. Please provide a valid argument.")


  signal <- tail_feature_list[[readname]][[4]]
  moves <- tail_feature_list[[readname]][[5]]

  #prevent from running the function on zero-only vectors
  if (sum(moves) == 0){
    stop("The move vector does not contain the value of 1. Please provide a valid move argument (use create_tail_chunk_list_moved() function or prefilter data manually).", call. =FALSE)
  }

  #initial coordinates (for all chunks)
  start_coordinates_total <- seq(1, length(moves), by=segment-overlap)
  end_coordinates_total   <- start_coordinates_total + segment - 1

  #extract indices of "moved" chunks
  moved_chunks_indices <- lapply(1:length(start_coordinates_total),
                                 function(i) 1 %in% moves[start_coordinates_total[i]:end_coordinates_total[i]])
  moved_chunks_indices <- which(unlist(moved_chunks_indices))

  #coordinates of selected "moved" chunks
  start_coordinates_selected <- start_coordinates_total[moved_chunks_indices]
  end_coordinates_selected <- start_coordinates_selected + segment - 1

  #extract ONLY "moved" signal chunks
  extracted_moved_signals <- lapply(1:length(start_coordinates_selected),
                                    function(i) signal[start_coordinates_selected[i]:end_coordinates_selected[i]])
  # replace NAs with 3 most frequent values (randomly sampled)
  # if all values would be equal, so the breaks would not be unique

  most_freq_vals <- as.numeric(names(sort(table(signal),decreasing=TRUE)[1:3]))
  extracted_moved_signals <- lapply(extracted_moved_signals,
                                    function(n) replace(n, is.na(n), sample(most_freq_vals,sum(is.na(n)),TRUE)))

  # naming chunks based on readnames & indices
  chunk_names <- paste0(rep(readname, length(extracted_moved_signals)), '_', unlist(moved_chunks_indices))

  #name_chunks
  names(extracted_moved_signals) <- chunk_names

  return(extracted_moved_signals)

}


#' Creates list of overlapping tail fragments extracted from
#' nanopolish output and fast5 files based on provided feature list.
#'
#' Only fragments containing move==1 are included
#' (as create_tail_feature_list function performs prefiltering according to this criterion).
#'
#' @param tail_feature_list list object produced by create_tail_feature_list
#' function.
#'
#' @param num_cores numeric [1]. Number of physical cores to use in processing
#' the data. Do not exceed 1 less than the number of cores at your disposal.
#'
#' @return A list of tail range chunks organized by the read ID
#' is returned. Only tails containing at least one move value equal to 1
#' are included. Always assign this returned list to a variable, otherwise
#' the long list will be printed to the console, which may crash your R session.
#'
#' @importFrom foreach %dopar%
#'
#' @export
#'
#' @examples
#'\dontrun{
#'
#' create_tail_chunk_list_moved(tail_feature_list = '/path/to/tail_feature_list', num_cores = 10)
#'
#'}
#'
#'
create_tail_chunk_list_moved <- function(tail_feature_list, num_cores){

  # variable binding (suppressing R CMD check from throwing an error)
  nam <- NULL

  # initial assertions
  if (missing(num_cores)) {
    stop("Number of declared cores is missing. Please provide a valid num_cores argument.", call. =FALSE)
  }

  if (missing(tail_feature_list)) {
    stop("List of features is missing. Please provide a valid tail_feature_list argument.", call. =FALSE)
  }

  assertthat::assert_that(assertive::is_numeric(num_cores),
                          msg = paste0("Declared core number must be numeric. Please provide a valid argument."))
  assertthat::assert_that(assertive::is_list(tail_feature_list),
                          msg = paste0("Given tail_feature_list is not a list (class). Please provide valid file format."))

  # creating cluster for parallel computing
  doParallel::registerDoParallel(cores = num_cores)

  # this is list of indexes required for parallel computing; the main list of reads is split for chunks
  index_list = split(1:length(names(tail_feature_list)), ceiling(1:length(names(tail_feature_list))/100))

  # header for progress bar
  cat(paste0('[', as.character(Sys.time()), '] ','Creating tail segmentation data...', '\n', sep=''))

  # progress bar
  pb <- utils::txtProgressBar(min = 0,
                              max = length(index_list),
                              style = 3,
                              width = 50,
                              char = "=",
                              file= stderr())

  #create empty list for extracted data
  tail_chunk_list = list()

  # loop for parallel extraction
  for (indx in 1:length(index_list)){

    # work on subsets of reads in parallel
    tail_chunk_list <- c(tail_chunk_list, foreach::foreach(nam = names(tail_feature_list)[index_list[[indx]]])
                         %dopar% split_with_overlaps_moved(nam, tail_feature_list, segment = 100, overlap = 50))

    utils::setTxtProgressBar(pb, indx)

  }

  close(pb)

  #stop cluster
  doParallel::stopImplicitCluster()

  #rename first level of nested list accordingly
  names(tail_chunk_list) <- names(tail_feature_list)

  # Done comm
  cat(paste0('[', as.character(Sys.time()), '] ','Done!', '\n', sep=''))

  return(tail_chunk_list)

}


#' Converting ONT signal to gramian angular summation field.
#'
#' @param tail_chunk character string. Name of the given signal chunk within the
#' analyzed dataset.
#'
#' @return an array (100,100,1) with values representing ONT signal.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' create_gasf(tail_chunk = "1234-jhgfdr54-io643gju-ouy4378989")
#'
#'}
#'
create_gasf <- function(tail_chunk){
  x <- NULL

  #assertions
  if (missing(tail_chunk)) {
    stop("Tail_chunk is missing. Please provide a valid tail_chunk argument.", call. =FALSE)
  }

  checkmate::assert_integerish(tail_chunk,
                               tol = sqrt(.Machine$double.eps),
                               lower = -Inf,
                               upper = Inf,
                               any.missing = TRUE,
                               all.missing = TRUE,
                               len = NULL,
                               min.len = 100L,
                               max.len = 100,
                               unique = FALSE,
                               sorted = FALSE,
                               names = NULL,
                               typed.missing = FALSE,
                               null.ok = FALSE,
                               coerce = FALSE,
                               .var.name = checkmate::vname(x),
                               add = NULL)



  # rescale values so that all of them fall in the interval [-1, 1]:
  tail_chunk <- (tail_chunk-max(tail_chunk)+(tail_chunk-min(tail_chunk)))/(max(tail_chunk)-min(tail_chunk))

  # calculate phi coefficient for interpolation to polar coordinates
  tail_chunk <- acos(tail_chunk)

  # create matrix by replicating vec
  tail_chunk <- cbind(replicate(length(tail_chunk), tail_chunk))

  #calculate sum of phi
  tail_chunk <- tail_chunk + t(tail_chunk)

  #calculate cosinus
  tail_chunk <- cos(tail_chunk)

  #reshape the data into new dimensions
  tail_chunk <- array(t(tail_chunk), c(100,100,1))

  # rescale values so that all of them fall in the interval [0, 1]:
  tail_chunk <- round((tail_chunk-min(tail_chunk))/(max(tail_chunk)-min(tail_chunk)), 4)


  return(tail_chunk)
}


#' Creates list of gramian angular summation matrices produced based on
#' list of splitted tails (tail chunks).
#'
#' @param tail_chunk_list character string. The list object produced
#' by create_chunk_list function.
#'
#' @param num_cores numeric [1]. Number of physical cores to use in processing
#' the data. Do not exceed 1 less than the number of cores at your disposal.
#'
#' @return A list of gasf matrices organized by the read ID_index
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
#' create_gasf_list(feature_list = '/path/to/chunk_list', num_cores = 10)
#'
#'}
#'
#'
create_gasf_list <- function(tail_chunk_list, num_cores){

  # Assertions
  if (missing(num_cores)) {
    stop("Number of declared cores is missing. Please provide a valid num_cores argument.", call. =FALSE)
  }

  if (missing(tail_chunk_list)) {
    stop("List of tail chunks is missing. Please provide a valid tail_chunk_list argument.", call. =FALSE)
  }

  assertthat::assert_that(assertive::is_list(tail_chunk_list),
                          msg = paste0("Given tail_chunk_list is not a list (class). Please provide valid file format."))
  assertthat::assert_that(assertive::is_numeric(num_cores),
                          msg=paste0("Declared core number must be numeric. Please provide a valid argument."))


  #register cores for parallelization
  doParallel::registerDoParallel(cores = num_cores)

  #retrieve chunknames
  chunknames <- gsub(".*?\\.","",names(rapply(tail_chunk_list, function(x) head(x, 1))))

  #create empty list for the data
  gasf_list = list()

  #set progressbar
  cat(paste0('[', as.character(Sys.time()), '] ','Computing gramian angular summation fields...', '\n', sep=''))
  pb <- utils::txtProgressBar(min = 0,
                              max = length(tail_chunk_list),
                              style = 3,
                              width = 50,
                              char = "=",
                              file= stderr())

  #loop through the nested list
  for (read in seq_along(tail_chunk_list)){
    for (chunk in seq_along(tail_chunk_list[[read]])){
      gasf_list <- c(gasf_list, foreach::foreach(chunk) %dopar% create_gasf(tail_chunk_list[[read]][[chunk]]))

      utils::setTxtProgressBar(pb, read)
    }
  }

  close(pb)

  #stop cluster
  doParallel::stopImplicitCluster()

  #restore names in the list
  names(gasf_list) <- chunknames

  # Done comm
  cat(paste0('[', as.character(Sys.time()), '] ','Done!', '\n', sep=''))

  return(gasf_list)

}

#' Classification of matrices produced from read chunks with CNN.
#'
#' @param gasf_list [list] A list of gasf matrices organized by the read ID_index.
#'
#' @return a list of gasfs predictions based on used model.
#'
#'
#' @export
#'
#' @examples
#'\dontrun{
#'
#' predict_classes(gasf_list = gasf_list, keras_model = "/path/to/the_model.h5")
#'
#'}
predict_classes <- function(gasf_list){

  #assertions
  if (missing(gasf_list)) {
    stop("List of transformed signal chunks is missing. Please provide a valid gasf_list argument.", call. =FALSE)
  }




  chunknames <- names(gasf_list)
  names(gasf_list)  <- NULL

  gasf_list <- simplify2array(gasf_list)
  gasf_list <- aperm(gasf_list, c(4,1,2,3))

  # Output info
  cat(paste0('[', as.character(Sys.time()), '] ','Classifying gramian angular summation fields...', '\n', sep=''))


  keras_model <- load_keras_model()

  #predict chunk class
  predicted_gasf_classes <- keras_model %>% stats::predict(gasf_list) %>% keras::k_argmax()
  predicted_gasf_classes <- as.numeric(predicted_gasf_classes)

  predicted_list = list() # creating empty list for the extracted  data

  predicted_list[["chunkname"]] <- chunknames
  predicted_list[["prediction"]] <- predicted_gasf_classes

  # Output info
  cat(paste0('[', as.character(Sys.time()), '] ','Done!', '\n', sep=''))


  return(predicted_list)
}

