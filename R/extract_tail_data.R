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
#'                       polya_summary = polya_summary_table,
#'                       workspace = '/path/to/folder/containing/multifast5s',
#'                       basecall_group = 'Basecall_1D_000')
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
                          msg = "Empty data frame provided as an input (polya_summary). Please provide a valid input table.")
  assertthat::assert_that(assertive::is_a_non_missing_nor_empty_string(workspace),
                          msg = "Empty string provided as an input. Please provide a valid path to basecalled fast5 files.")
  assertthat::assert_that(assertive::is_character(workspace),
                          msg = "Path to basecalled fast5 files is not a character string. Please provide a valid path to basecalled fast5 files.")
  assertthat::assert_that(assertive::is_character(readname),
                          msg = "Given readname is not a character string. Please provide a valid readname argument.")

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
  read_id <- which(polya_summary$readname == readname) # this is to retrieve row number because tibbles cannot have rownames
  polya_start_position <- polya_summary$polya_start[read_id]
  transcript_start_position <- polya_summary$transcript_start[read_id]
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

  extracted_data_single_list[["fast5_filename"]] <- polya_summary$filename[read_id]
  extracted_data_single_list[["polya_start_pos"]] <- polya_start_position
  extracted_data_single_list[["transcript_start_pos"]] <- transcript_start_position
  extracted_data_single_list[["tail_signal"]] <- signal
  extracted_data_single_list[["tail_moves"]] <- moves_tail_range


  return(extracted_data_single_list)

}
