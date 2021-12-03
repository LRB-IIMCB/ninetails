#' Extracts features of single RNA read from respective basecalled multi-fast5 file.
#'
#' This function extracts metadata of RNA read from multi-fast5 file basecalled
#' by guppy. As default. Filenames are taken from the sequencing
#' summary file.
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
#' extract_squiggle_data(readname = 'abc123de-fg45-6789-0987-6543hijk2109',
#'                       polya_summary = polya_summary_table,
#'                       workspace = '/path/to/folder/containing/multifast5s',
#'                       basecall_group = 'Basecall_1D_000')
#'
#'}
#'
#'

# EXTRACTING READ METADATA
extract_squiggle_data <- function(readname, polya_summary, workspace, basecall_group){


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

  assertthat::assert_that(assertive::has_rows(polya_summary),msg = "Empty data frame provided as an input (polya_summary). Please provide valid input table.")
  assertthat::assert_that(assertive::is_character(workspace), msg = "Path to basecalled fast5 files is not a character string. Please provide valid path to basecalled fast5 files.")

  # Extract data from fast5 file
  fast5_filenames <- polya_summary$filename # this is the list with file names & corresponding readnames (named vec)
  fast5_readname <- paste0("read_",readname) # in fast5 file structure each readname has suffix "read_"
  fast5_file <- fast5_filenames[readname] # this is a particular file name (eg. FAO12345.fast5)
  fast5_file_path <-file.path(workspace, fast5_file) #this is path to browsed fast5 file


  raw_signal <- rhdf5::h5read(file.path(fast5_file_path),paste0(fast5_readname,"/Raw/Signal"))
  extracted_raw <- as.vector(raw_signal) # initially signal is stored as an array, I prefer to vectorize it for potential further manipulations with dframes


  # close all handled instances (groups, attrs) of fast5 file
  rhdf5::h5closeAll()

  # Extract data from polya summary
  #extract polya start and end positions from nanopolish data (start and end positions of polyA tail predicted by nanopolish polya function)
  read_id <- which(polya_summary$readname == readname) # this is to retrieve row number because tibbles cannot have rownames

  extracted_data_single_list = list() # creating empty list for the extracted fast5 data

  extracted_data_single_list[["fast5_filename"]] <- polya_summary$filename[read_id]
  extracted_data_single_list[["polya_start_pos"]] <- polya_summary$polya_start[read_id]
  extracted_data_single_list[["transcript_start_pos"]] <- polya_summary$transcript_start[read_id]
  extracted_data_single_list[["raw_signal"]] <- extracted_raw


  return(extracted_data_single_list)

}
