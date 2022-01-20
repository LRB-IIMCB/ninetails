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
#' @export
#'
#' @examples
#'\dontrun{
#'
#' create_tail_feature_list(nanopolish = '/path/to/file',
#'                          sequencing_summary = '/path/to/file',
#'                          workspace = '/path/to/guppy/workspace',
#'                          num_cores = 10,
#'                          basecalled_group = 'Basecall_1D_000')
#'
#'}
#'
#'

create_tail_feature_list <- function(nanopolish, sequencing_summary, workspace, num_cores, basecall_group, pass_only=TRUE){

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

  assertthat::assert_that(assertive::is_numeric(num_cores), msg=paste("Declared core number must be numeric. Please provide a valid argument."))

  # Extracting and processing polya & sequencing summary data
  polya_summary <- extract_polya_data(nanopolish, sequencing_summary, pass_only)

  # this is list of indexes required for parallel computing; the main list of reads is split for chunks (does not accept readname)
  index_list = split(1:length(names(polya_summary$filename)), ceiling(1:length(names(polya_summary$filename))/100))


  #checking data format
  check_fast5_filetype(workspace, basecall_group)


  #create empty list for extracted fast5 data
  extracted_data_multiple_list = list()


  # header for progress bar
  cat(paste('Extracting features of provided reads...', '\n', sep=''))

  # progress bar
  pb <- utils::txtProgressBar(min = 0, max = length(index_list), style = 3, width = 50, char = "=")


  # loop for parallel extraction
  for (indx in 1:length(index_list)){

    # use selected number of cores
    doParallel::registerDoParallel(cores = num_cores)

    # work on subsets of reads in parallel
    extracted_data_multiple_list <- c(extracted_data_multiple_list, foreach::foreach(nam = names(polya_summary$filename)[index_list[[indx]]]) %dopar% extract_tail_data(nam,polya_summary,workspace,basecall_group))

    utils::setTxtProgressBar(pb, indx)

  }

  close(pb)

  #label each signal according to corresponding read name to avoid confusion
  squiggle_names <- polya_summary$readname
  names(extracted_data_multiple_list) <- squiggle_names


  return(extracted_data_multiple_list)

}
