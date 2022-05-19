#' Wrapper function for ninetails package.
#'
#' This function allows to perform all of the steps required to discover
#' nonadenosine nucleotides within the given dataset using ninetails.
#' Please keep in mind, that during computations the function creates large
#' segmentation data. Therefore it may be wise to split the nanopolish table
#' beforehand and then run this function on nanopolish table chunks.
#'
#' The output of this function is a list of 3 dataframes, containing: a) binary
#' classification of reads satisfying filtering criteria*, b) detailed positional
#' info regarding all potential nonadenosine residues detected, c) info
#' regarding discarded reads (the reason why the following read has been omitted
#' from the analysis).
#'
#' The filtering criteria applied by ninetails are as follows: move of value =1
#' present, qc_tag = "PASS" or "PASS" & "SUFFCLIP", length estimated
#' by nanopolish >=10 nt).
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
#' @param save_dir character string. Full path of the directory where the output
#' files containing the tail composition information should be stored.
#'
#' @export
#'
#' @return A list containing tail information organized by the read ID
#' is returned. Always assign this returned list to a variable, otherwise
#' the long list will be printed to the console, which may crash your R session.
#'
#' @importFrom foreach %dopar%
#' @importFrom utils head
#'
#' @examples
#'\dontrun{
#'
#' check_tails(nanopolish = '/path/to/file',
#'             sequencing_summary = '/path/to/file',
#'             workspace = '/path/to/guppy/workspace',
#'             num_cores = 10,
#'             basecalled_group = 'Basecall_1D_000',
#'             pass_only=TRUE,
#'             save_dir = '/directory/where/output/shall/be/stored')
#'
#' }

check_tails <- function(nanopolish, sequencing_summary, workspace, num_cores, basecall_group, pass_only, save_dir){

  # variable binding (suppressing R CMD check from throwing an error)
  nam <- NULL

  #Show console message
  cat(paste0('Welcome to Ninetails ', as.character(utils::packageVersion("ninetails")), '\n', '\n', sep = ""))

  # Create a log file
  if (dir.exists(file.path(save_dir))) {
    log_filename <- paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_ninetails.log", sep = "")
    log_filepath <- file.path(save_dir, log_filename, fsep = .Platform$file.sep)
    log_file <- file(log_filepath, open = "a")
    sink(log_file, append=TRUE, split = TRUE, type='output')
    on.exit(sink(file=NULL, type = 'output'))
  }


  # Handle if save dir does not exist
  if (!dir.exists(file.path(save_dir))) {
    cat(paste0('[', as.character(Sys.time()), '] ', 'Defined save directory does not exist. Creating...','\n', sep=''))
    tryCatch({dir.create(file.path(save_dir, fsep = .Platform$file.sep))
      cat('Done!\n')
    },
    error=function(e){
      cat(paste0('[', as.character(Sys.time()), '] ', 'Failed to create the save dir. Results will be stored in the current working directory.\n', sep=''))
      save_dir <- getwd()
    })
    log_filename <- paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_ninetails.log", sep = "")
    log_filepath <- file.path(save_dir, log_filename, fsep = .Platform$file.sep)
    log_file <- file(log_filepath, open = "a")
    sink(log_file, append=TRUE, split = TRUE, type='output')
    on.exit(sink(file=NULL, type = 'output'))
  }


  # user-specified options
  cat(paste(' Ninetails was launched with following options:', '\n', sep=''))
  cat(paste(' nanopolish output file:       ', nanopolish, '\n', sep=''))
  cat(paste(' sequencing sumary file:       ', sequencing_summary, '\n', sep=''))
  cat(paste(' fast5 files directory:        ', workspace, '\n', sep=''))
  cat(paste(' number of cores:              ', num_cores, '\n', sep=''))
  cat(paste(' basecall group:               ', basecall_group, '\n', sep=''))
  cat(paste(' only "PASS" reads included:   ', pass_only, '\n', sep=''))
  cat(paste(' output directory:             ', save_dir, '\n', '\n', '\n'))

  cat(paste0('[', as.character(Sys.time()), '] ', 'Pipeline initialized','\n','\n'))


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

  if (missing(save_dir)) {
    stop("A save dir for the output files is missing. Please provide a valid save_dir argument.", call. =FALSE)
  }

  assertthat::assert_that(assertive::is_numeric(num_cores), msg=paste0("Declared core number must be numeric. Please provide a valid argument."))
  assertthat::assert_that(assertive::is_a_non_missing_nor_empty_string(workspace),
                          msg = paste0("Empty string provided as an input. Please provide a valid path to basecalled fast5 files."))
  assertthat::assert_that(assertive::is_character(workspace),
                          msg = paste0("Path to basecalled fast5 files is not a character string. Please provide a valid path to basecalled fast5 files."))
  assertthat::assert_that(assertive::is_character(save_dir),
                          msg = paste0("Path to output files is not a character string. Please provide a valid save_dir path."))


  #####################################################
  # CREATE TAIL FEATURE LIST
  #####################################################
  # Extracting and processing polya & sequencing summary data
  polya_summary <- extract_polya_data(nanopolish, sequencing_summary, pass_only)

  # this is list of indexes required for parallel computing; the main list of reads is split for chunks (does not accept readname)
  index_list = split(1:length(names(polya_summary$filename)), ceiling(1:length(names(polya_summary$filename))/100))

  #checking data format
  check_fast5_filetype(workspace, basecall_group)

  # use selected number of cores
  doParallel::registerDoParallel(cores = num_cores)

  #create empty list for extracted fast5 data
  tail_feature_list = list()

  # header for progress bar
  cat(paste0('[', as.character(Sys.time()), '] ','Extracting features of provided reads...', '\n', sep=''))

  #SINK #1
  sink(file=NULL, type = 'output')
  close(log_file)

  # progress bar
  pb <- utils::txtProgressBar(min = 0, max = length(index_list), style = 3, width = 50, char = "=")


  # loop for parallel extraction
  for (indx in 1:length(index_list)){

    # work on subsets of reads in parallel
    tail_feature_list <- c(tail_feature_list, foreach::foreach(nam = names(polya_summary$filename)[index_list[[indx]]]) %dopar% extract_tail_data(nam,polya_summary,workspace,basecall_group))

    utils::setTxtProgressBar(pb, indx)

  }

  close(pb)

  #SINK #2
  log_file <- file(log_filepath, open = "a")
  sink(log_file, append=TRUE, split = TRUE, type='output')

  #label each signal according to corresponding read name to avoid confusion
  squiggle_names <- polya_summary$readname
  names(tail_feature_list) <- squiggle_names

  # remove reads with only zero moved tails
  tail_feature_list <- Filter(function(x) sum(x$tail_moves) !=0, tail_feature_list)

  # Done comm
  cat(paste0('[', as.character(Sys.time()), '] ','Done!', '\n', sep=''))

  #####################################################
  # CREATE TAIL CHUNK LIST
  #####################################################


  # this is list of indexes required for parallel computing; the main list of reads is split for chunks
  index_list = split(1:length(names(tail_feature_list)), ceiling(1:length(names(tail_feature_list))/100))

  # header for progress bar
  cat(paste0('[', as.character(Sys.time()),'] ','Creating tail segmentation data...', '\n', sep=''))

  #SINK #1
  sink(file=NULL, type = 'output')
  close(log_file)


  # progress bar
  pb <- utils::txtProgressBar(min = 0, max = length(index_list), style = 3, width = 50, char = "=")

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

  #SINK #2
  log_file <- file(log_filepath, open = "a")
  sink(log_file, append=TRUE, split = TRUE, type='output')


  #rename first level of nested list accordingly
  names(tail_chunk_list) <- names(tail_feature_list)

  # Done comm
  cat(paste0('[', as.character(Sys.time()), '] ','Done!', '\n', sep=''))


  #####################################################
  # CREATE GASF LIST
  #####################################################

  #retrieve chunknames
  chunknames <- gsub(".*?\\.","",names(rapply(tail_chunk_list, function(x) head(x, 1))))

  #create empty list for the data
  gasf_list = list()

  #header for progressbar
  cat(paste0('[', as.character(Sys.time()), '] ','Computing gramian angular summation fields...', '\n', sep=''))

  #SINK #1
  sink(file=NULL, type = 'output')
  close(log_file)

  #set progressbar
  pb <- utils::txtProgressBar(min = 0, max = length(tail_chunk_list), style = 3, width = 50, char = "=")

  #loop through the nested list
  for (read in seq_along(tail_chunk_list)){
    for (chunk in seq_along(tail_chunk_list[[read]])){
      gasf_list <- c(gasf_list, foreach::foreach(chunk) %dopar% create_gasf(tail_chunk_list[[read]][[chunk]]))

      utils::setTxtProgressBar(pb, read)
    }
  }

  close(pb)

  #SINK #2
  log_file <- file(log_filepath, open = "a")
  sink(log_file, append=TRUE, split = TRUE, type='output')


  #restore names in the list
  names(gasf_list) <- chunknames

  #stop cluster
  doParallel::stopImplicitCluster()

  # Done comm
  cat(paste0('[', as.character(Sys.time()), '] ','Done!', '\n', sep=''))

  #####################################################
  # PREDICT CLASSES
  #####################################################

  predicted_list <- tryCatch({predict_classes(gasf_list)}) #suppress tensorflow console messages


  #####################################################
  # CREATE COORD DF
  #####################################################

  # creating cluster for parallel computing
  doParallel::registerDoParallel(cores = num_cores)


  # OVERLAP COUNT LIST
  # header for progress bar
  cat(paste0('[', as.character(Sys.time()), '] ','Retrieving segmentation data...', '\n', sep=''))

  #SINK #1
  sink(file=NULL, type = 'output')
  close(log_file)

  # progress bar
  pb <- utils::txtProgressBar(min = 0, max = length(index_list), style = 3, width = 50, char = "=")

  #create empty lists for extracted data
  overlap_count_list = list()

  # loop for parallel extraction
  for (indx in 1:length(index_list)){

    # work on subsets of reads in parallel
    overlap_count_list <- c(overlap_count_list, foreach::foreach(nam = names(tail_feature_list)[index_list[[indx]]])
                            %dopar% count_overlaps(signal = tail_feature_list[[nam]][[4]], segment = 100, overlap = 50))

    utils::setTxtProgressBar(pb, indx)

  }

  close(pb)

  #SINK #2
  log_file <- file(log_filepath, open = "a")
  sink(log_file, append=TRUE, split = TRUE, type='output')

  # Done comm
  cat(paste0('[', as.character(Sys.time()), '] ','Done!', '\n', sep=''))

  ### TAIL LENGTH LIST

  # header for progress bar
  cat(paste0('[', as.character(Sys.time()), '] ','Retrieving position calibrating data...', '\n', sep=''))


  #SINK #1
  sink(file=NULL, type = 'output')
  close(log_file)


  # progress bar
  pb <- utils::txtProgressBar(min = 0, max = length(index_list), style = 3, width = 50, char = "=")

  #create empty list for extracted data
  tail_length_list = list()

  # loop for parallel extraction
  for (indx in 1:length(index_list)){

    # work on subsets of reads in parallel
    tail_length_list <- c(tail_length_list, foreach::foreach(nam = names(tail_feature_list)[index_list[[indx]]])
                          %dopar% tail_feature_list[[nam]][[6]])

    utils::setTxtProgressBar(pb, indx)

  }

  close(pb)

  #SINK #2
  log_file <- file(log_filepath, open = "a")
  sink(log_file, append=TRUE, split = TRUE, type='output')

  #stop cluster
  doParallel::stopImplicitCluster()

  # coerce to df, add names
  overlap_count_list <- do.call("rbind.data.frame", overlap_count_list)
  overlap_count_list$readname <- names(tail_feature_list)
  colnames(overlap_count_list) <- c("total_chunk", "readname")

  tail_length_list <- do.call("rbind.data.frame", tail_length_list)
  tail_length_list$readname <- names(tail_feature_list)
  colnames(tail_length_list) <- c("tail_length", "readname")

  #produce final output
  coordinate_df <- merge(overlap_count_list, tail_length_list)

  # Done comm
  cat(paste0('[', as.character(Sys.time()), '] ','Done!', '\n', sep=''))

  #####################################################
  # ANALYZE RESULTS
  #####################################################

  #header for function
  cat(paste0('[', as.character(Sys.time()), '] ','Merging predictions...', '\n', sep=''))

  results <- analyze_results(nanopolish, coordinate_df, predicted_list, pass_only=pass_only)


  #dump output to files:
  names(results) <- c("binary_classified_reads", "detailed_positional_nonadenosine_residues")
  mapply(function (x,y) utils::write.table(x, file = file.path(save_dir, paste0(as.character(Sys.time()), '_', y, '.tsv')),
                                    row.names = F, sep="\t", quote = F), results, names(results))

  # Done comm
  cat(paste0('[', as.character(Sys.time()), '] ','Done!', '\n', sep=''))

  #display console messages
  cat(paste0('[', as.character(Sys.time()), '] ','Pipeline finished','\n'))
  cat(paste0('[', as.character(Sys.time()), '] ','The output files have been saved in: ','\n'))
  cat(paste0('  ', save_dir, '\n'))
  cat(paste0('[', as.character(Sys.time()), '] ','A logfile has been saved in: ','\n'))
  cat(paste0('  ', log_filepath, '\n'))

  cat(paste0('\n','Ninetails exited successfully.','\n','\n'))
  cat(paste0('Thank you for using Ninetails.'))



  return(results)

  # close logfile connection
  close(log_file)

}
