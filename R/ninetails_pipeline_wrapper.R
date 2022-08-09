#' Wrapper function for ninetails package.
#'
#' This function allows to perform all of the steps required to discover
#' nonadenosine nucleotides within the given dataset using ninetails.
#' Please keep in mind, that during computations the function creates large
#' segmentation data. Therefore it may be wise to split the nanopolish table
#' beforehand and then run this function on nanopolish table chunks.
#'
#' The output of this function is a list of 2 dataframes, containing: a) binary
#' classification of reads satisfying filtering criteria, b) detailed positional
#' info regarding all potential nonadenosine residues detected.
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
#' This parameter is set to 1 by default.
#'
#' @param basecall_group character string ["Basecall_1D_000"]. Name of the
#' level in the Fast5 file hierarchy from which the data should be extracted.
#' This parameter is set to 'Basecall_1D_000' by default.
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
#'results <- ninetails::check_tails(nanopolish = system.file('extdata', 'test_data', 'nanopolish_output.tsv',
#'                                                           package = 'ninetails'),
#'                                  sequencing_summary = system.file('extdata', 'test_data', 'sequencing_summary.txt',
#'                                                                   package = 'ninetails'),
#'                                  workspace = system.file('extdata', 'test_data', 'basecalled_fast5',
#'                                                          package = 'ninetails'),
#'                                  num_cores = 2,
#'                                  basecall_group = 'Basecall_1D_000',
#'                                  pass_only=TRUE,
#'                                  save_dir = '~/Downloads')
#'
#' }

check_tails <- function(nanopolish, sequencing_summary, workspace, num_cores=1, basecall_group="Basecall_1D_000", pass_only=TRUE, save_dir){

  # variable binding (suppressing R CMD check from throwing an error)
  i <- readname <- polya_length <- qc_tag <- chunkname <-  NULL

  # Create a log file
  if (dir.exists(file.path(save_dir))) {
    log_filename <- paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_ninetails.log", sep = "")
    log_filepath <- file.path(save_dir, log_filename, fsep = .Platform$file.sep)
    log_file <- file(log_filepath, open = "a")
    sink(log_file, append=TRUE, split = TRUE, type='output')
    on.exit(sink(file=NULL, type = 'output'))
  }

  #Show console message
  cat(paste0('Welcome to Ninetails ', as.character(utils::packageVersion("ninetails")), '\n', '\n', sep = ""))


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


  # Assertions
  if (missing(num_cores)) {
    sink(log_file, type='message')
    cat(paste0('[', as.character(Sys.time()), '] ','Ninetails encountered an error. Aborted.\n'))
    on.exit(sink())
    stop("Number of declared cores is missing. Please provide a valid num_cores argument.", call. =FALSE)
  }

  if (missing(basecall_group)) {
    sink(log_file, type='message')
    cat(paste0('[', as.character(Sys.time()), '] ','Ninetails encountered an error. Aborted.\n'))
    on.exit(sink())
    stop("Basecall group is missing. Please provide a valid basecall_group argument.", call. =FALSE)
  }

  if (missing(workspace)) {
    sink(log_file, type='message')
    cat(paste0('[', as.character(Sys.time()), '] ','Ninetails encountered an error. Aborted.\n'))
    on.exit(sink())
    stop("Directory with basecalled fast5s (guppy workspace) is missing. Please provide a valid workspace argument.", call. =FALSE)
  }

  if (missing(save_dir)) {
    sink(log_file, type='message')
    cat(paste0('[', as.character(Sys.time()), '] ','Ninetails encountered an error. Aborted.\n'))
    on.exit(sink())
    stop("A save dir for the output files is missing. Please provide a valid save_dir argument.", call. =FALSE)
  }

  sink(log_file, append=TRUE, split = F, type='message')
  assertthat::assert_that(assertive::is_numeric(num_cores),
                          msg=paste0("Declared core number must be numeric. Please provide a valid argument."))
  assertthat::assert_that(assertive::is_a_non_missing_nor_empty_string(workspace),
                          msg = paste0("Empty string provided as an input. Please provide a valid path to basecalled fast5 files."))
  assertthat::assert_that(assertive::is_character(workspace),
                          msg = paste0("Path to basecalled fast5 files is not a character string. Please provide a valid path to basecalled fast5 files."))
  assertthat::assert_that(assertive::is_character(save_dir),
                          msg = paste0("Path to output files is not a character string. Please provide a valid save_dir path."))
  sink(type='message')



  #####################################################
  # FAST5 QUICK CHECKUP
  #####################################################

  #checking data format
  ninetails::check_fast5_filetype(workspace, basecall_group)

  #####################################################
  # CREATING TAIL FEATURE LIST
  #####################################################

  # Extracting and processing polya & sequencing summary data
  polya_summary <- ninetails::extract_polya_data(nanopolish, sequencing_summary, pass_only)


  #create empty list for extracted fast5 data
  tail_features_list = list()

  # creating cluster for parallel computing
  my_cluster <- parallel::makeCluster(num_cores)
  doSNOW::registerDoSNOW(my_cluster)
  `%dopar%` <- foreach::`%dopar%`
  `%do%` <- foreach::`%do%`

  # header for progress bar
  cat(paste0('[', as.character(Sys.time()), '] ','Extracting features of provided reads...', '\n', sep=''))

  #SINK #1
  sink(file=NULL, type = 'output')
  close(log_file)

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
                                         .combine = c, .inorder = TRUE,
                                         .errorhandling = 'pass',
                                         .options.snow = opts) %dopar% {lapply(polya_summary$readname[i], function(x) ninetails::extract_tail_data(x,polya_summary,workspace,basecall_group))}

  close(pb)

  #SINK #2
  log_file <- file(log_filepath, open = "a")
  sink(log_file, append=TRUE, split = TRUE, type='output')

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

  #####################################################
  # CREATING TAIL CHUNK LIST
  #####################################################

  # header for progress bar
  cat(paste0('[', as.character(Sys.time()), '] ','Creating tail segmentation data...', '\n', sep=''))

  #SINK #1
  sink(file=NULL, type = 'output')
  close(log_file)

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
                                      .combine = c, .inorder = TRUE,
                                      .errorhandling = 'pass',
                                      .options.snow = opts) %dopar% {
                                        lapply(names(tail_feature_list[[1]][i]), function(x) ninetails::split_tail_centered(x,tail_feature_list))
                                      }

  close(pb)

  #SINK #2
  log_file <- file(log_filepath, open = "a")
  sink(log_file, append=TRUE, split = TRUE, type='output')

  #rename first level of nested list accordingly
  names(tail_chunk_list) <- names(tail_feature_list[[1]])

  # Done comm
  cat(paste0('[', as.character(Sys.time()), '] ','Done!', '\n', sep=''))

  #####################################################
  # CREATING GAF LIST
  #####################################################

  #create empty list for the data
  gaf_list = list()

  #progressbar header
  cat(paste0('[', as.character(Sys.time()), '] ','Computing angular summation fields...', '\n', sep=''))

  #SINK #1
  sink(file=NULL, type = 'output')
  close(log_file)

  # progress bar
  pb <- utils::txtProgressBar(min = 0,
                              max = length(tail_chunk_list),
                              style = 3,
                              width = 50,
                              char = "=",
                              file = stderr())
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)


  gaf_list <- foreach::foreach(i = seq_along(tail_chunk_list), .combine = c, .inorder = TRUE,
                               .errorhandling = 'pass',
                               .options.snow = opts) %dopar% {
                                 lapply(tail_chunk_list[[i]], function(x_ij) ninetails::combine_gafs(x_ij[['chunk_sequence']]))
                                 }

  close(pb)

  #SINK #2
  log_file <- file(log_filepath, open = "a")
  sink(log_file, append=TRUE, split = TRUE, type='output')

  # Done comm
  cat(paste0('[', as.character(Sys.time()), '] ','Done!', '\n', sep=''))

  #stop the cluster before tensorflow
  parallel::stopCluster(my_cluster)

  #####################################################
  # PREDICT CLASSES
  #####################################################

  predicted_list <- tryCatch({ninetails::predict_gaf_classes(gaf_list)}) #suppress tensorflow console messages

  #####################################################
  # CREATE OUTPUT
  #####################################################

  #read nanopolish data
  nanopolish_polya_table <- vroom::vroom(nanopolish,
                                         col_select=c(readname, polya_length, qc_tag),
                                         show_col_types = FALSE)

  # assertions #2
  assertthat::assert_that(assertive::is_a_non_missing_nor_empty_string(nanopolish),
                          msg = "Empty string provided as an input. Please provide a nanopolish as a string")
  assertthat::assert_that(assertive::is_existing_file(nanopolish),
                          msg=paste0("File ",nanopolish," does not exist",sep=""))
  assertthat::assert_that(assertive::has_rows(nanopolish_polya_table),
                          msg = "Empty data frame provided as an input (nanopolish). Please provide valid input")

  # HANDLE LENGTH/POSITION CALIBRATING DATA
  # reinstantiate the cluster
  my_cluster <- parallel::makeCluster(num_cores)
  on.exit(parallel::stopCluster(my_cluster))
  doSNOW::registerDoSNOW(my_cluster)
  `%dopar%` <- foreach::`%dopar%`
  `%do%` <- foreach::`%do%`

  #create empty list for the data
  tail_length_list = list()

  # header for progress bar
  cat(paste0('[', as.character(Sys.time()), '] ','Retrieving estimated length data...', '\n', sep=''))

  #SINK #1
  sink(file=NULL, type = 'output')
  close(log_file)

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
                                       .options.snow = opts) %dopar% {
                                         lapply(names(tail_feature_list[[1]][i]), function(x) length(tail_feature_list[[1]][[x]][[2]]))
                                       }

  close(pb)

  #SINK #2
  log_file <- file(log_filepath, open = "a")
  sink(log_file, append=TRUE, split = TRUE, type='output')

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

  #SINK #1
  sink(file=NULL, type = 'output')
  close(log_file)

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
                                          .options.snow = opts) %dopar% {
                                            lapply(tail_chunk_list[[i]], function(x) x[['chunk_start_pos']]+50)
                                          }
  close(pb)

  #SINK #2
  log_file <- file(log_filepath, open = "a")
  sink(log_file, append=TRUE, split = TRUE, type='output')

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
  moved_chunks_table$prediction[moved_chunks_table$prediction ==0] <- "A"
  moved_chunks_table$prediction[moved_chunks_table$prediction ==1] <- "C"
  moved_chunks_table$prediction[moved_chunks_table$prediction ==2] <- "G"
  moved_chunks_table$prediction[moved_chunks_table$prediction ==3] <- "U"

  #extract reads with move==1 and modification absent (not detected)
  moved_unmodified_readnames <- names(which(with(moved_chunks_table, tapply(prediction, readname, unique) == 'A')))

  # cleaned chunks_table
  moved_chunks_table <- subset(moved_chunks_table, !(readname %in% moved_unmodified_readnames))
  # delete A-containing rows
  moved_chunks_table <- moved_chunks_table[!(moved_chunks_table$prediction=="A"),]
  #merge data from feats & predictions
  moved_chunks_table <- dplyr::left_join(moved_chunks_table, non_a_position_list, by=c("readname", "chunkname"))

  #estimate non-A nucleotide position
  moved_chunks_table$est_nonA_pos <- round(moved_chunks_table$polya_length-((moved_chunks_table$polya_length*moved_chunks_table$centr_signal_pos)/moved_chunks_table$signal_length), digits=2)

  #clean up the output nonA table:
  moved_chunks_table <- moved_chunks_table[,c(3,2,8,6,7)]

  # Handle other (discarded) reads:
  discarded_reads <- nanopolish_polya_table[!nanopolish_polya_table$readname %in% moved_chunks_table$readname,]

  # Add filtering criterion: select only pass or pass $ suffclip
  if(pass_only == TRUE){
    discarded_reads <- discarded_reads %>%
      dplyr::filter(!readname %in% moved_chunks_table$readname) %>%
      dplyr::mutate(comments = dplyr::case_when(polya_length < 10 ~ "insufficient read length",
                                                qc_tag == "SUFFCLIP" ~ "not included in the analysis (pass only = T)",
                                                qc_tag == "ADAPTER" ~ "nanopolish qc failed",
                                                qc_tag == "NOREGION" ~ "nanopolish qc failed",
                                                qc_tag == "READ_FAILED_LOAD" ~ "nanopolish qc failed",
                                                readname %in% moved_unmodified_readnames ~ "move transition present, nonA residue undetected",
                                                TRUE ~ "move transition absent, nonA residue undetected"),
                    class = dplyr::case_when(polya_length < 10 ~ "unclassified",
                                             readname %in% moved_unmodified_readnames ~ "unmodified",
                                             comments == "move transition absent, nonA residue undetected" ~ "unmodified",
                                             TRUE ~ "unclassified"))
  } else {
    discarded_reads <- discarded_reads %>%
      dplyr::filter(!readname %in% moved_chunks_table$readname) %>%
      dplyr::mutate(comments = dplyr::case_when(polya_length < 10 ~ "insufficient read length",
                                                qc_tag == "ADAPTER" ~ "nanopolish qc failed",
                                                qc_tag == "NOREGION" ~ "nanopolish qc failed",
                                                qc_tag == "READ_FAILED_LOAD" ~ "nanopolish qc failed",
                                                readname %in% moved_unmodified_readnames ~ "move transition present, nonA residue undetected",
                                                TRUE ~ "move transition absent, nonA residue undetected"),
                    class = dplyr::case_when(polya_length < 10 ~ "unclassified",
                                             readname %in% moved_unmodified_readnames ~ "unmodified",
                                             comments == "move transition absent, nonA residue undetected" ~ "unmodified",
                                             TRUE ~ "unclassified"))
  }


  modified_reads <- nanopolish_polya_table[nanopolish_polya_table$readname %in% moved_chunks_table$readname,]
  modified_reads <- modified_reads %>% dplyr::mutate(class = "modified",
                                                     comments = "move transition present, nonA residue detected")

  #merge second tabular output:
  nanopolish_polya_table <- rbind(modified_reads, discarded_reads)

  #CREATE FINAL OUTPUT
  ninetails_output <- list()
  ninetails_output[['read_classes']] <- nanopolish_polya_table
  ninetails_output[['nonadenosine_residues']] <- moved_chunks_table


  #dump output to files:
  mapply(function (x,y) utils::write.table(x, file = file.path(save_dir, paste0(as.character(Sys.time()), '_', y, '.tsv')),
                                           row.names = F, sep="\t", quote = F), ninetails_output, names(ninetails_output))

  # Done comm
  cat(paste0('[', as.character(Sys.time()), '] ','Done!', '\n', '\n', '\n', sep=''))

  #display console messages
  cat(paste0('[', as.character(Sys.time()), '] ','Pipeline finished','\n'))
  cat(paste0('[', as.character(Sys.time()), '] ','The output files have been saved in: ','\n'))
  cat(paste0('  ', save_dir, '\n'))
  cat(paste0('[', as.character(Sys.time()), '] ','A logfile has been saved in: ','\n'))
  cat(paste0('  ', log_filepath, '\n'))

  cat(paste0('\n','Ninetails exited successfully.','\n','\n'))
  cat(paste0('Thank you for using Ninetails.'))

  # close logfile connection
  on.exit(close(log_file))


  return(ninetails_output)

}
