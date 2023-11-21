#' Wrapper function for ninetails package.
#'
#' This function allows to perform all of the steps required to discover
#' nonadenosine nucleotides within the given dataset using ninetails.
#' Please keep in mind, that during computations the function creates large
#' segmentation data. Therefore it may be wise to split the nanopolish table
#' beforehand and then run this function on nanopolish table chunks.
#'
#' The output of this function is a list of 2 dataframes, containing:\itemize{
#' \item read_classes - classification of reads based on applied criteria
#' \item nonadenosine_residues - detailed positional info regarding all
#' potential nonadenosine residues detected.
#' }
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
#' @param save_dir character string. Full path of the directory where the output
#' files containing the tail composition information should be stored.
#'
#' @export
#'
#' @return A list containing tail information organized by the read ID
#' is returned. Always assign this returned list to a variable, otherwise
#' the long list will be printed to the console, which may crash your R session.
#' Also a log file as well as read_classes & nonadenosine_residues files
#' are created in the user-specified directory.
#'
#' @importFrom foreach %dopar%
#' @importFrom utils head
#'
#' @examples
#'\dontrun{
#'
#'results <- ninetails::check_tails(
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
#'  pass_only=TRUE,
#'  qc=TRUE,
#'  save_dir = '~/Downloads')
#'
#' }
check_tails <- function(nanopolish,
                        sequencing_summary,
                        workspace,
                        num_cores=1,
                        basecall_group="Basecall_1D_000",
                        pass_only=TRUE,
                        qc=TRUE,
                        save_dir){

  # variable binding (suppressing R CMD check from throwing an error)
  i <- readname <- polya_length <- qc_tag <- chunkname <-contig <- est_nonA_pos <-  NULL

  tryCatch({
    #show console message
    cat(paste0('\n', 'Welcome to Ninetails ', as.character(utils::packageVersion("ninetails")), '\n', '\n', sep = ""))

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


    #Assertions
    sink(log_file, append=TRUE, split = F, type='message')
    assertthat::assert_that(is.numeric(num_cores),
                            msg=paste0("Declared core number must be numeric. Please provide a valid argument.",
                                       '\n','[', as.character(Sys.time()), '] ',
                                       "Ninetails aborted"))
    assertthat::assert_that(is.character(workspace),
                            msg = paste0("Path to basecalled fast5 files is not a character string. Please provide a valid path to basecalled fast5 files.",
                                         '\n','[', as.character(Sys.time()), '] ',
                                         "Ninetails aborted"))
    assertthat::assert_that(checkmate::test_string(workspace, null.ok=F, min.chars=1),
                            msg = paste0("Empty string provided as an input. Please provide a valid path to basecalled fast5 files.",
                                         '\n','[', as.character(Sys.time()), '] ',
                                         "Ninetails aborted"))
    assertthat::assert_that(is.character(save_dir),
                            msg = paste0("Path to output files is not a character string. Please provide a valid save_dir path.",
                                         '\n','[', as.character(Sys.time()), '] ',
                                         "Ninetails aborted"))
    assertthat::assert_that(is.logical(pass_only),
                            msg = paste0("The pass_only is not a logical [TRUE/FALSE]. Please provide a valid argument.",
                                         '\n','[', as.character(Sys.time()), '] ',
                                         "Ninetails aborted"))
    assertthat::assert_that(is.logical(qc),
                            msg = paste0("The qc variable is not a logical [TRUE/FALSE]. Please provide a valid argument.",
                                         '\n','[', as.character(Sys.time()), '] ',
                                         "Ninetails aborted"))
    sink(type='message')


    # user-specified options
    cat(paste0('[', as.character(Sys.time()), '] ',' Ninetails was launched with following options:', '\n', sep=''))

    if (!is.object(nanopolish)){
      cat(paste0(' nanopolish output file:       ', nanopolish, '\n', sep=''))
    } else {
      cat(paste0(' nanopolish output file:       ', deparse(substitute(nanopolish)), '\n', sep=''))
    }
    if (!is.object(sequencing_summary)){
      cat(paste0(' sequencing sumary file:       ', sequencing_summary, '\n', sep=''))
    } else {
      cat(paste0(' sequencing sumary file:       ', deparse(substitute(sequencing_summary)), '\n', sep=''))
    }
    cat(paste0(' fast5 files directory:        ', workspace, '\n', sep=''))
    cat(paste0(' number of cores:              ', num_cores, '\n', sep=''))
    cat(paste0(' basecall group:               ', basecall_group, '\n', sep=''))
    cat(paste0(' only "PASS" reads included:   ', pass_only, '\n', sep=''))
    cat(paste0(' output quality control:       ', qc, '\n', sep=''))
    cat(paste0(' output directory:             ', save_dir, '\n', '\n', '\n'))

    cat(paste0('[', as.character(Sys.time()), '] ', 'Pipeline initialized','\n','\n'))


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
    on.exit(parallel::stopCluster(my_cluster))
    doSNOW::registerDoSNOW(my_cluster)
    `%dopar%` <- foreach::`%dopar%`
    `%do%` <- foreach::`%do%`
    mc_options <- list(preschedule = TRUE, set.seed = FALSE, cleanup = TRUE)

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
                                           .options.snow = opts,
                                           .options.multicore = mc_options) %dopar% {lapply(polya_summary$readname[i], function(x) ninetails::extract_tail_data(x,polya_summary,workspace,basecall_group))}

    close(pb)

    #SINK #2
    log_file <- file(log_filepath, open = "a")
    sink(log_file, append=TRUE, split = TRUE, type='output')

    #label each signal according to corresponding read name to avoid confusion
    #squiggle_names <- polya_summary$readname
    #names(tail_features_list) <- squiggle_names
    squiggle_names <- as.vector(sapply(tail_features_list, function(x) attributes(x[[1]])$names))
    tail_features_list <- stats::setNames(tail_features_list, squiggle_names)


    ####### SANITY CHECK #########################################################
    # slows down pipeline a little bit, but facilitates debugging in case of error

    # check whether a list of lists was created - remove unwanted items
    tail_features_list <- Filter(function(x) is.list(x), tail_features_list)
    # check the sublists' structure - remove unwanted items
    tail_features_list <- Filter(function(x) is.numeric(x$tail_signal)& is.numeric(x$tail_moves), tail_features_list)

    if (length(tail_features_list) == 0){
      stop("Produced feature list is of length 0. Check your input (nanopolish, fast5 files) for integrity.", call. =FALSE)
    }

    warn_message <- FALSE

    if (length(squiggle_names) != length(unique(names(tail_features_list)))) {
      warn_message <- TRUE
    }

    ####### SANITY CHECK ########################################################
    # remove reads with only zero moved tails
    tail_features_list <- Filter(function(x) sum(x$tail_moves) !=0, tail_features_list)
    zeromoved_readnames <- squiggle_names[!(squiggle_names %in% names(tail_features_list))]

    # prevent from running on reads which do not fulfill the pseudomove condition
    tail_features_list <- Filter(function(x)any(with(rle(x$tail_pseudomoves), lengths[values!=0]>=5)), tail_features_list)

    if (length(tail_features_list) == 0){
      sink(log_file, type='message')
      cat(paste0('[', as.character(Sys.time()), '] ','Ninetails encountered an error. Aborted.\n'))
      on.exit(sink())
      stop("None of the provided signals fulfilled filtering criteria. Check your fast5 files. Re-basecall could be a good idea.", call. =FALSE)
    }

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
                                        .combine = c,
                                        .inorder = TRUE,
                                        .errorhandling = 'pass',
                                        .options.snow = opts,
                                        .options.multicore = mc_options) %dopar% {
                                          lapply(names(tail_feature_list[[1]][i]), function(x) ninetails::split_tail_centered(x,tail_feature_list))
                                        }

    close(pb)

    #SINK #2
    log_file <- file(log_filepath, open = "a")
    sink(log_file, append=TRUE, split = TRUE, type='output')

    #rename first level of nested list accordingly
    names(tail_chunk_list) <- names(tail_feature_list[[1]])

    #drop moves variable - recursive fn for prunning moves vector from the output
    .prune_moves <- function(i)
      lapply(i, function(x)
        if (is.list(x)) {
          if(!is.null(names(x))) .prune_moves(x[names(x)!="chunk_moves"]) else .prune_moves(x)
        } else x
      )

    tail_chunk_list <- .prune_moves(tail_chunk_list)

    # Done comm
    cat(paste0('[', as.character(Sys.time()), '] ','Done!', '\n', sep=''))

    #####################################################
    # CREATING GAF LIST
    #####################################################

    #create empty list for the data
    gaf_list = list()

    #progressbar header
    cat(paste0('[', as.character(Sys.time()), '] ','Computing gramian angular fields...', '\n', sep=''))

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


    gaf_list <- foreach::foreach(i = seq_along(tail_chunk_list),
                                 .combine = c,
                                 .inorder = TRUE,
                                 .errorhandling = 'pass',
                                 .options.snow = opts,
                                 .options.multicore = mc_options) %dopar% {
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
    # Accept either path to file or in-memory file - PK's GH issue
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
    # reinstantiate the cluster
    my_cluster <- parallel::makeCluster(num_cores)
    on.exit(parallel::stopCluster(my_cluster))
    doSNOW::registerDoSNOW(my_cluster)
    `%dopar%` <- foreach::`%dopar%`
    `%do%` <- foreach::`%do%`
    mc_options <- list(preschedule = TRUE, set.seed = FALSE, cleanup = TRUE)

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
                                         .combine = c,
                                         .inorder = TRUE,
                                         .errorhandling = 'pass',
                                         .options.snow = opts,
                                         .options.multicore = mc_options) %dopar% {
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
                                            .combine = c,
                                            .inorder = TRUE,
                                            .errorhandling = 'pass',
                                            .options.snow = opts,
                                            .options.multicore = mc_options) %dopar% {
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
    decorated_reads <- decorated_reads %>% dplyr::mutate(class = "decorated",
                                                         comments = "YAY")

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
      decorated_reads_edited <- unique(decorated_reads_edited)
      moved_chunks_table_qc <- unique(moved_chunks_table_qc)

      ninetails_output[['read_classes']] <- decorated_reads_edited
      ninetails_output[['nonadenosine_residues']] <- moved_chunks_table_qc

    } else{
      #CREATE FINAL OUTPUT
      nanopolish_polya_table <- unique(nanopolish_polya_table)
      moved_chunks_table <- unique(moved_chunks_table)

      ninetails_output[['read_classes']] <- nanopolish_polya_table
      ninetails_output[['nonadenosine_residues']] <- moved_chunks_table
    }

    #dump output to files:
    mapply(function (x,y) utils::write.table(x, file = file.path(save_dir, paste0(as.character(format(Sys.time(), "%Y-%m-%d_%H-%M-%S")), '_', y, '.tsv')),
                                             row.names = F, sep="\t", quote = F), ninetails_output, names(ninetails_output))

    # Done comm
    cat(paste0('[', as.character(Sys.time()), '] ','Done!', '\n', '\n', '\n', sep=''))

    #display console messages
    cat(paste0('[', as.character(Sys.time()), '] ','Pipeline finished','\n'))
    cat(paste0('[', as.character(Sys.time()), '] ','The output files have been saved in: ','\n'))
    cat(paste0('  ', save_dir, '\n'))
    cat(paste0('[', as.character(Sys.time()), '] ','A logfile has been saved in: ','\n'))
    cat(paste0('  ', log_filepath, '\n'))


    # WARN OR FULL SUCCESS #######################################################

    if (warn_message == TRUE) {
      cat(paste0('\n','Ninetails exited with WARNING.','\n','\n','Check your Nanopolish and Guppy outputs for consistency. It seems like some reads passing quality criteria, which were provided in nanopolish output file (*.tsv) are absent from provided fast5 files. They were omitted in this analysis. However they might constitute a significant fraction of your data. If you are sure that everything was provided correctly, just ignore that info. ','\n','\n'))

    } else {

      cat(paste0('\n','Ninetails exited successfully.','\n','\n'))

    }



    cat(paste0('Thank you for using Ninetails.'))

    on.exit(closeAllConnections()) #fixed err

    return(ninetails_output)

  },error=function(e){ cat(paste0(e$message,'\n','[', as.character(Sys.time()), '] ', 'Ninetails aborted.\n', sep=''))})

}
