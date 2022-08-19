#' Creates the list object containing tabular outputs of ninetails pipeline.
#'
#' @param tail_feature_list list object produced by create_tail_feature_list
#' function.
#'
#' @param tail_chunk_list list object produced by create_tail_chunk_list
#' function.
#'
#' @param nanopolish character string. Full path of the .tsv file produced
#' by nanopolish polya function.
#'
#' @param predicted_list a list object produced by predict_classes function.
#'
#' @param num_cores numeric [1]. Number of physical cores to use in processing
#' the data. Do not exceed 1 less than the number of cores at your disposal.
#'
#' @param pass_only logical [TRUE/FALSE]. If TRUE, only reads tagged by
#' nanopolish as "PASS" would be taken into consideration. Otherwise, reads
#' tagged as "PASS" & "SUFFCLIP" will be taken into account in analysis.
#' As a default, "TRUE" value is set.
#'
#' @return This function returns a list object containing two fataframes:
#' "read_classes" and "nonadenosine_residues" with the final output.
#' First dataframe contains initial indications, whether the given read was
#' classified or omitted (with reason) and if classified, whether read was
#' recognized as modified (containing non-adenosine residue) or not.
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
#'                pass_only=TRUE)
#'}
#'
#'
create_outputs <- function(tail_feature_list,
                           tail_chunk_list,
                           nanopolish,
                           predicted_list,
                           num_cores,
                           pass_only=TRUE){

  #variable binding
  readname <- polya_length <- qc_tag<- i <- chunkname <- NULL

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

  assertthat::assert_that(assertive::is_list(tail_feature_list),
                          msg = paste0("Given tail_feature_list is not a list (class). Please provide valid object."))
  assertthat::assert_that(assertive::is_list(tail_chunk_list),
                          msg = paste0("Given tail_chunk_list is not a list (class). Please provide valid object."))
  assertthat::assert_that(assertive::is_list(predicted_list),
                          msg = paste0("Given predicted_list is not a list (class). Please provide valid object."))


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
  # creating cluster for parallel computing
  my_cluster <- parallel::makeCluster(num_cores)
  on.exit(parallel::stopCluster(my_cluster))

  doSNOW::registerDoSNOW(my_cluster)
  `%dopar%` <- foreach::`%dopar%`
  `%do%` <- foreach::`%do%`

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
                                       .options.snow = opts) %dopar% {
                                         lapply(names(tail_feature_list[[1]][i]), function(x) length(tail_feature_list[[1]][[x]][[2]]))
                                       }

  close(pb)

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
                                          .options.snow = opts) %dopar% {
                                            lapply(tail_chunk_list[[i]], function(x) x[['chunk_start_pos']]+50)
                                          }

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

  return(ninetails_output)

}
