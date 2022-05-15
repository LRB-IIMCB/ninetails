#' Creates a dataframe with segmentation and tail length info required for
#' non-adenosine position estimation.
#'
#' @param tail_feature_list list object produced by create_tail_feature_list function.
#'
#' @param num_cores numeric [1]. Number of physical cores to use in processing
#' the data. Do not exceed 1 less than the number of cores at your disposal.
#'
#' @return a dataframe with 3 columns: read ID (readname), total number of segments
#' per given tail (total_chunk) and poly(A) tail length (tail_length) estimated by
#' nanopolish.
#'
#' @importFrom foreach %dopar%
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' create_coordinate_dataframe(tail_feature_list = '/path/to/tail_feature_list',
#'                             num_cores = 10)
#'
#'}
#'
#'
create_coordinate_dataframe <- function(tail_feature_list, num_cores){

  # variable binding (suppressing R CMD check from throwing an error)
  nam <- NULL

  #assertions
  if (missing(num_cores)) {
    stop("Number of declared cores is missing. Please provide a valid num_cores argument.", call. =FALSE)
  }

  if (missing(tail_feature_list)) {
    stop("List of features is missing. Please provide a valid tail_feature_list argument.", call. =FALSE)
  }

  assertthat::assert_that(assertive::is_numeric(num_cores), msg = paste0("Declared core number must be numeric. Please provide a valid argument."))
  assertthat::assert_that(assertive::is_list(tail_feature_list),
                          msg = paste0("Given tail_feature_list is not a list (class). Please provide valid file format."))

  # creating cluster for parallel computing
  doParallel::registerDoParallel(cores = num_cores)

  # this is list of indexes required for parallel computing; the main list of reads is split for chunks
  index_list = split(1:length(names(tail_feature_list)), ceiling(1:length(names(tail_feature_list))/100))

  # OVERLAP COUNT LIST
  # header for progress bar
  cat(paste0('[', as.character(Sys.time()), '] ','Retrieving segmentation data...', '\n', sep=''))

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

  # Done comm
  cat(paste0('[', as.character(Sys.time()), '] ','Done!', '\n', sep=''))

  ### TAIL LENGTH LIST
  # header for progress bar
  cat(paste0('[', as.character(Sys.time()), '] ','Retrieving position calibrating data...', '\n', sep=''))

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


  return(coordinate_df)

}



#' Summarizes classification results into list of 2 dataframes containing detailed
#' (per each position) and preliminary output (binary classification of reads for
#' modified - containing non-adenosine residue(s) and unmodified - with adenosines
#' exclusively).
#'
#' @param coordinate_df a data frame object produced by create_coordinate_df
#' function.
#' @param predicted_list a list object produced by predict_classes function.
#'
#' @return a list of 2 dataframes with prediction results. First dataframe contains
#' detailed info on type and estimated positions of non-adenosine residues detected.
#' The second dataframe contains initial indications, whether the given read was
#' recognized as modified (containing non-adenosine residue) or not.
#'
#' @importFrom foreach %dopar%
#' @importFrom dplyr %>%
#'
#' @export
#'
#' @examples
#'\dontrun{
#'
#' analyze_results(coordinate_df = '/path/to/coordinate_df',
#'                 predicted_list = '/path/to/predicted_list')
#'
#'}
#'
analyze_results <- function(coordinate_df, predicted_list){

  # variable binding (suppressing R CMD check from throwing an error)
  chunk <- total_chunk <- tail_length <- init <- NULL

  #assertions
  if (missing(coordinate_df)) {
    stop("Data frame of coordinate calibration values is missing. Please provide a valid coordinate_df argument.", call. =FALSE)
  }

  if (missing(predicted_list)) {
    stop("List of predictions is missing. Please provide a valid predicted_list argument.", call. =FALSE)
  }



  assertthat::assert_that(assertive::is_data.frame(coordinate_df),
                          msg = paste0("Given coordinate_df is not a data frame (class). Please provide valid file format."))
  assertthat::assert_that(assertive::is_list(predicted_list),
                          msg = paste0("Given predicted_list is not a list (class). Please provide valid file format."))


  #FIRST LIST
  moved_chunks_table <- data.frame(t(Reduce(rbind, predicted_list)))

  moved_chunks_table <- tidyr::separate(moved_chunks_table, init, into = c("readname", "chunk"), sep = "\\_")
  moved_chunks_table$chunk <- as.numeric(moved_chunks_table$chunk)
  colnames(moved_chunks_table)[3] <- "prediction"

  moved_chunks_table$prediction[moved_chunks_table$prediction ==0] <- "A"
  moved_chunks_table$prediction[moved_chunks_table$prediction ==1] <- "C"
  moved_chunks_table$prediction[moved_chunks_table$prediction ==2] <- "G"
  moved_chunks_table$prediction[moved_chunks_table$prediction ==3] <- "U"

  #extract reads with move==1 and modification absent (not detected)
  moved_unmodified_readnames <- names(which(with(moved_chunks_table, tapply(prediction, readname, unique) == 'A')))

  # cleaned chunks_table
  moved_chunks_table <- subset(moved_chunks_table, !(readname %in% moved_unmodified_readnames))
  moved_chunks_table <- dplyr::left_join(moved_chunks_table, coordinate_df,by="readname")

  # delete A-containing rows
  moved_chunks_table<-moved_chunks_table[!(moved_chunks_table$prediction=="A"),]

  # total chunks table
  readname <- coordinate_df$readname

  ## bin chunks according to the position:
  moved_chunks_table <- moved_chunks_table %>%
    dplyr::mutate(interval = dplyr::case_when(total_chunk==1 ~ "3'end",
                                              total_chunk==2 ~ dplyr::case_when(chunk==1 ~ "3'end",
                                                                                TRUE~ "5'end"),
                                              total_chunk==3 ~ dplyr::case_when(chunk==1 ~ "3'end",
                                                                                chunk==2 ~ "center",
                                                                                TRUE~ "5'end"),
                                              TRUE ~ dplyr::case_when(chunk==1 ~ "3'end",
                                                                      chunk < (total_chunk/3) ~ "3'distal",
                                                                      chunk < 2*(total_chunk/3) ~ "center",
                                                                      chunk==total_chunk ~ "5'end",
                                                                      TRUE ~ "5'distal")))
  #estimated hit centered_position
  moved_chunks_table <- moved_chunks_table %>% dplyr::mutate(centered_position = round(((tail_length/total_chunk)*chunk),2))

  ## SECOND LIST
  # handle unmodified IDs
  unmodified_readnames <- readname[!(readname %in% moved_chunks_table$readname)]

  preliminary_classification  <- data.frame(readname)
  preliminary_classification <- preliminary_classification %>%
    dplyr::mutate(class = ifelse((readname %in% unmodified_readnames), "unmodified", "modified"))

  # factorize variables:
  moved_chunks_table$interval <- factor(moved_chunks_table$interval, levels= c("3'end", "3'distal", "center", "5'distal", "5'end"))
  moved_chunks_table$prediction <- factor(moved_chunks_table$prediction, levels= c("C", "G", "U")) # add "A" if drop A would be commented out

  preliminary_classification$class <- factor(preliminary_classification$class, levels=c("unmodified", "modified"))

  #create final output
  analyzed_results_list <- list(moved_chunks_table, preliminary_classification)

  return(analyzed_results_list)
}


#' Classifies a dataframe containing classification of discarded reads.
#'
#' The output dataframe includes reads which were not taken into consideration
#' in ninetails classification pipeline (i.e. reads, which did not fulfill
#' the analysis requirements).
#'
#' The output dataframe consists of 2 columns: readname and class. Class column
#' contains short explanation of main reason why the given read was omitted.
#' The reason behind this could be (1) insufficient read quality according to the
#' nanopolish polya function, (2) insufficient polyA tail length, which makes
#' segmentation and classification impossible, (3) filtering criterion applied
#' by the user (reads marked by nanopolish as "SUFFCLIP" may be included/excluded,
#' depending on user's preference).
#'
#' @param analyzed_results_list a list object produced by analyze_results function.
#' @param nanopolish character string. Full path of the .tsv file produced
#' by nanopolish polya function.
#'
#' @return a dataframe with characteristics of reads excluded from ninetails
#' analysis. The object contains readnames and classification - column containing
#' the reason why the particular read was discarded.
#'
#' @importFrom dplyr %>%
#'
#' @export
#'
#' @examples
#'\dontrun{
#'
#' handle_discarded_reads(analyzed_results_list = list_object,
#'                        nanopolish = "path/to/nanopolish.tsv/file")
#'
#'}
#'
handle_discarded_reads<- function(analyzed_results_list, nanopolish){

  # variable binding (suppressing R CMD check from throwing an error)
  readname <- polya_length <- qc_tag  <- predicted_list <- NULL

  #assertions
  if (missing(nanopolish)) {
    stop("Nanopolish polya output is missing. Please provide a valid nanopolish argument.", .call = FALSE)
  }

  if (missing(analyzed_results_list)) {
    stop("Analyzed_results_list object is missing. Please provide a valid analyzed_results_list argument.", .call = FALSE)
  }

  assertthat::assert_that(assertive::is_list(analyzed_results_list),
                          msg = paste0("Given analyzed_results_list is not a list (class). Please provide valid object."))

  #read nanopolish polya
  nanopolish_polya_table <- vroom::vroom(nanopolish, col_select=c(readname, polya_length, qc_tag), show_col_types = FALSE)


  #assertions
  assertthat::assert_that(assertive::is_a_non_missing_nor_empty_string(nanopolish),msg = "Empty string provided as an input. Please provide a nanopolish as a string")
  assertthat::assert_that(assertive::is_existing_file(nanopolish), msg=paste0("File ",nanopolish," does not exist",sep=""))
  assertthat::assert_that(assertive::has_rows(nanopolish_polya_table), msg = "Empty data frame provided as an input (nanopolish). Please provide valid input")


  preliminary_classification <- predicted_list[[2]]

  #reads with sufficient tail length and move==1
  classified_reads <- preliminary_classification$readname

  discarded_reads <- nanopolish_polya_table %>%
    dplyr::filter(!readname %in% classified_reads) %>%
    dplyr::mutate(class = dplyr::case_when(polya_length < 10 ~ "unclassified - insufficient read length",
                                           qc_tag == "SUFFCLIP" ~ "unclassified - nanopolish qc failed (suffclip)",
                                           qc_tag == "ADAPTER" ~ "unclassified - nanopolish qc failed (adapter)",
                                           qc_tag == "NOREGION" ~ "unclassified - nanopolish qc failed (noregion)",
                                           qc_tag == "READ_FAILED_LOAD" ~ "unclassified - nanopolish qc failed (read_failed_load)",
                                           TRUE ~ "unclassified - required move not found")) %>%
    dplyr::select(-c(qc_tag, polya_length))

  return(discarded_reads)
}
