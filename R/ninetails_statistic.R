#' Performs Fisher's exact test on ninetails' merged output data.
#'
#' It runs Fisher's exact test for testing the null of independence of rows and
#' columns in a contingency table representing given transcript in ninetails
#' output data. This is a wrapper for fisher.test function from stats package
#' with additional features to facilitate data wrangling.
#'
#' It is suitable only for the pairwise comparisons (i.e. for 2x2 contingency
#' table), where 2 conditions (e.g. WT vs KO) are compared at once.
#'
#' The function allows the user to set a cutoff number of reads required for
#' the analysis.
#'
#' This function is intended to work under the calculate_fisher function.
#'
#' The function was inspired by the Nanotail package written & maintained by
#' Pawel Krawczyk (smaegol): https://github.com/LRB-IIMCB/nanotail/blob/dev/R/polya_stats.R
#'
#' Many thanks to the developer of original source code.
#'
#'
#' @param ninetails_data dataframe - the output of \code{\link{merge_nonA_tables}}
#' function (merged tabular output containing read classification &
#' non-A position data).
#'
#' @param grouping_factor [character string] the name of factor variable defining
#' groups/conditions (needs to have 2 levels!)
#'
#' @param base [character string] letter representing particular non-A nucleotide,
#' for which the statistics are meant to be computed. Currently function accepts
#' "C", "G", "U" and "all" arguments. The "C" value is set by default. The "all"
#' value is for the nonA residues alltogether.
#'
#' @param min_reads [numeric] minimum number of reads representing given
#' transcript to include it in the analysis
#'
#' @param transcript_id_column [character string] name of the column in which
#' the identifiers of the transcripts are stored. It is set to NA
#' by default.
#'
#' @return a tibble with results for given transcript, including pvalue,
#' adjusted pvalue, stats_code (the variable describing whether conditions
#' are met) and significance (FDR based on padj).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' test <- ninetails::nonA_fisher(ninetails_data=merged_nonA_tables,
#'                                grouping_factor = "sample_name",
#'                                base="C",
#'                                min_reads=100)
#'
#' }
#'
nonA_fisher <- function(ninetails_data,grouping_factor, base, min_reads=0, transcript_id_column=NA) {

  # var binding
  counts_unmod <- NULL

  # Assertions
  if (missing(ninetails_data)) {
    stop("Ninetails data are missing. Please provide a valid ninetails_data argument",
         call. = FALSE)
  }
  if (missing(base)) {
    stop("Base is missing. Please provide 'base' argument as character string (C, G or U).",
         call. =FALSE)
  }
  if (missing(transcript_id_column)) {
    stop("Transcript_id_column is missing. Please provide 'transcript_id_column' argument as character string.",
         call. =FALSE)
  }

  assertthat::assert_that(assertive::has_rows(ninetails_data),
                          msg = "Empty data.frame provided as an input")
  assertthat::assert_that(assertive::is_numeric(min_reads),
                          msg = "Non-numeric parameter provided (min_reads)")
  assertthat::assert_that(grouping_factor %in% colnames(ninetails_data),
                          msg=paste0(grouping_factor," is not a column of input dataset"))


  # if grouping factor has more than two levels
  if (length(levels(ninetails_data[[grouping_factor]]))>2) {
    if(is.na(condition1) && is.na(condition2)) {
      #throw error when no conditions for comparison are specified
      stop(paste0("grouping_factor ",grouping_factor," has more than 2 levels. Please specify condtion1 and condition2 to select comparison pairs"))
    } else {
      # filter input data leaving only specified conditions, dropping other factor levels
      assertthat::assert_that(condition1 %in% levels(ninetails_data[[grouping_factor]]),
                              msg=paste0(condition1," is not a level of ",grouping_factor," (grouping_factor)"))
      assertthat::assert_that(condition2 %in% levels(ninetails_data[[grouping_factor]]),
                              msg=paste0(condition2," is not a level of ",grouping_factor," (grouping_factor)"))
      assertthat::assert_that(condition2 != condition1,
                              msg="condition2 should be different than condition1")
      ninetails_data <- ninetails_data %>% dplyr::filter(!!rlang::sym(grouping_factor) %in% c(condition1,condition2)) %>% droplevels()
    }
  } else if (length(levels(ninetails_data[[grouping_factor]]))==1) {
    stop("Only 1 level present for grouping factor. Choose another groping factor for comparison")
  } else {
    condition1 = levels(ninetails_data[[grouping_factor]])[1]
    condition2 = levels(ninetails_data[[grouping_factor]])[2]
  }

  # initial status code
  stats_code = codes_stats = "OK"
  # calculate group counts
  group_counts = ninetails_data %>% dplyr::group_by(!!!rlang::syms(c(grouping_factor))) %>% dplyr::count()

  stats <- NA

  if (base=="C") {
    count_column <- "counts_C"
  } else if (base=="G") {
    count_column <- "counts_G"
  } else if (base=="U") {
    count_column <- "counts_U"
  } else if (base=="all") {
    #ninetails_data <- ninetails_data %>% dplyr::mutate(counts_nonA=)
    count_column <- "counts_nonA"
  } else {
    stop("Wrong non-A nucleotide defined. To compute statistics, please provide 'base' argument as character string (C, G, U or all).")
  }


  if (nrow(group_counts)==2) {
    if (group_counts[1,]$n < min_reads) {
      if (group_counts[2,]$n < min_reads) {
        stats_code = "B_LC"
      } else {
        stats_code = "G_LC"
      }
    } else if (group_counts[2,]$n < min_reads) {
      stats_code = "G_LC"
    } else {
      options(scipen = 999)

      # summarize nonAs
      contingency_table <- ninetails::summarize_nonA(merged_nonA_tables = ninetails_data,
                                                     summary_factors=grouping_factor,
                                                     transcript_id_column=transcript_id_column) %>%
        dplyr::select(!!rlang::sym(grouping_factor),
                      counts_unmod,
                      !!rlang::sym(count_column))
      contingency_table <- as.data.frame(contingency_table) # coerce tibble to df as setting names to tibble is deprecated
      row.names(contingency_table) <- contingency_table[[grouping_factor]] # set rownames
      contingency_table[[grouping_factor]] <- NULL # drop grouping col
      stats <- suppressWarnings(stats::fisher.test(contingency_table))$p.value


    }
  } else if (nrow(group_counts)==1) {
    stats_code = "G_NA"
  } else if (nrow(group_counts)==0) {
    stats_code = "B_NA"
  } else {
    stats_code = "ERR"
  }

  # create output
  stats <- tibble::tibble(p.value=stats,stats_code=as.character(stats_code))

  return(stats)

}

stat_codes_list = list(OK = "OK",
                       G1_NA = "GROUP1_NA",
                       G2_NA = "GROUP2_NA",
                       G1_LC = "G1_LOW_COUNT",
                       G2_LC = "G2_LOW_COUNT",
                       B_NA = "DATA FOR BOTH GROUPS NOT AVAILABLE",
                       B_LC = "LOW COUNTS FOR BOTH GROUPS",
                       G_LC = "LOW COUNT FOR ONE GROUP",
                       G_NA = "DATA FOR ONE GROUP NOT AVAILABLE",
                       ERR = "OTHER ERROR")



#' Performs Fisher's exact test for each transcript in ninetails
#' output data. Then, it performs the Benjamini-Hochberg procedure
#' (BH step-up procedure) to control the FDR.
#'
#' This is a wrapper for fisher.test function from stats package and p.adjust
#' functions with additional features to facilitate data wrangling.
#'
#' The function was inspired by the Nanotail package written & maintained by
#' Pawel Krawczyk (smaegol): https://github.com/LRB-IIMCB/nanotail/blob/dev/R/polya_stats.R
#'
#' Many thanks to the developer of original source code.
#'
#' @param ninetails_data dataframe - the output of ninetails::merge_nonA_tables
#' function (merged tabular output containing read classification &
#' non-A position data).
#'
#' @param transcript_id_column [character string] column with transcript id data
#' (default: "ensembl_transcript_id_short"; can be changed by the user)
#'
#' @param min_reads [numeric] minimum number of reads representing given
#' transcript to include it in the analysis. This parameter is set by default
#' to 0. Please keep in mind that taking into account many transcripts with low
#' coverage increases the risk of reject true null hypothesis
#' (Benjamini-Hochberg procedure).
#'
#' @param min_nonA_reads [numeric] minimum number of reads containing
#' nonadenosine residues (summary for C, G, U alltogether) per given
#' transcript to include it in the analysis. This parameter prevents from
#' considering too many observations as nonsignificant in the further pvalue
#' adjustation. In general, the non-A containing reads are a small fraction
#' of the total pool of reads. As a rule of thumb, additional filtering
#' criteria can provide more valuable information regarding the samples
#' (prevent from rejecting true null hypothesis).
#' This is set by default to 0.
#'
#' @param grouping_factor [character string] grouping variable
#' (e.g. "sample_name" - default)
#'
#' @param condition1 [character string] first level of `grouping_factor`
#' to use for comparison
#'
#' @param condition2 [character string] second level of `grouping_factor`
#' to use for comparison
#'
#' @param alpha [numeric] an alpha value to consider a hit significant.
#' Default: 0.05.
#'
#' @param base [character string] letter representing particular non-A nucleotide,
#' for which the statistics are meant to be computed. Currently function accepts
#' "C", "G", "U" and "all" arguments. The "C" value is set by default. The "all"
#' value is for the nonA residues alltogether.
#'
#' @param ... additional parameters to pass to nonA_fisher (under development)
#'
#' @return a summary table with pvalues, padj and significance levels
#' for each transcript  (tibble)
#'
#' @importFrom rlang :=
#'
#' @export
#'
#' @examples
#' \dontrun{
#' test <- ninetails::calculate_fisher(ninetails_data=merged_nonA_tables,
#'                                     transcript_id_column = "ensembl_transcript_id_short",
#'                                     min_reads=100,
#'                                     min_nonA_reads=10,
#'                                     grouping_factor = "sample_name",
#'                                     condition1="WT",
#'                                     condition2="KO",
#'                                     alpha=0.05,
#'                                     base="C")
#' }
calculate_fisher <- function(ninetails_data,
                             transcript_id_column = "ensembl_transcript_id_short",
                             min_reads = 0,
                             min_nonA_reads=0,
                             grouping_factor = "sample_name",
                             condition1=NA,
                             condition2=NA,
                             alpha=0.05,
                             base="C",
                             ...)
{

  # vr binding
  merged_nonA_tables<- sum_nonA <- counts_nonA <- contig <- data<- stats<- padj<- NULL

  # Assertions
  if (missing(ninetails_data)) {
    stop("Ninetails data are missing. Please provide a valid ninetails_data argument",
         call. = FALSE)
  }
  if (missing(transcript_id_column)) {
    stop("Transcript_id_column is missing. Please provide a valid transcript_id_column argument",
         call. = FALSE)
  }
  if (missing(min_reads)) {
    stop("Min_reads are missing. Please provide a valid min_reads argument",
         call. = FALSE)
  }
  if (missing(min_nonA_reads)) {
    stop("Min_nonA_reads data are missing. Please provide a valid min_nonA_reads argument",
         call. = FALSE)
  }
  if (missing(base)) {
    stop("Base definition is missing. Please provide a valid base argument",
         call. = FALSE)
  }

  assertthat::assert_that(assertive::is_numeric(min_reads),
                          msg=paste0("Min_reads must be numeric. Please provide a valid argument."))
  assertthat::assert_that(assertive::is_numeric(min_nonA_reads),
                          msg=paste0("Min_nonA_reads must be numeric. Please provide a valid argument."))
  assertthat::assert_that(assertive::is_numeric(alpha),
                          msg=paste0("Alpha must be numeric. Please provide a valid argument."))

  assertthat::assert_that(assertive::has_rows(ninetails_data),
                          msg = "Empty data.frame provided as an input")


  # if grouping factor has more than two levels
  if (length(levels(ninetails_data[[grouping_factor]]))>2) {
    if(is.na(condition1) && is.na(condition2)) {
      #throw error when no conditions for comparison are specified
      stop(paste0("grouping_factor ",grouping_factor," has more than 2 levels. Please specify condtion1 and condition2 to select comparison pairs"))
    }
    else {
      # filter input data leaving only specified conditions, dropping other factor levels
      assertthat::assert_that(condition1 %in% levels(ninetails_data[[grouping_factor]]),
                              msg=paste0(condition1," is not a level of ",grouping_factor," (grouping_factor)"))
      assertthat::assert_that(condition2 %in% levels(ninetails_data[[grouping_factor]]),
                              msg=paste0(condition2," is not a level of ",grouping_factor," (grouping_factor)"))
      assertthat::assert_that(condition2 != condition1,msg="condition2 should be different than condition1")

      ninetails_data <- ninetails_data %>% dplyr::filter(!!rlang::sym(grouping_factor) %in% c(condition1,condition2)) %>%
        dplyr::mutate() %>%
        droplevels()

    }
  }
  else if (length(levels(ninetails_data[[grouping_factor]]))==1) {
    stop("Only 1 level present for grouping factor. Choose another groping factor for comparison")
  }
  else {
    condition1 = levels(ninetails_data[[grouping_factor]])[1]
    condition2 = levels(ninetails_data[[grouping_factor]])[2]
  }

  # filter out transcripts with not enough amount of non-A reads among the whole pool of reads:
  mod_summarized <- ninetails_data %>% dplyr::ungroup() %>%
    dplyr::mutate(sum_nonA = rowSums(dplyr::across(dplyr::starts_with('prediction_')))) %>%
    dplyr::group_by(!!!rlang::syms(c(transcript_id_column,grouping_factor))) %>%
    dplyr::summarise(dplyr::across(c(sum_nonA), list(counts = ~ sum(.x != 0))), .groups= 'drop') %>%
    dplyr::rename_with(~stringr::str_replace(.x, '^\\w+_(\\w+)_(\\w+)', '\\2_\\1'), 3:dplyr::last_col())
  # apply filtering criterion (minimal nonA read content)
  mod_summarized_filtered <- mod_summarized %>% dplyr::filter(counts_nonA>=min_nonA_reads)
  #extract filtered trans
  contig <- as.name(transcript_id_column)
  mod_summarized_filtered <-  unique(mod_summarized_filtered[[contig]])

  ninetails_data <- ninetails_data %>% dplyr::filter(!!rlang::sym(contig) %in% mod_summarized_filtered)

  #ninetails_data <- ninetails_data[ninetails_data$contig %in% mod_summarized_filtered,]
  #ninetails_data <- ninetails_data[transcript_id_column %in% mod_summarized_filtered,]


  ninetails_data_stat <- ninetails_data %>%
    dplyr::mutate(transcript_id=get(c(transcript_id_column))) %>%
    dplyr::group_by(get(c(transcript_id_column))) %>%
    tidyr::nest()

  message("calculating statistics")
  ninetails_data_stat <- ninetails_data_stat %>%
    dplyr::mutate(stats=purrr::map(data,ninetails::nonA_fisher,grouping_factor=grouping_factor,min_reads=min_reads, base=base,transcript_id_column = transcript_id_column)) %>%
    dplyr::select(-data) %>% tidyr::unnest(cols = c(stats)) %>% dplyr::rename(!!transcript_id_column := "get(c(transcript_id_column))")


  message("Finished")

  message("Adjusting p.value")
  ninetails_data_stat$padj <- stats::p.adjust(ninetails_data_stat$p.value, method = "BH")

  # create significance factor
  ninetails_data_stat <- ninetails_data_stat %>%
    dplyr::mutate(significance = dplyr::case_when(is.na(padj)  ~ "NotSig",
                                                  (padj < alpha) ~ paste0("FDR<", alpha),
                                                  TRUE ~ "NotSig"))
  ninetails_data_stat$stats_code <- sapply(ninetails_data_stat$stats_code,
                                           FUN = function(x) {stat_codes_list[[x]]},simplify = "vector",USE.NAMES = FALSE) %>% unlist()



  ninetails_data_stat <- ninetails_data_stat %>% dplyr::arrange(padj)

  return(ninetails_data_stat)
}
