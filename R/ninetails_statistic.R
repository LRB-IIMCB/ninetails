#' Perform Fisher's exact test on a single transcript in ninetails output
#'
#' Runs Fisher's exact test for testing the null hypothesis of independence
#' of rows and columns in a 2x2 contingency table representing a given
#' transcript in ninetails output data. This is a wrapper for
#' \code{\link[stats]{fisher.test}} with additional data wrangling features.
#'
#' The function is suitable only for pairwise comparisons (2x2 contingency
#' tables), where two conditions (e.g. WT vs KO) are compared at once.
#' The user can set a cutoff number of reads required for the analysis.
#'
#' This function is intended to be called internally by
#' \code{\link{calculate_fisher}} and is not typically used standalone.
#'
#' @section Acknowledgements:
#' Inspired by the Nanotail package written and maintained by
#' Pawel Krawczyk (smaegol):
#' \url{https://github.com/LRB-IIMCB/nanotail/blob/dev/R/polya_stats.R}.
#' Many thanks to the developer of the original source code.
#'
#' @param ninetails_data Data frame. The output of
#'   \code{\link{merge_nonA_tables}} (merged tabular output containing read
#'   classification and non-A position data).
#'
#' @param grouping_factor Character string. The name of the factor variable
#'   defining groups/conditions (must have exactly 2 levels).
#'
#' @param base Character string. Letter representing the non-A nucleotide
#'   for which statistics are computed. Accepted values: \code{"C"},
#'   \code{"G"}, \code{"U"}, or \code{"all"} (for all non-A residues
#'   combined). Default: \code{"C"}.
#'
#' @param min_reads Numeric. Minimum number of reads representing a given
#'   transcript to include it in the analysis. Default: 0.
#'
#' @param transcript_id_column Character string. Name of the column in
#'   which transcript identifiers are stored. Default: \code{NA}.
#'
#' @return A tibble with results for the given transcript, containing:
#'   \describe{
#'     \item{p.value}{Numeric. Raw p-value from Fisher's exact test, or
#'       \code{NA} if conditions were not met.}
#'     \item{stats_code}{Character. Status code describing whether the test
#'       conditions were met (see \code{stat_codes_list} for definitions).}
#'   }
#'
#' @seealso
#' \code{\link{calculate_fisher}} for the per-transcript wrapper with
#' multiple-testing correction,
#' \code{\link{merge_nonA_tables}} for preparing the input,
#' \code{\link{summarize_nonA}} for contingency table computation
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' test <- ninetails::nonA_fisher(
#'   ninetails_data = merged_nonA_tables,
#'   grouping_factor = "sample_name",
#'   base = "C",
#'   min_reads = 100
#' )
#'
#' }
#'
nonA_fisher <- function(ninetails_data,
                        grouping_factor,
                        base,
                        min_reads = 0,
                        transcript_id_column = NA) {

  # Assertions
  if (missing(ninetails_data)) {
    stop(
      "Ninetails data are missing. Please provide a valid ninetails_data argument",
      call. = FALSE
    )
  }
  if (missing(base)) {
    stop(
      "Base is missing. Please provide 'base' argument as character string (C, G or U).",
      call. = FALSE
    )
  }
  if (missing(transcript_id_column)) {
    stop(
      "Transcript_id_column is missing. Please provide 'transcript_id_column' argument as character string.",
      call. = FALSE
    )
  }

  if (!is.data.frame(ninetails_data) || nrow(ninetails_data) == 0) {
    stop(
      "Empty data frame provided as an input (ninetails_data). Please provide valid input"
    )
  }

  assert_condition(
    is.numeric(min_reads),
    "Non-numeric parameter provided (min_reads)"
  )
  assert_condition(
    grouping_factor %in% colnames(ninetails_data),
    paste0(grouping_factor, " is not a column of input dataset")
  )

  # if grouping factor has more than two levels
  if (length(levels(ninetails_data[[grouping_factor]])) > 2) {
    if (is.na(condition1) && is.na(condition2)) {
      #throw error when no conditions for comparison are specified
      stop(paste0(
        "grouping_factor ",
        grouping_factor,
        " has more than 2 levels. Please specify condtion1 and condition2 to select comparison pairs"
      ))
    } else {
      # filter input data leaving only specified conditions, dropping other factor levels
      assert_condition(
        condition1 %in% levels(ninetails_data[[grouping_factor]]),
        paste0(
          condition1,
          " is not a level of ",
          grouping_factor,
          " (grouping_factor)"
        )
      )
      assert_condition(
        condition2 %in% levels(ninetails_data[[grouping_factor]]),
        paste0(
          condition2,
          " is not a level of ",
          grouping_factor,
          " (grouping_factor)"
        )
      )
      assert_condition(
        condition2 != condition1,
        "condition2 should be different than condition1"
      )
      ninetails_data <- ninetails_data %>%
        dplyr::filter(
          !!rlang::sym(grouping_factor) %in% c(condition1, condition2)
        ) %>%
        droplevels()
    }
  } else if (length(levels(ninetails_data[[grouping_factor]])) == 1) {
    stop(
      "Only 1 level present for grouping factor. Choose another groping factor for comparison"
    )
  } else {
    condition1 = levels(ninetails_data[[grouping_factor]])[1]
    condition2 = levels(ninetails_data[[grouping_factor]])[2]
  }

  # initial status code
  stats_code = codes_stats = "OK"
  # calculate group counts
  group_counts = ninetails_data %>%
    dplyr::group_by(!!!rlang::syms(c(grouping_factor))) %>%
    dplyr::count()

  stats <- NA

  if (base == "C") {
    count_column <- "counts_C"
  } else if (base == "G") {
    count_column <- "counts_G"
  } else if (base == "U") {
    count_column <- "counts_U"
  } else if (base == "all") {
    #ninetails_data <- ninetails_data %>% dplyr::mutate(counts_nonA=)
    count_column <- "counts_nonA"
  } else {
    stop(
      "Wrong non-A nucleotide defined. To compute statistics, please provide 'base' argument as character string (C, G, U or all)."
    )
  }

  if (nrow(group_counts) == 2) {
    if (group_counts[1, ]$n < min_reads) {
      if (group_counts[2, ]$n < min_reads) {
        stats_code = "B_LC"
      } else {
        stats_code = "G_LC"
      }
    } else if (group_counts[2, ]$n < min_reads) {
      stats_code = "G_LC"
    } else {
      options(scipen = 999)

      # summarize nonAs
      contingency_table <- ninetails::summarize_nonA(
        merged_nonA_tables = ninetails_data,
        summary_factors = grouping_factor,
        transcript_id_column = transcript_id_column
      ) %>%
        dplyr::select(
          !!rlang::sym(grouping_factor),
          counts_blank,
          !!rlang::sym(count_column)
        )
      contingency_table <- as.data.frame(contingency_table) # coerce tibble to df as setting names to tibble is deprecated
      row.names(contingency_table) <- contingency_table[[grouping_factor]] # set rownames
      contingency_table[[grouping_factor]] <- NULL # drop grouping col
      stats <- suppressWarnings(stats::fisher.test(contingency_table))$p.value
    }
  } else if (nrow(group_counts) == 1) {
    stats_code = "G_NA"
  } else if (nrow(group_counts) == 0) {
    stats_code = "B_NA"
  } else {
    stats_code = "ERR"
  }

  # create output
  stats <- tibble::tibble(
    p.value = stats,
    stats_code = as.character(stats_code)
  )

  return(stats)
}

stat_codes_list = list(
  OK = "OK",
  G1_NA = "GROUP1_NA",
  G2_NA = "GROUP2_NA",
  G1_LC = "G1_LOW_COUNT",
  G2_LC = "G2_LOW_COUNT",
  B_NA = "DATA FOR BOTH GROUPS NOT AVAILABLE",
  B_LC = "LOW COUNTS FOR BOTH GROUPS",
  G_LC = "LOW COUNT FOR ONE GROUP",
  G_NA = "DATA FOR ONE GROUP NOT AVAILABLE",
  ERR = "OTHER ERROR"
)


#' Perform Fisher's exact test per transcript with BH p-value adjustment
#'
#' Runs Fisher's exact test for each transcript in ninetails output data
#' and applies the Benjamini-Hochberg (BH) step-up procedure to control
#' the false discovery rate. This is a high-level wrapper for
#' \code{\link[stats]{fisher.test}} and \code{\link[stats]{p.adjust}} with
#' additional data wrangling features.
#'
#' @section Acknowledgements:
#' Inspired by the Nanotail package written and maintained by
#' Pawel Krawczyk (smaegol):
#' \url{https://github.com/LRB-IIMCB/nanotail/blob/dev/R/polya_stats.R}.
#' Many thanks to the developer of the original source code.
#'
#' @param ninetails_data Data frame. The output of
#'   \code{\link{merge_nonA_tables}} (merged tabular output containing read
#'   classification and non-A position data).
#'
#' @param transcript_id_column Character string. Column with transcript ID
#'   data. Default: \code{"ensembl_transcript_id_short"}; can be changed by
#'   the user.
#'
#' @param min_reads Numeric. Minimum number of reads representing a given
#'   transcript to include it in the analysis. Default: 0. Keep in mind that
#'   including many transcripts with low coverage increases the risk of
#'   rejecting the true null hypothesis (Benjamini-Hochberg procedure).
#'
#' @param min_nonA_reads Numeric. Minimum number of reads containing
#'   non-adenosine residues (sum of C, G, U) per transcript to include it
#'   in the analysis. This prevents considering too many observations as
#'   non-significant during p-value adjustment. Non-A-containing reads are
#'   typically a small fraction of the total pool, so additional filtering
#'   can provide more meaningful results. Default: 0.
#'
#' @param grouping_factor Character string. Name of the grouping variable.
#'   Default: \code{"sample_name"}.
#'
#' @param condition1 Character string. First level of \code{grouping_factor}
#'   to use for comparison. Required when \code{grouping_factor} has more
#'   than 2 levels.
#'
#' @param condition2 Character string. Second level of \code{grouping_factor}
#'   to use for comparison. Required when \code{grouping_factor} has more
#'   than 2 levels.
#'
#' @param alpha Numeric. Significance threshold for FDR. Default: 0.05.
#'
#' @param base Character string. Letter representing the non-A nucleotide
#'   for which statistics are computed. Accepted values: \code{"C"},
#'   \code{"G"}, \code{"U"}, or \code{"all"} (for all non-A residues
#'   combined). Default: \code{"C"}.
#'
#' @param ... Additional parameters passed to \code{\link{nonA_fisher}}
#'   (under development).
#'
#' @return A tibble with per-transcript test results, sorted by adjusted
#'   p-value. Columns include:
#'   \describe{
#'     \item{<transcript_id_column>}{Transcript identifier}
#'     \item{p.value}{Numeric. Raw p-value from Fisher's exact test}
#'     \item{stats_code}{Character. Descriptive status code (see
#'       \code{stat_codes_list} for full definitions)}
#'     \item{padj}{Numeric. BH-adjusted p-value}
#'     \item{significance}{Character. \code{"FDR<alpha"} if significant,
#'       \code{"NotSig"} otherwise}
#'   }
#'
#' @seealso
#' \code{\link{nonA_fisher}} for the per-transcript Fisher's exact test,
#' \code{\link{merge_nonA_tables}} for preparing the input,
#' \code{\link{summarize_nonA}} for contingency table computation
#'
#' @importFrom rlang :=
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' test <- ninetails::calculate_fisher(
#'   ninetails_data = merged_nonA_tables,
#'   transcript_id_column = "ensembl_transcript_id_short",
#'   min_reads = 100,
#'   min_nonA_reads = 10,
#'   grouping_factor = "sample_name",
#'   condition1 = "WT",
#'   condition2 = "KO",
#'   alpha = 0.05,
#'   base = "C"
#' )
#'
#' }
calculate_fisher <- function(ninetails_data,
                             transcript_id_column = "ensembl_transcript_id_short",
                             min_reads = 0,
                             min_nonA_reads = 0,
                             grouping_factor = "sample_name",
                             condition1 = NA,
                             condition2 = NA,
                             alpha = 0.05,
                             base = "C",
                             ...) {
  # Assertions
  if (missing(ninetails_data)) {
    stop(
      "Ninetails data are missing. Please provide a valid ninetails_data argument",
      call. = FALSE
    )
  }
  if (missing(transcript_id_column)) {
    stop(
      "Transcript_id_column is missing. Please provide a valid transcript_id_column argument",
      call. = FALSE
    )
  }
  if (missing(min_reads)) {
    stop(
      "Min_reads are missing. Please provide a valid min_reads argument",
      call. = FALSE
    )
  }
  if (missing(min_nonA_reads)) {
    stop(
      "Min_nonA_reads data are missing. Please provide a valid min_nonA_reads argument",
      call. = FALSE
    )
  }
  if (missing(base)) {
    stop(
      "Base definition is missing. Please provide a valid base argument",
      call. = FALSE
    )
  }

  assert_condition(
    is.numeric(min_reads),
    "Min_reads must be numeric. Please provide a valid argument."
  )
  assert_condition(
    is.numeric(min_nonA_reads),
    "Min_nonA_reads must be numeric. Please provide a valid argument."
  )
  assert_condition(
    is.numeric(alpha),
    "Alpha must be numeric. Please provide a valid argument."
  )

  if (!is.data.frame(ninetails_data) || nrow(ninetails_data) == 0) {
    stop(
      "Empty data frame provided as an input (ninetails_data). Please provide valid input"
    )
  }

  # if grouping factor has more than two levels
  if (length(levels(ninetails_data[[grouping_factor]])) > 2) {
    if (is.na(condition1) && is.na(condition2)) {
      #throw error when no conditions for comparison are specified
      stop(paste0(
        "grouping_factor ",
        grouping_factor,
        " has more than 2 levels. Please specify condtion1 and condition2 to select comparison pairs"
      ))
    } else {
      # filter input data leaving only specified conditions, dropping other factor levels
      assert_condition(
        condition1 %in% levels(ninetails_data[[grouping_factor]]),
        paste0(
          condition1,
          " is not a level of ",
          grouping_factor,
          " (grouping_factor)"
        )
      )
      assert_condition(
        condition2 %in% levels(ninetails_data[[grouping_factor]]),
        paste0(
          condition2,
          " is not a level of ",
          grouping_factor,
          " (grouping_factor)"
        )
      )
      assert_condition(
        condition2 != condition1,
        "condition2 should be different than condition1"
      )

      ninetails_data <- ninetails_data %>%
        dplyr::filter(
          !!rlang::sym(grouping_factor) %in% c(condition1, condition2)
        ) %>%
        dplyr::mutate() %>%
        droplevels()
    }
  } else if (length(levels(ninetails_data[[grouping_factor]])) == 1) {
    stop(
      "Only 1 level present for grouping factor. Choose another groping factor for comparison"
    )
  } else {
    condition1 = levels(ninetails_data[[grouping_factor]])[1]
    condition2 = levels(ninetails_data[[grouping_factor]])[2]
  }

  # filter out transcripts with not enough amount of non-A reads among the whole pool of reads:
  mod_summarized <- ninetails_data %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      sum_nonA = rowSums(dplyr::across(dplyr::starts_with('prediction_')))
    ) %>%
    dplyr::group_by(
      !!!rlang::syms(c(transcript_id_column, grouping_factor))
    ) %>%
    dplyr::summarise(
      dplyr::across(c(sum_nonA), list(counts = ~ sum(.x != 0))),
      .groups = 'drop'
    ) %>%
    dplyr::rename_with(
      ~ stringr::str_replace(.x, '^\\w+_(\\w+)_(\\w+)', '\\2_\\1'),
      3:dplyr::last_col()
    )
  # apply filtering criterion (minimal nonA read content)
  mod_summarized_filtered <- mod_summarized %>%
    dplyr::filter(counts_nonA >= min_nonA_reads)
  #extract filtered trans
  contig <- as.name(transcript_id_column)
  mod_summarized_filtered <- unique(mod_summarized_filtered[[contig]])

  ninetails_data <- ninetails_data %>%
    dplyr::filter(!!rlang::sym(contig) %in% mod_summarized_filtered)

  #ninetails_data <- ninetails_data[ninetails_data$contig %in% mod_summarized_filtered,]
  #ninetails_data <- ninetails_data[transcript_id_column %in% mod_summarized_filtered,]

  ninetails_data_stat <- ninetails_data %>%
    dplyr::mutate(transcript_id = get(c(transcript_id_column))) %>%
    dplyr::group_by(get(c(transcript_id_column))) %>%
    tidyr::nest()

  message("calculating statistics")
  ninetails_data_stat <- ninetails_data_stat %>%
    dplyr::mutate(
      stats = purrr::map(
        data,
        ninetails::nonA_fisher,
        grouping_factor = grouping_factor,
        min_reads = min_reads,
        base = base,
        transcript_id_column = transcript_id_column
      )
    ) %>%
    dplyr::select(-data) %>%
    tidyr::unnest(cols = c(stats)) %>%
    dplyr::rename(!!transcript_id_column := "get(c(transcript_id_column))")

  message("Finished")

  message("Adjusting p.value")
  ninetails_data_stat$padj <- stats::p.adjust(
    ninetails_data_stat$p.value,
    method = "BH"
  )

  # create significance factor
  ninetails_data_stat <- ninetails_data_stat %>%
    dplyr::mutate(
      significance = dplyr::case_when(
        is.na(padj) ~ "NotSig",
        (padj < alpha) ~ paste0("FDR<", alpha),
        TRUE ~ "NotSig"
      )
    )
  ninetails_data_stat$stats_code <- sapply(
    ninetails_data_stat$stats_code,
    FUN = function(x) {
      stat_codes_list[[x]]
    },
    simplify = "vector",
    USE.NAMES = FALSE
  ) %>%
    unlist()

  ninetails_data_stat <- ninetails_data_stat %>% dplyr::arrange(padj)

  return(ninetails_data_stat)
}
