#' Reads ninetails read_classes data frame from file.
#'
#' Imports a single \code{read_classes} output file produced by
#' \code{\link{check_tails_guppy}} (or other ninetails pipelines) into R.
#' Optionally attaches a sample identifier and parses GENCODE-style contig
#' names into Ensembl transcript IDs.
#'
#' @details
#' After loading, the function applies
#' \code{\link{correct_labels}} to ensure backward compatibility with the
#' label changes introduced in v0.9 (decorated/blank/unclassified). If the
#' \code{contig} column contains GENCODE-formatted identifiers, three
#' additional columns are created: \code{transcript} (gene symbol),
#' \code{ensembl_transcript_id_full} (with version), and
#' \code{ensembl_transcript_id_short} (without version).
#'
#' @section Acknowledgements:
#' Function based on \code{read_polya_single} from the NanoTail package
#' by P. Krawczyk (smaegol):
#' \url{https://github.com/LRB-IIMCB/nanotail/}.
#'
#' @param class_path Character string. Path to the ninetails
#'   \code{read_classes} output file (\code{.tsv}).
#'
#' @param sample_name Character string (optional, default \code{NA}). If
#'   specified, added as a factor column \code{sample_name}. Overwrites any
#'   existing \code{sample_name} column with a warning.
#'
#' @return A tibble containing read_classes predictions with the following
#'   additional columns (when GENCODE contigs are present):
#' \describe{
#'   \item{transcript}{Character. Gene symbol parsed from GENCODE contig.}
#'   \item{ensembl_transcript_id_full}{Character. Ensembl transcript ID
#'     with version number.}
#'   \item{ensembl_transcript_id_short}{Character. Ensembl transcript ID
#'     without version number.}
#' }
#'
#' @seealso \code{\link{read_class_multiple}} for batch loading,
#'   \code{\link{correct_labels}} for the backward-compatibility step,
#'   \code{\link{check_tails_guppy}} for the pipeline that produces the
#'   input file.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' class_path <- "/directory/with/ninetails/read_class_output.tsv"
#' class_data <- ninetails::read_class_single(class_path)
#'
#' }
read_class_single <- function(class_path, sample_name = NA) {
  # assertions
  if (missing(class_path)) {
    stop(
      "The path to class predictions (argument class_path) is missing",
      call. = FALSE
    )
  }

  assert_condition(
    !is.na(class_path) && nchar(class_path) > 0,
    paste("File path ", class_path, " is missing or invalid", sep = "")
  )

  #missing empty file
  assert_condition(
    file.exists(class_path) && file.info(class_path)$size > 0,
    paste("File ", class_path, " is empty or does not exist", sep = "")
  )

  # load class data
  message(paste0("Loading data from ", class_path))
  class_data <- vroom::vroom(class_path, show_col_types = FALSE) %>%
    dplyr::as_tibble()

  # resolve backcompatibility issues (from v.0.9)
  class_data <- ninetails::correct_labels(class_data)

  if (!is.na(sample_name)) {
    # set sample_name (if was set)
    if (!"sample_name" %in% colnames(class_data)) {
      warning(
        "The sample_name was provided in the input file. Overwriting with the provided one"
      )
    }
    class_data$sample_name = sample_name
    class_data$sample_name <- as.factor(class_data$sample_name)
  }

  ## correct annotation <- if gencode format: create new columns with ensembl IDs
  # else created columns could be easily dropped
  # this code chunk was originally written by PaweÅ‚ Krawczyk (smaegol) & incorporated in NanoTail package
  transcript_names <- gsub(
    ".*?\\|.*?\\|.*?\\|.*?\\|.*?\\|(.*?)\\|.*",
    "\\1",
    class_data$contig
  )
  class_data$transcript <- transcript_names
  ensembl_transcript_ids <- gsub("^(.*?)\\|.*\\|.*", "\\1", class_data$contig)
  ensembl_transcript_ids_short <- gsub(
    "(.*)\\..*",
    "\\1",
    ensembl_transcript_ids
  ) # without version number
  class_data$ensembl_transcript_id_full <- ensembl_transcript_ids
  class_data$ensembl_transcript_id_short <- ensembl_transcript_ids_short

  return(class_data)
}


#' Reads multiple ninetails read_classes outputs at once.
#'
#' Batch-loads any number of \code{read_classes} prediction files with a
#' single invocation, attaching sample-level metadata from the provided
#' \code{samples_table}. Each file is loaded via
#' \code{\link{read_class_single}} and the results are combined into a
#' single long-format tibble.
#'
#' @section Acknowledgements:
#' Function based on \code{read_polya_multiple} from the NanoTail package
#' by P. Krawczyk (smaegol):
#' \url{https://github.com/LRB-IIMCB/nanotail/}.
#'
#' @param samples_table Data frame or tibble containing sample metadata and
#'   file paths. Must have at least two columns:
#' \describe{
#'   \item{class_path}{Character. Path to the \code{read_classes}
#'     prediction file for each sample.}
#'   \item{sample_name}{Character/factor. Unique sample identifier.}
#' }
#'   Additional metadata columns (e.g. \code{group}, \code{condition})
#'   are preserved and propagated to the output.
#'
#' @param ... Additional parameters passed to
#'   \code{\link{read_class_single}} (currently accepts
#'   \code{sample_name}).
#'
#' @return A tibble containing read_classes data for all specified
#'   samples, with metadata from \code{samples_table} stored as
#'   additional columns.
#'
#' @seealso \code{\link{read_class_single}} for the per-file loader,
#'   \code{\link{read_residue_multiple}} for the analogous residue data
#'   loader,
#'   \code{\link{merge_nonA_tables}} for combining class and residue
#'   outputs.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' samples_table <- data.frame(
#'   class_path = c(path_1, path_2, path_3, path_4),
#'   sample_name = c("wt_1", "mut_1", "wt_2", "mut_2"),
#'   group = c("wt", "mut", "wt", "mut"))
#'
#' classes_data <- ninetails::read_class_multiple(samples_table)
#'
#' }
read_class_multiple <- function(samples_table, ...) {

  if (missing(samples_table)) {
    stop("Samples table argument is missing", call. = FALSE)
  }

  if (!is.data.frame(samples_table) || nrow(samples_table) == 0) {
    stop(
      "Empty data frame provided as an input (samples_table). Please provide samples_table describing data to load"
    )
  }

  assert_condition(
    "class_path" %in% colnames(samples_table),
    "The samples_table should contain at least class_path and sample_name columns"
  )
  assert_condition(
    "sample_name" %in% colnames(samples_table),
    "The samples_table should contain at least class_path and sample_name columns"
  )

  samples_data <- samples_table %>%
    tibble::as_tibble() %>%
    dplyr::mutate_if(is.character, as.factor) %>%
    dplyr::mutate(class_path = as.character(class_path)) %>%
    dplyr::group_by(sample_name) %>%
    dplyr::mutate(
      class_contents = purrr::map(class_path, function(x) {
        ninetails::read_class_single(x)
      })
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-class_path)

  class_data <- tidyr::unnest(samples_data, cols = c(class_contents))

  return(class_data)
}

#' Counts read classes found in a read_classes data frame produced by the
#' ninetails pipeline.
#'
#' Tabulates the prediction information contained in the \code{read_classes}
#' data frame. Counts can be computed at two levels of granularity
#' (detailed or crude) and optionally stratified by a grouping variable.
#'
#' @details
#' When \code{detailed = TRUE}, counts are based on the \code{comments}
#' column, which carries fine-grained labels such as \code{"YAY"} (tail
#' with non-A residues detected). When \code{detailed = FALSE}, counts
#' are based on the \code{class} column with three broad categories:
#' \code{"decorated"}, \code{"blank"}, and \code{"unclassified"}.
#'
#' @param class_data Data frame or tibble containing read_classes
#'   predictions produced by the ninetails pipeline.
#'
#' @param grouping_factor Character string (default \code{NA}). Name of a
#'   column in \code{class_data} to use as a grouping variable
#'   (e.g. \code{"sample_name"}).
#'
#' @param detailed Logical \code{[TRUE]}. If \code{TRUE}, counts are
#'   provided based on the \code{comments} column (fine-grained). If
#'   \code{FALSE}, counts are provided based on the \code{class} column
#'   (crude: decorated / blank / unclassified only).
#'
#' @return A tibble with columns for the grouping variable (if provided),
#'   the classification label (\code{comments} or \code{class}), and
#'   \code{n} (the count).
#'
#' @seealso \code{\link{read_class_single}} and
#'   \code{\link{read_class_multiple}} for loading class data,
#'   \code{\link{count_residues}} for the analogous residue-level counts,
#'   \code{\link{merge_nonA_tables}} for combining class and residue
#'   data.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' class_counted <- ninetails::count_class(
#'   class_data = out[[1]],
#'   grouping_factor = NA,
#'   detailed = TRUE)
#'
#' }
count_class <- function(class_data, grouping_factor = NA, detailed = TRUE) {

  #assertions
  if (!is.data.frame(class_data) || nrow(class_data) == 0) {
    stop(
      "Empty data frame provided as an input (class_data). Please provide valid input"
    )
  }

  if (detailed == TRUE) {
    if (!is.na(grouping_factor)) {
      assert_condition(
        grouping_factor %in% colnames(class_data),
        paste0(grouping_factor, " is not a column of input dataset")
      )
      class_counts <- class_data %>%
        dplyr::mutate(
          comments = forcats::fct_relevel(comments, "YAY", after = Inf)
        ) %>%
        dplyr::group_by(!!rlang::sym(grouping_factor), comments) %>%
        dplyr::count()
    } else {
      class_counts <- class_data %>%
        dplyr::mutate(
          comments = forcats::fct_relevel(comments, "YAY", after = Inf)
        ) %>%
        dplyr::group_by(comments) %>%
        dplyr::count()
    }
  } else {
    if (!is.na(grouping_factor)) {
      assert_condition(
        grouping_factor %in% colnames(class_data),
        paste0(grouping_factor, " is not a column of input dataset")
      )
      class_counts <- class_data %>%
        dplyr::mutate(
          class = forcats::fct_relevel(class, "decorated", after = Inf)
        ) %>%
        dplyr::group_by(!!rlang::sym(grouping_factor), class) %>%
        dplyr::count()
    } else {
      class_counts <- class_data %>%
        dplyr::mutate(
          class = forcats::fct_relevel(class, "decorated", after = Inf)
        ) %>%
        dplyr::group_by(class) %>%
        dplyr::count()
    }
  }

  return(class_counts)
}


#' Reads ninetails nonadenosine_residues data from file.
#'
#' Imports a single \code{nonadenosine_residues} output file produced by
#' \code{\link{check_tails_guppy}} (or other ninetails pipelines) into R.
#' Optionally attaches a sample identifier and parses GENCODE-style contig
#' names into Ensembl transcript IDs.
#'
#' @details
#' If the \code{contig} column contains GENCODE-formatted identifiers,
#' three additional columns are created: \code{transcript},
#' \code{ensembl_transcript_id_full}, and
#' \code{ensembl_transcript_id_short} (see \code{\link{read_class_single}}
#' for the same logic applied to class data).
#'
#' @section Acknowledgements:
#' Function based on \code{read_polya_single} from the NanoTail package
#' by P. Krawczyk (smaegol):
#' \url{https://github.com/LRB-IIMCB/nanotail/}.
#'
#' @param residue_path Character string. Path to the ninetails
#'   \code{nonadenosine_residues} output file (\code{.tsv}).
#'
#' @param sample_name Character string (optional, default \code{NA}). If
#'   specified, added as a factor column \code{sample_name}. Overwrites any
#'   existing \code{sample_name} column with a warning.
#'
#' @return A tibble containing non-A residue predictions with the
#'   following additional columns (when GENCODE contigs are present):
#' \describe{
#'   \item{transcript}{Character. Gene symbol parsed from GENCODE contig.}
#'   \item{ensembl_transcript_id_full}{Character. Ensembl transcript ID
#'     with version number.}
#'   \item{ensembl_transcript_id_short}{Character. Ensembl transcript ID
#'     without version number.}
#' }
#'
#' @seealso \code{\link{read_residue_multiple}} for batch loading,
#'   \code{\link{read_class_single}} for the analogous class data loader,
#'   \code{\link{check_tails_guppy}} for the pipeline that produces the
#'   input file.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' residue_path <- "/directory/with/ninetails/nonadenosine_residues_output.tsv"
#' residue_data <- ninetails::read_residue_single(residue_path)
#'
#' }
read_residue_single <- function(residue_path, sample_name = NA) {

  #assertions
  assert_condition(
    !is.na(residue_path) && nchar(residue_path) > 0,
    "Input is either missing or an empty string. Please provide a residue_path as a string"
  )

  assert_condition(
    file.exists(residue_path) && file.info(residue_path)$size > 0,
    paste0("Input is either missing or an empty file: ", residue_path)
  )

  # load the data
  message(paste0("Loading non-A residue data from ", residue_path))
  residue_data <- vroom::vroom(residue_path, show_col_types = FALSE) %>%
    dplyr::as_tibble()

  if (!is.na(sample_name)) {
    # set sample_name (if was set)
    if (!"sample_name" %in% colnames(residue_data)) {
      warning(
        "sample_name was provided in the input file. Overwriting with the provided one"
      )
    }
    residue_data$sample_name = sample_name
    residue_data$sample_name <- as.factor(residue_data$sample_name)
  }

  ## correct annotation <- if gencode format: create new columns with ensembl IDs
  # else created columns could be easily dropped
  # this code chunk was originally written by PaweÅ‚ Krawczyk (smaegol) & incorporated in NanoTail package

  transcript_names <- gsub(
    ".*?\\|.*?\\|.*?\\|.*?\\|.*?\\|(.*?)\\|.*",
    "\\1",
    residue_data$contig
  )
  residue_data$transcript <- transcript_names
  ensembl_transcript_ids <- gsub("^(.*?)\\|.*\\|.*", "\\1", residue_data$contig)
  ensembl_transcript_ids_short <- gsub(
    "(.*)\\..*",
    "\\1",
    ensembl_transcript_ids
  ) # without version number
  residue_data$ensembl_transcript_id_full <- ensembl_transcript_ids
  residue_data$ensembl_transcript_id_short <- ensembl_transcript_ids_short

  return(residue_data)
}


#' Reads multiple ninetails nonadenosine_residues outputs at once.
#'
#' Batch-loads any number of \code{nonadenosine_residues} prediction files
#' with a single invocation, attaching sample-level metadata from the
#' provided \code{samples_table}. Each file is loaded via
#' \code{\link{read_residue_single}} and the results are combined into a
#' single long-format tibble.
#'
#' @section Acknowledgements:
#' Function based on \code{read_polya_multiple} from the NanoTail package
#' by P. Krawczyk (smaegol):
#' \url{https://github.com/LRB-IIMCB/nanotail/}.
#'
#' @param samples_table Data frame or tibble containing sample metadata and
#'   file paths. Must have at least two columns:
#' \describe{
#'   \item{residue_path}{Character. Path to the
#'     \code{nonadenosine_residues} prediction file for each sample.}
#'   \item{sample_name}{Character/factor. Unique sample identifier.}
#' }
#'   Additional metadata columns are preserved and propagated to the
#'   output.
#'
#' @param ... Additional parameters passed to
#'   \code{\link{read_residue_single}} (currently accepts
#'   \code{sample_name}).
#'
#' @return A tibble containing non-A residue data for all specified
#'   samples, with metadata from \code{samples_table} stored as
#'   additional columns.
#'
#' @seealso \code{\link{read_residue_single}} for the per-file loader,
#'   \code{\link{read_class_multiple}} for the analogous class data
#'   loader,
#'   \code{\link{merge_nonA_tables}} for combining class and residue
#'   outputs.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' samples_table <- data.frame(
#'   residue_path = c(path_1, path_2, path_3, path_4),
#'   sample_name = c("wt_1", "mut_1", "wt_2", "mut_2"),
#'   group = c("wt", "mut", "wt", "mut"))
#'
#' residues_data <- ninetails::read_residue_multiple(samples_table)
#'
#' }
read_residue_multiple <- function(samples_table, ...) {

  #assertions
  if (missing(samples_table)) {
    stop("The samples_table argument is missing", call. = FALSE)
  }

  if (!is.data.frame(samples_table) || nrow(samples_table) == 0) {
    stop(
      "Empty data frame provided as an input (samples_table). Please provide valid input"
    )
  }

  assert_condition(
    "residue_path" %in% colnames(samples_table),
    "The samples_table should contain at least residue_path and sample_name columns."
  )
  assert_condition(
    "sample_name" %in% colnames(samples_table),
    "The samples_table should contain at least residue_path and sample_name columns."
  )

  samples_data <- samples_table %>%
    tibble::as_tibble() %>%
    dplyr::mutate_if(is.character, as.factor) %>%
    dplyr::mutate(residue_path = as.character(residue_path)) %>%
    dplyr::group_by(sample_name) %>%
    dplyr::mutate(
      residue_contents = purrr::map(residue_path, function(x) {
        ninetails::read_residue_single(x)
      })
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-residue_path)
  residue_data <- tidyr::unnest(samples_data, cols = c(residue_contents))

  return(residue_data)
}

#' Counts non-A residues found in a nonadenosine_residues data frame
#' produced by the ninetails pipeline.
#'
#' Tabulates the \code{prediction} column of the
#' \code{nonadenosine_residues} data frame, counting occurrences of each
#' non-A nucleotide type (C, G, U). Counts represent total hits across
#' all reads, \emph{not} per-read summaries (a single read may contribute
#' multiple hits).
#'
#' @param residue_data Data frame or tibble containing non-A residue
#'   predictions produced by the ninetails pipeline.
#'
#' @param grouping_factor Character string (default \code{NA}). Name of a
#'   column in \code{residue_data} to use as a grouping variable
#'   (e.g. \code{"sample_name"}).
#'
#' @return A tibble with columns for the grouping variable (if provided),
#'   \code{prediction} (nucleotide type), and \code{n} (the count).
#'
#' @seealso \code{\link{count_class}} for counting read-level classes,
#'   \code{\link{read_residue_single}} and
#'   \code{\link{read_residue_multiple}} for loading residue data,
#'   \code{\link{summarize_nonA}} for transcript-level summaries.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' residue_counted <- ninetails::count_residues(
#'   residue_data = out[[2]],
#'   grouping_factor = NA)
#'
#' }
#'
count_residues <- function(residue_data, grouping_factor = NA) {

  # assertions
  if (!is.data.frame(residue_data) || nrow(residue_data) == 0) {
    stop(
      "Empty data frame provided as an input (residue_data). Please provide valid input"
    )
  }

  if (!is.na(grouping_factor)) {
    assert_condition(
      grouping_factor %in% colnames(residue_data),
      paste0(grouping_factor, " is not a column of input dataset")
    )

    residue_counts <- residue_data %>%
      dplyr::mutate(
        prediction = forcats::fct_relevel(prediction, "U", after = Inf)
      ) %>%
      dplyr::group_by(!!rlang::sym(grouping_factor), prediction) %>%
      dplyr::count()
  } else {
    residue_counts <- residue_data %>%
      dplyr::mutate(
        prediction = forcats::fct_relevel(prediction, "U", after = Inf)
      ) %>%
      dplyr::group_by(prediction) %>%
      dplyr::count()
  }

  return(residue_counts)
}


#' Reshapes nonadenosine_residues data frame to wide format.
#'
#' Pivots the \code{prediction} column of the
#' \code{nonadenosine_residues} data frame into three separate columns
#' (\code{prediction_C}, \code{prediction_G}, \code{prediction_U}), each
#' containing the count of respective non-A residues per read. An
#' additional \code{nonA_residues} column provides a CIGAR-like summary
#' string listing all non-A positions from 5' to 3', separated by
#' \code{":"}.
#'
#' @details
#' This function supports data from both the legacy Guppy/nanopolish
#' pipeline and the Dorado pipeline. The \code{qc_tag} column may be
#' either character (Guppy: \code{"PASS"}, \code{"SUFFCLIP"}, etc.) or
#' numeric (Dorado: MAPQ score). Both types are handled transparently
#' and preserved in the output.
#'
#' @param residue_data Data frame or tibble containing non-A residue
#'   predictions produced by the ninetails pipeline (filename typically
#'   ends with \code{_nonadenosine_residues.tsv}). Required columns:
#'   \code{readname}, \code{group}, \code{prediction}, \code{est_nonA_pos}.
#'   The \code{qc_tag} column may be character (Guppy/nanopolish) or
#'   numeric (Dorado MAPQ).
#'
#' @return A tibble in wide format with one row per read. In addition to
#'   the original columns (minus \code{est_nonA_pos} and
#'   \code{prediction}), the following are added:
#' \describe{
#'   \item{prediction_C}{Integer. Number of C residues detected in the
#'     read.}
#'   \item{prediction_G}{Integer. Number of G residues detected in the
#'     read.}
#'   \item{prediction_U}{Integer. Number of U residues detected in the
#'     read.}
#'   \item{nonA_residues}{Character. CIGAR-like string summarising all
#'     non-A positions (e.g. \code{"C12:G25:U40"}).}
#' }
#'
#' @seealso \code{\link{merge_nonA_tables}} which calls this function
#'   internally,
#'   \code{\link{read_residue_single}} for loading residue data.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' spread_table <- ninetails::spread_nonA_residues(
#'   residue_data = residue_data)
#'
#' }
spread_nonA_residues <- function(residue_data) {


  # Assertions
  if (missing(residue_data)) {
    stop(
      "A dataframe with non-A residue data is missing. Please provide a valid residue_data argument",
      call. = FALSE
    )
  }

  if (!is.data.frame(residue_data) || nrow(residue_data) == 0) {
    stop(
      "Empty data frame provided as an input (residue_data). Please provide valid input"
    )
  }

  # Validate required columns
  required_cols <- c("readname", "group", "prediction", "est_nonA_pos")
  missing_cols <- setdiff(required_cols, colnames(residue_data))
  if (length(missing_cols) > 0) {
    stop(
      paste0("Required columns missing from residue_data: ",
             paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  # Create column for non-A residue summary (cigar-like string)
  cigar <- residue_data %>%
    dplyr::group_by(group, readname) %>%
    dplyr::arrange(est_nonA_pos, .by_group = TRUE) %>%
    dplyr::summarise(
      nonA_residues = paste0(prediction, est_nonA_pos, collapse = ':'),
      .groups = 'drop'
    )

  # Create contingency table with C, G, U counts per read
  contingency <- residue_data %>%
    dplyr::select(-est_nonA_pos) %>%
    tidyr::pivot_wider(
      names_from = prediction,
      names_sort = TRUE,
      names_prefix = 'prediction_',
      values_from = prediction,
      values_fn = length,
      values_fill = 0
    )

  # Merge both tables
  spread_table <- contingency %>%
    dplyr::left_join(cigar, by = c("readname", "group"))

  return(spread_table)
}



#' Merges ninetails tabular outputs (read classes and nonadenosine residue
#' data) to produce one concise table.
#'
#' Combines the \code{read_classes} and \code{nonadenosine_residues} data
#' frames into a single wide-format tibble with one row per read. The
#' \code{prediction} column from the residue data is spread into three
#' columns (\code{prediction_C}, \code{prediction_G}, \code{prediction_U})
#' via \code{\link{spread_nonA_residues}}, and a CIGAR-like
#' \code{nonA_residues} summary column is appended.
#'
#' @details
#' Unclassified reads are excluded before merging. The quality-tag filter
#' behaviour depends on the pipeline that produced the input data:
#'
#' \strong{Guppy/nanopolish pipeline} (character \code{qc_tag}):
#' The \code{pass_only} parameter controls filtering. When \code{TRUE},
#' only reads tagged as \code{"PASS"} are retained. When \code{FALSE},
#' reads tagged as \code{"PASS"} or \code{"SUFFCLIP"} are retained.
#'
#' \strong{Dorado pipeline} (numeric \code{qc_tag}):
#' The \code{qc_tag} column contains MAPQ scores. Filtering is based on
#' mapping quality: only reads with \code{qc_tag > 0} are retained. The
#' \code{pass_only} parameter is ignored for numeric \code{qc_tag}.
#'
#' After the full join, any \code{NA} values in numeric columns are
#' replaced with 0.
#'
#' @section Pipeline Compatibility:
#' This function automatically detects whether data originates from the
#' legacy Guppy/nanopolish pipeline (character \code{qc_tag}) or the
#' Dorado pipeline (numeric \code{qc_tag}) and applies appropriate
#' filtering logic. Both \code{class_data} and \code{residue_data} must
#' originate from the same pipeline to ensure consistent \code{qc_tag}
#' types during the merge operation.
#'
#' @param class_data Data frame or tibble containing read_classes
#'   predictions produced by the ninetails pipeline. The \code{qc_tag}
#'   column may be character (Guppy/nanopolish: \code{"PASS"},
#'   \code{"SUFFCLIP"}, etc.) or numeric (Dorado: MAPQ score).
#'
#' @param residue_data Data frame or tibble containing non-A residue
#'   predictions produced by the ninetails pipeline. The \code{qc_tag}
#'   column type must match that of \code{class_data}.
#'
#' @param pass_only Logical \code{[TRUE]}. Applies only when
#'   \code{qc_tag} is character (Guppy/nanopolish pipeline). If
#'   \code{TRUE}, only reads tagged as \code{"PASS"} are included. If
#'   \code{FALSE}, reads tagged as \code{"PASS"} or \code{"SUFFCLIP"}
#'   are included. \strong{Ignored when \code{qc_tag} is numeric}
#'   (Dorado pipeline), in which case MAPQ > 0 filtering is applied
#'   instead.
#'
#' @return A tibble with summarised information from both ninetails
#'   outputs: all columns from \code{class_data} plus
#'   \code{prediction_C}, \code{prediction_G}, \code{prediction_U}, and
#'   \code{nonA_residues}.
#'
#' @seealso \code{\link{spread_nonA_residues}} for the reshaping step,
#'   \code{\link{summarize_nonA}} for transcript-level summaries from the
#'   merged table,
#'   \code{\link{calculate_fisher}} for statistical testing on the merged
#'   table,
#'   \code{\link{read_class_single}} and \code{\link{read_residue_single}}
#'   for loading the input data.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # Guppy/nanopolish data (character qc_tag)
#' merged_tables <- ninetails::merge_nonA_tables(
#'   class_data = class_data,
#'   residue_data = residue_data,
#'   pass_only = TRUE)
#'
#' # Dorado data (numeric qc_tag) - pass_only is ignored
#' merged_tables <- ninetails::merge_nonA_tables(
#'   class_data = dorado_class_data,
#'   residue_data = dorado_residue_data,
#'   pass_only = TRUE)
#'
#' }
merge_nonA_tables <- function(class_data,
                              residue_data,
                              pass_only = TRUE) {

  # Assertions
  if (missing(class_data)) {
    stop(
      "A dataframe with class data is missing. Please provide a valid class_data argument",
      call. = FALSE
    )
  }
  if (missing(residue_data)) {
    stop(
      "A dataframe with non-A residue data is missing. Please provide a valid residue_data argument",
      call. = FALSE
    )
  }

  if (!is.data.frame(class_data) || nrow(class_data) == 0) {
    stop(
      "Empty data frame provided as an input (class_data). Please provide valid input"
    )
  }

  if (!is.data.frame(residue_data) || nrow(residue_data) == 0) {
    stop(
      "Empty data frame provided as an input (residue_data). Please provide valid input"
    )
  }

  assert_condition(
    is.logical(pass_only),
    "Please provide TRUE/FALSE values for pass_only parameter"
  )

  # Validate qc_tag column exists
  if (!"qc_tag" %in% colnames(class_data)) {
    stop(
      "Column 'qc_tag' is missing from class_data. Please provide valid ninetails output.",
      call. = FALSE
    )
  }

  # Drop all unclassified reads

  class_data <- class_data[!(class_data$class == "unclassified"), ]

  # Filter class_data according to qc_tag type

  # Character qc_tag: Guppy/nanopolish pipeline ("PASS", "SUFFCLIP", etc.)
  # Numeric qc_tag: Dorado pipeline (MAPQ score)
  if (is.character(class_data$qc_tag)) {
    # Guppy/nanopolish pipeline: filter by qc_tag category
    if (pass_only == TRUE) {
      class2 <- class_data[class_data$qc_tag == "PASS", ]
    } else {
      class2 <- class_data[class_data$qc_tag %in% c("PASS", "SUFFCLIP"), ]
    }
  } else if (is.numeric(class_data$qc_tag)) {
    # Dorado pipeline: filter by MAPQ > 0
    # pass_only parameter is ignored for numeric qc_tag
    class2 <- class_data[class_data$qc_tag > 0, ]
  } else {
    stop(
      paste0("Unexpected qc_tag type: ", class(class_data$qc_tag),
             ". Expected character (Guppy/nanopolish) or numeric (Dorado)."),
      call. = FALSE
    )
  }

  # Check if filtering resulted in empty data frame
  if (nrow(class2) == 0) {
    warning(
      "No reads passed quality filtering. Returning empty merged table.",
      call. = FALSE
    )
  }

  # Spread residue_data
  spread_table <- ninetails::spread_nonA_residues(residue_data)

  # Merge the data (spreaded residue + classes)
  # Join on all common columns
  merged_tables <- class2 %>% dplyr::full_join(spread_table)

  # Replace NA with 0 in numeric columns
  # tidyselect::where() issue solved as in

  # https://stackoverflow.com/questions/62459736/how-do-i-use-tidyselect-where-in-a-custom-package
  merged_nonA_tables <- merged_tables %>%
    dplyr::mutate(dplyr::across(
      tidyselect::vars_select_helpers$where(is.numeric),
      ~ ifelse(is.na(.), 0, .)
    ))

  return(merged_nonA_tables)
}



#' Produces summary table of non-A occurrences within an analyzed dataset.
#'
#' Creates a per-transcript summary table with read counts, non-A residue
#' counts and hits, and poly(A) tail length statistics, grouped by
#' user-defined factors (e.g. sample, condition).
#'
#' @details
#' The distinction between \strong{counts} and \strong{hits}:
#' \itemize{
#'   \item \strong{counts} — number of reads containing at least one
#'     occurrence of a given non-A residue type.
#'   \item \strong{hits} — total number of occurrences of a given non-A
#'     residue type across all reads (a single read may contribute
#'     multiple hits).
#' }
#'
#' @param merged_nonA_tables Data frame or tibble. Output of
#'   \code{\link{merge_nonA_tables}}.
#'
#' @param summary_factors Character string or vector of strings. Column
#'   name(s) used for grouping (default: \code{"group"}).
#'
#' @param transcript_id_column Character string. Column containing the
#'   transcript identifier (default:
#'   \code{"ensembl_transcript_id_short"}, as added during data
#'   pre-processing; can be changed by the user).
#'
#' @return A tibble with one row per transcript per group, containing:
#' \describe{
#'   \item{polya_median}{Numeric. Median poly(A) tail length for the
#'     transcript.}
#'   \item{polya_mean}{Numeric. Mean poly(A) tail length for the
#'     transcript.}
#'   \item{counts_total}{Integer. Total number of reads mapped to the
#'     transcript.}
#'   \item{counts_blank}{Integer. Number of reads with no non-A
#'     residues.}
#'   \item{counts_nonA / hits_nonA}{Integer. Reads with any non-A /
#'     total non-A hits.}
#'   \item{counts_C / hits_C}{Integer. Reads with C / total C hits.}
#'   \item{counts_G / hits_G}{Integer. Reads with G / total G hits.}
#'   \item{counts_U / hits_U}{Integer. Reads with U / total U hits.}
#' }
#'
#' @seealso \code{\link{merge_nonA_tables}} for preparing the input,
#'   \code{\link{calculate_fisher}} for statistical testing on merged
#'   data,
#'   \code{\link{annotate_with_biomart}} for adding gene-level
#'   annotation.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' summarized <- ninetails::summarize_nonA(
#'   merged_nonA_tables = merged_nonA_tables,
#'   summary_factors = "group",
#'   transcript_id_column = "ensembl_transcript_id_short")
#'
#' }
summarize_nonA <- function(merged_nonA_tables,
                           summary_factors = c("group"),
                           transcript_id_column = c("ensembl_transcript_id_short")) {

  #Assertions
  if (missing(merged_nonA_tables)) {
    stop(
      "Ninetails' merged_nonA_tables output is missing. Please provide a valid merged_nonA_tables argument",
      call. = FALSE
    )
  }

  if (!is.data.frame(merged_nonA_tables) || nrow(merged_nonA_tables) == 0) {
    stop(
      "Empty data frame provided as an input (merged_nonA_tables). Please provide valid input"
    )
  }

  assert_condition(
    is.character(summary_factors),
    "Non-character argument is not alowed for `summary_factors`. Please provide either string or vector of strings"
  )
  assert_condition(
    all(summary_factors %in% colnames(merged_nonA_tables)),
    "Non-existent column name provided as the argument (summary_factors)"
  )

  # this function is slow; TODO: rethink the syntax
  # this syntax is slightly faster than without creating new variable according to microbenchmark
  nonA_data_summarized <- merged_nonA_tables %>%
    # drop all previous grouping vars
    dplyr::ungroup() %>%
    # count all nonA reads
    dplyr::mutate(
      sum_nonA = rowSums(dplyr::across(dplyr::starts_with('prediction_')))
    ) %>%
    # group by provided vars
    dplyr::group_by(
      !!!rlang::syms(c(transcript_id_column, summary_factors))
    ) %>%
    #provide summaries
    dplyr::summarise(
      polya_median = stats::median(polya_length), #add median polya length
      polya_mean = mean(polya_length), # add mean polya length
      counts_total = dplyr::n(), # add transcript count
      #summarize counts (number of reads with given non-A residues)
      #and hits (occurences of given residues)
      dplyr::across(
        c(sum_nonA, prediction_C, prediction_G, prediction_U),
        list(counts = ~ sum(.x != 0), hits = ~ sum(.x))
      ),
      .groups = 'drop'
    ) %>%
    # rename columns according to desired convention
    dplyr::rename_with(
      ~ stringr::str_replace(.x, '^\\w+_(\\w+)_(\\w+)', '\\2_\\1'),
      4:dplyr::last_col()
    )

  # add counts of blank reads - it has to be assigned to new var, because otherwise throws an error
  zeromod_summarized <- merged_nonA_tables %>%
    dplyr::ungroup() %>%
    dplyr::group_by(
      !!!rlang::syms(c(transcript_id_column, summary_factors))
    ) %>%
    dplyr::summarise(
      counts_blank = sum(dplyr::if_all(
        tidyselect::starts_with('prediction_'),
        ~ .x == 0
      )),
      .groups = 'drop'
    )

  nonA_data_summarized <- nonA_data_summarized %>%
    # add counts of blank reads - this syntax is because otherwise it throws an error
    # also: avoiding assigning new variable to temporary table
    dplyr::left_join(
      zeromod_summarized,
      by = c(transcript_id_column, summary_factors)
    ) %>%
    # move blank reads count near total counts - for table clarity
    dplyr::relocate(counts_blank, .after = counts_total)

  return(nonA_data_summarized)
}


#' Aggregates nanopolish polya quality control information.
#'
#' Returns read counts for each \code{qc_tag} category (e.g.
#' \code{"PASS"}, \code{"SUFFCLIP"}, \code{"NOREGION"}, etc.), optionally
#' stratified by a user-defined grouping variable.
#'
#' @details
#' \strong{Caution:} do not use the output of
#' \code{\link{merge_nonA_tables}} as input, because that table is already
#' filtered for quality. Use the raw \code{class_data} output from the
#' ninetails pipeline instead.
#'
#' @section Acknowledgements:
#' This is the ninetails implementation of the
#' \code{get_nanopolish_processing_info} function originally written by
#' P. Krawczyk (smaegol) and incorporated in the NanoTail package:
#' \url{https://github.com/LRB-IIMCB/nanotail/blob/master/R/polya_stats.R}.
#' Variable names were adjusted to match the ninetails naming convention.
#'
#' @param class_data Data frame or tibble containing the raw
#'   \code{class_data} output from the ninetails pipeline (not the merged
#'   output).
#'
#' @param grouping_factor Character string (default \code{NA}). Variable
#'   used for grouping (e.g. \code{"sample_name"}).
#'
#' @return A tibble with columns for the grouping variable (if provided),
#'   \code{qc_tag}, and \code{n} (the count).
#'
#' @seealso \code{\link{read_class_single}} for loading class data,
#'   \code{\link{count_class}} for counting read classification labels.
#'
#' @export
#'
nanopolish_qc <- function(class_data, grouping_factor = NA) {

  #assertions
  if (!is.data.frame(class_data) || nrow(class_data) == 0) {
    stop(
      "Empty data frame provided as an input (class_data). Please provide valid input"
    )
  }

  if (!is.na(grouping_factor)) {
    assert_condition(
      grouping_factor %in% colnames(class_data),
      paste0(grouping_factor, " is not a column of input dataset")
    )
    processing_info <- class_data %>%
      dplyr::mutate(
        qc_tag = forcats::fct_relevel(qc_tag, "PASS", after = Inf)
      ) %>%
      dplyr::group_by(!!rlang::sym(grouping_factor), qc_tag) %>%
      dplyr::count()
  } else {
    processing_info <- class_data %>%
      dplyr::mutate(
        qc_tag = forcats::fct_relevel(qc_tag, "PASS", after = Inf)
      ) %>%
      dplyr::group_by(qc_tag) %>%
      dplyr::count()
  }

  return(processing_info)
}


#' Marks uncertain positions of non-A residues in ninetails output data.
#'
#' Annotates each non-A residue position with a quality flag
#' (\code{qc_pos}) indicating whether the position is likely genuine
#' (\code{"Y"}) or a potential nanopolish segmentation artefact
#' (\code{"N"}). The assessment is based on quantile thresholds of
#' poly(A) tail length, modal non-A position per transcript, and
#' optional species-specific whitelists of transcripts with A-rich
#' 3' UTRs.
#'
#' @details
#' Nanopolish segmentation can misidentify nucleotides from A-rich
#' 3' UTR regions as part of the poly(A) tail. For such transcripts, a
#' peak of non-A positions accumulates near the transcript body. This
#' function flags those positions as ambiguous using a combination of:
#' \itemize{
#'   \item The modal position of non-A residues per transcript.
#'   \item The 0.05 quantile of poly(A) tail lengths per transcript
#'     (segmentation error boundary).
#'   \item The 0.05 quantile of non-A positions per transcript and
#'     prediction type.
#'   \item Species-specific whitelists of transcripts with hybrid tails
#'     (3' UTRs with > 80\% A in the last 20 positions).
#' }
#'
#' @param class_data Data frame or tibble containing read_classes
#'   predictions from the ninetails pipeline.
#'
#' @param residue_data Data frame or tibble containing non-A residue
#'   predictions from the ninetails pipeline.
#'
#' @param grouping_factor Character string or \code{NULL} (default). A
#'   grouping variable (e.g. \code{"sample_name"}, \code{"group"}).
#'
#' @param transcript_column Character string. Name of the column
#'   containing transcript identifiers (e.g.
#'   \code{"ensembl_transcript_id_short"}).
#'
#' @param ref Character string, character vector, or \code{NULL}
#'   (default). Whitelist of transcripts with hybrid tails. Built-in
#'   options:
#' \describe{
#'   \item{\code{"athaliana"}}{\emph{Arabidopsis thaliana}}
#'   \item{\code{"hsapiens"}}{\emph{Homo sapiens}}
#'   \item{\code{"mmusculus"}}{\emph{Mus musculus}}
#'   \item{\code{"scerevisiae"}}{\emph{Saccharomyces cerevisiae}}
#'   \item{\code{"celegans"}}{\emph{Caenorhabditis elegans}}
#'   \item{\code{"tbrucei"}}{\emph{Trypanosoma brucei}}
#' }
#'   A custom character vector of transcript IDs may also be provided.
#'   Must be consistent with the content of \code{transcript_column}.
#'
#' @return A tibble based on \code{residue_data} with the following
#'   additional columns:
#' \describe{
#'   \item{mode_pos}{Integer. Most frequent non-A position reported for
#'     the transcript.}
#'   \item{mode_len}{Integer. Most frequent tail length reported for the
#'     transcript.}
#'   \item{seg_err_quart}{Numeric. 0.05 quantile of tail length for the
#'     transcript.}
#'   \item{qc_pos}{Character. Quality flag: \code{"Y"} for likely
#'     genuine, \code{"N"} for ambiguous.}
#'   \item{pos_err_quart}{Numeric. 0.05 quantile of non-A position for
#'     the transcript and prediction type.}
#'   \item{count_nonA}{Integer. Number of non-A-containing reads for
#'     the transcript.}
#'   \item{count}{Integer. Total number of reads for the transcript.}
#' }
#'
#' @seealso \code{\link{correct_class_data}} for the companion function
#'   that reclassifies reads based on this output,
#'   \code{\link{reclassify_ninetails_data}} for the high-level wrapper,
#'   \code{\link{check_tails_guppy}} and \code{\link{create_outputs}} for
#'   the pipeline that produces the input data.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' residue_data_edited <- ninetails::correct_residue_data(
#'   class_data = results[[1]],
#'   residue_data = results[[2]],
#'   transcript_column = "contig")
#'
#' }
#'
correct_residue_data <- function(class_data,
                                 residue_data,
                                 grouping_factor = NULL,
                                 transcript_column,
                                 ref = NULL) {

  # assertions
  if (missing(class_data)) {
    stop(
      "The class_data argument is missing. Please provide the valid class prediction dataframe.",
      call. = FALSE
    )
  }

  if (missing(residue_data)) {
    stop(
      "The residue_data argument is missing. Please provide the valid residue prediction dataframe.",
      call. = FALSE
    )
  }

  if (missing(transcript_column)) {
    stop(
      "The transcript_column argument is missing. Please provide the name of the column that stores the transcript IDs.",
      call. = FALSE
    )
  }

  if (!is.data.frame(class_data) || nrow(class_data) == 0) {
    stop(
      "Empty data frame provided as an input (class_data). Please provide valid input"
    )
  }

  if (!is.data.frame(residue_data) || nrow(residue_data) == 0) {
    stop(
      "Empty data frame provided as an input (residue_data). Please provide valid input"
    )
  }

  #prevent bugs
  class_data <- unique(class_data)
  residue_data <- unique(residue_data)

  if (is.null(grouping_factor)) {
    # mark positions located in first quantile of length (i.e. in close proximity to the transcript body;
    # those are ambiguous as they can arise as the segmentation artifacts inherited from nanopolish)
    class_data_filtered <- class_data %>%
      dplyr::group_by(!!rlang::sym(transcript_column)) %>%
      dplyr::mutate(
        seg_err_quart = stats::quantile(polya_length, probs = 0.05),
        mode_len = which.max(tabulate(polya_length)),
        count = dplyr::n()
      ) %>%
      dplyr::ungroup()

    # mark most frequent position of nonA residue; the mode will be then used to filter out positions
    # which are most likely artifacts
    residue_data_filtered <- residue_data %>%
      dplyr::group_by(!!rlang::sym(transcript_column), prediction) %>%
      dplyr::mutate(
        mode_pos = which.max(tabulate(est_nonA_pos)),
        pos_err_quart = stats::quantile(est_nonA_pos, probs = 0.05),
        count_nonA = dplyr::n()
      ) %>%
      dplyr::ungroup()
  } else {
    # add new columns to class data
    class_data_filtered <- class_data %>%
      dplyr::group_by(
        !!rlang::sym(grouping_factor),
        !!rlang::sym(transcript_column)
      ) %>%
      dplyr::mutate(
        seg_err_quart = stats::quantile(polya_length, probs = 0.05),
        mode_len = which.max(tabulate(polya_length)),
        count = dplyr::n()
      ) %>%
      dplyr::ungroup()

    # add new columns to residue data
    residue_data_filtered <- residue_data %>%
      dplyr::group_by(
        !!rlang::sym(grouping_factor),
        !!rlang::sym(transcript_column),
        prediction
      ) %>%
      dplyr::mutate(
        mode_pos = which.max(tabulate(est_nonA_pos)),
        pos_err_quart = stats::quantile(est_nonA_pos, probs = 0.05),
        count_nonA = dplyr::n()
      ) %>%
      dplyr::ungroup()
  }

  # load whitelists
  path_to_builtin_whitelists <- system.file(
    "extdata",
    "whitelists",
    "whitelist.RData",
    package = "ninetails"
  )
  load(path_to_builtin_whitelists)

  #whitelists
  if (ref == "mmusculus") {
    whitelist = mouse_whitelist
  } else if (ref == "hsapiens") {
    whitelist = human_whitelist
  } else if (ref == "scerevisiae") {
    whitelist = saccer_whitelist
  } else if (ref == "celegans") {
    whitelist = celegans_whitelist
  } else if (ref == "athaliana") {
    whitelist = arabidopsis_whitelist
  } else if (ref == "tbrucei") {
    whitelist = trypa_whitelist
  } else {
    whitelist = ref
  }

  # merge the filtered data
  residue_data_edited <- residue_data_filtered %>%
    dplyr::left_join(class_data_filtered)
  # mark ambiguous positions
  residue_data_edited <- residue_data_edited %>%
    dplyr::mutate(
      qc_pos = dplyr::case_when(
        !!rlang::sym(transcript_column) %in% whitelist ~ "Y",
        count_nonA > 10 &
          (est_nonA_pos > mode_pos & mode_pos < pos_err_quart) ~ "Y",
        count_nonA > 10 &
          (est_nonA_pos > mode_pos & mode_pos > pos_err_quart) &
          est_nonA_pos - pos_err_quart > 4 ~ "Y",
        count_nonA > 10 & mode_pos > seg_err_quart ~ "Y",
        count_nonA > 10 &
          est_nonA_pos > pos_err_quart &
          est_nonA_pos > seg_err_quart ~ "Y",
        count_nonA < 10 & est_nonA_pos > seg_err_quart ~ "Y",
        count_nonA > 10 &
          polya_length < 50 &
          (est_nonA_pos / polya_length) * 100 < mode_len &
          est_nonA_pos > pos_err_quart &
          pos_err_quart > 10 ~ "Y",
        count < 10 & count_nonA < 10 ~ "Y",
        count_nonA < 10 ~ "Y",
        TRUE ~ "N"
      )
    )

  return(residue_data_edited)
}

#' Corrects the classification of reads contained in the class_data table.
#'
#' Introduces two additional columns (\code{corr_class} and
#' \code{corr_comments}) to the \code{class_data} table, reclassifying
#' reads based on the positional quality flags computed by
#' \code{\link{correct_residue_data}}. Reads whose \emph{all} non-A
#' positions were flagged as ambiguous (\code{qc_pos == "N"}) are
#' downgraded from \code{"decorated"} to \code{"blank"}.
#'
#' @details
#' The reclassification logic: for each read, if every non-A residue
#' position has \code{qc_pos == "N"}, the read is reclassified as
#' \code{"blank"} with comment \code{"MPU"} (to maintain compatibility
#' with the tag system used in plotting functions). Reads with at least
#' one \code{qc_pos == "Y"} position retain their original class and
#' comment.
#'
#' @section Caution:
#' It is recommended to use \code{\link{reclassify_ninetails_data}}
#' for downstream analyses, as it wraps this function together with
#' \code{\link{correct_residue_data}} and performs the necessary column
#' renaming. If using this function directly, rename \code{corr_class}
#' and \code{corr_comments} to \code{class} and \code{comments} before
#' plotting.
#'
#' @param residue_data_edited Data frame or tibble. Corrected non-A
#'   residue predictions produced by \code{\link{correct_residue_data}}.
#'   Must contain columns \code{mode_pos}, \code{seg_err_quart}, and
#'   \code{qc_pos}.
#'
#' @param class_data Data frame or tibble containing read_classes
#'   predictions from the ninetails pipeline.
#'
#' @return A tibble identical to \code{class_data} with two additional
#'   columns:
#' \describe{
#'   \item{corr_class}{Character. Corrected classification: original
#'     \code{class} value, or \code{"blank"} for artefact-only reads.}
#'   \item{corr_comments}{Character. Corrected comment: original
#'     \code{comments} value, or \code{"MPU"} for reclassified reads.}
#' }
#'
#' @seealso \code{\link{correct_residue_data}} for the preceding step,
#'   \code{\link{reclassify_ninetails_data}} for the high-level wrapper,
#'   \code{\link{create_outputs}} for the pipeline that produces the
#'   input data.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' class_data_corrected <- ninetails::correct_class_data(
#'   residue_data_edited = residue_data_edited,
#'   class_data = results[[1]])
#'
#' }
#'
correct_class_data <- function(residue_data_edited, class_data) {

  # assertions
  if (missing(residue_data_edited)) {
    stop(
      "The residue_data_edited argument is missing. Please provide the output of correct_residue_data function.",
      call. = FALSE
    )
  }

  if (missing(class_data)) {
    stop(
      "The class_data argument is missing. Please provide the valid class prediction dataframe.",
      call. = FALSE
    )
  }

  if (!is.data.frame(residue_data_edited) || nrow(residue_data_edited) == 0) {
    stop(
      "Empty data frame provided as an input (residue_data_edited). Please provide valid input"
    )
  }

  if (!is.data.frame(class_data) || nrow(class_data) == 0) {
    stop(
      "Empty data frame provided as an input (class_data). Please provide valid input"
    )
  }

  basic_colnames = c("mode_pos", "seg_err_quart", "qc_pos")

  assert_condition(
    basic_colnames[1] %in% colnames(residue_data_edited),
    "mode_pos column is missing in the input residue_data_edited. Is that valid output of correct_residue_data()?"
  )
  assert_condition(
    basic_colnames[2] %in% colnames(residue_data_edited),
    "seg_err_quart column is missing in the input residue_data_edited. Is that valid output of correct_residue_data()?"
  )
  assert_condition(
    basic_colnames[3] %in% colnames(residue_data_edited),
    "qc_pos column is missing in the input residue_data_edited. Is that valid output of correct_residue_data()?"
  )

  #prevent bugs
  class_data <- unique(class_data)
  residue_data_edited <- unique(residue_data_edited)

  # prepare residue summary
  residue_data_summarized <- residue_data_edited %>%
    dplyr::group_by(readname) %>%
    dplyr::summarize(n_resid = dplyr::n(), no_qc_pos_N = sum(qc_pos == "N"))

  # deal with class data
  class_data <- class_data %>%
    dplyr::left_join(residue_data_summarized, by = "readname") %>%
    dplyr::mutate(
      corr_class = dplyr::case_when(
        n_resid == no_qc_pos_N ~ "blank",
        n_resid > no_qc_pos_N ~ "decorated",
        TRUE ~ class
      ),
      corr_comments = dplyr::case_when(
        class == corr_class ~ comments,
        TRUE ~ "MPU"
      )
    ) %>%
    dplyr::select(-c(n_resid, no_qc_pos_N))

  return(class_data)
}

#' Reclassifies ambiguous non-A residues to mitigate potential errors
#' inherited from nanopolish segmentation.
#'
#' High-level wrapper that combines
#' \code{\link{correct_residue_data}} and
#' \code{\link{correct_class_data}} into a single call, producing
#' cleaned class and residue data frames ready for downstream analysis
#' and visualisation.
#'
#' @details
#' Nanopolish segmentation can misidentify nucleotides from A-rich
#' 3' UTR regions as part of the poly(A) tail. When tail boundaries are
#' recognised incorrectly, non-A positions accumulate near the 3' end of
#' the transcript, significantly affecting analysis results. This function
#' flags and removes those ambiguous positions and reclassifies affected
#' reads accordingly.
#'
#' The procedure:
#' \enumerate{
#'   \item \code{\link{correct_residue_data}} annotates each non-A
#'     position with a quality flag (\code{qc_pos}).
#'   \item \code{\link{correct_class_data}} reclassifies reads whose
#'     \emph{all} non-A positions are flagged as ambiguous.
#'   \item Ambiguous positions (\code{qc_pos == "N"}) are dropped from
#'     the residue table.
#'   \item Corrected columns (\code{corr_class}, \code{corr_comments})
#'     are renamed to \code{class} and \code{comments}.
#' }
#'
#' @section Caution:
#' Reads containing only ambiguous non-A positions are reclassified as
#' \code{"blank"} in the \code{class} column, and their \code{comments}
#' are changed from \code{"YAY"} to \code{"MPU"}.
#'
#' @param residue_data Data frame or tibble containing non-A residue
#'   predictions from the ninetails pipeline.
#'
#' @param class_data Data frame or tibble containing read_classes
#'   predictions from the ninetails pipeline.
#'
#' @param grouping_factor Character string or \code{NULL} (default). A
#'   grouping variable (e.g. \code{"sample_name"}).
#'
#' @param transcript_column Character string. Name of the column
#'   containing transcript identifiers (e.g. \code{"contig"},
#'   \code{"ensembl_transcript_id_short"}).
#'
#' @param ref Character string, character vector, or \code{NULL}
#'   (default). Whitelist of transcripts with hybrid tails. Built-in
#'   options:
#' \describe{
#'   \item{\code{"athaliana"}}{\emph{Arabidopsis thaliana}}
#'   \item{\code{"hsapiens"}}{\emph{Homo sapiens}}
#'   \item{\code{"mmusculus"}}{\emph{Mus musculus}}
#'   \item{\code{"scerevisiae"}}{\emph{Saccharomyces cerevisiae}}
#'   \item{\code{"celegans"}}{\emph{Caenorhabditis elegans}}
#'   \item{\code{"tbrucei"}}{\emph{Trypanosoma brucei}}
#' }
#'   A custom character vector may also be provided. Must be consistent
#'   with the content of \code{transcript_column}. Using a whitelist is
#'   optional but allows retrieval of more true positive data.
#'
#' @return A named list with two data frames:
#' \describe{
#'   \item{class_data}{Data frame. Corrected read classifications with
#'     \code{class} and \code{comments} columns updated. Compatible with
#'     plotting functions.}
#'   \item{residue_data}{Data frame. Filtered non-A residue predictions
#'     with ambiguous positions removed. Intermediate QC columns are
#'     dropped.}
#' }
#'
#' @seealso \code{\link{correct_residue_data}} and
#'   \code{\link{correct_class_data}} for the underlying steps,
#'   \code{\link{check_tails_guppy}} and \code{\link{create_outputs}} for
#'   the pipeline that produces the input data.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' rec_results <- ninetails::reclassify_ninetails_data(
#'   residue_data = results[[2]],
#'   class_data = results[[1]],
#'   transcript_column = "contig")
#'
#' }
#'
reclassify_ninetails_data <- function(residue_data,
                                      class_data,
                                      grouping_factor = NULL,
                                      transcript_column,
                                      ref = NULL) {

  # assertions
  if (missing(residue_data)) {
    stop(
      "The residue_data argument is missing. Please provide the valid residue_data (output of core ninetails pipeline).",
      call. = FALSE
    )
  }

  if (missing(class_data)) {
    stop(
      "The class_data argument is missing. Please provide the valid class_data (output of core ninetails pipeline).",
      call. = FALSE
    )
  }

  if (missing(transcript_column)) {
    stop(
      "The transcript_column argument is missing. Please provide the name of the column that stores the transcript IDs.",
      call. = FALSE
    )
  }

  if (!is.data.frame(class_data) || nrow(class_data) == 0) {
    stop(
      "Empty data frame provided as an input (class_data). Please provide valid input"
    )
  }

  if (!is.data.frame(residue_data) || nrow(residue_data) == 0) {
    stop(
      "Empty data frame provided as an input (residue_data). Please provide valid input"
    )
  }

  #prevent bugs
  class_data <- unique(class_data)
  residue_data <- unique(residue_data)

  # load whitelists
  path_to_builtin_whitelists <- system.file(
    "extdata",
    "whitelists",
    "whitelist.RData",
    package = "ninetails"
  )
  load(path_to_builtin_whitelists)

  # deal with ambiguous positions
  residue_data_edited <- ninetails::correct_residue_data(
    class_data = class_data,
    residue_data = residue_data,
    grouping_factor = grouping_factor,
    transcript_column = transcript_column,
    ref = ref
  )

  # deal with classification
  class_data <- ninetails::correct_class_data(
    residue_data_edited = residue_data_edited,
    class_data = class_data
  ) %>%
    dplyr::select(-c(class, comments)) %>%
    dplyr::rename(class = corr_class, comments = corr_comments)

  # filter residue data - drop ambiguous positions
  residue_data_edited <- residue_data_edited %>%
    dplyr::filter(qc_pos == "Y") %>%
    dplyr::select(
      -c(
        mode_pos,
        seg_err_quart,
        qc_pos,
        pos_err_quart,
        count_nonA,
        mode_len,
        count
      )
    )

  #reassure that tibbles are coerced to df
  class_data <- as.data.frame(class_data)
  residue_data_edited <- as.data.frame(residue_data_edited)

  # produce the output
  output <- list()
  output[["class_data"]] <- class_data
  output[["residue_data"]] <- residue_data_edited

  return(output)
}

#' Counts reads by number of non-A occurrence instances.
#'
#' Categorises reads into three abundance bins based on how many
#' separate non-A residue occurrences (instances) were detected per
#' read: \code{"single"} (1), \code{"two"} (2), or \code{"more"} (>= 3).
#' Returns the count of reads in each bin, optionally stratified by a
#' grouping variable.
#'
#' @param residue_data Data frame or tibble containing non-A residue
#'   predictions from the ninetails pipeline.
#'
#' @param grouping_factor Character string (default \code{NA}). Name of a
#'   column to use as a grouping variable (e.g. \code{"sample_name"}).
#'
#' @return A tibble with columns for the grouping variable (if provided),
#'   \code{instances} (one of \code{"single"}, \code{"two"},
#'   \code{"more"}), and \code{count} (the number of reads in each bin).
#'
#' @seealso \code{\link{count_residues}} for counting total residue hits,
#'   \code{\link{count_class}} for counting read-level classes,
#'   \code{\link{read_residue_single}} for loading residue data.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' nonA_abundance <- ninetails::count_nonA_abundance(
#'   residue_data = residue_data,
#'   grouping_factor = "sample_name")
#'
#' }
count_nonA_abundance <- function(residue_data, grouping_factor = NA) {

  # assertions
  if (!is.data.frame(residue_data) || nrow(residue_data) == 0) {
    stop(
      "Empty data frame provided as an input (residue_data). Please provide valid input"
    )
  }

  if (!is.na(grouping_factor)) {
    assert_condition(
      grouping_factor %in% colnames(residue_data),
      paste0(grouping_factor, " is not a column of input dataset")
    )
  }

  nonA_counts <- residue_data %>%
    dplyr::select(!!rlang::sym(grouping_factor), readname) %>%
    dplyr::group_by(!!rlang::sym(grouping_factor), readname) %>%
    dplyr::mutate(
      instances = ifelse(
        dplyr::n() == 1,
        "single",
        ifelse(dplyr::n() == 2, "two", "more")
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(!!rlang::sym(grouping_factor), instances) %>%
    dplyr::mutate(count = dplyr::n()) %>%
    dplyr::select(-readname) %>%
    dplyr::distinct() %>%
    dplyr::ungroup()

  return(nonA_counts)
}
