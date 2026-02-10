################################################################################
# ANNOTATE WITH BIOMART
################################################################################
#' Annotate ninetails output data with biomaRt
#'
#' Retrieves gene-level annotation from Ensembl for transcripts in ninetails
#' output data. This is a convenience wrapper for \code{\link[biomaRt]{getBM}}
#' with built-in organism presets and support for custom mart objects.
#'
#' Requires the \pkg{biomaRt} package version >= 2.40.
#'
#' Exactly one of \code{organism} or \code{mart_to_use} must be provided.
#' The two arguments are mutually exclusive; if both are declared the function
#' throws an error.
#'
#' @section Acknowledgements:
#' Based on PK (smaegol) NanoTail \code{annotate_with_biomart()}:
#' \url{https://github.com/LRB-IIMCB/nanotail/}.
#' Many thanks to the NanoTail developer for help and kind advice.
#'
#' @param input_data Data frame. Tabular output of the ninetails pipeline.
#'   Must contain a column named \code{ensembl_transcript_id_short}.
#'
#' @param attributes_to_get Character vector. Annotation attributes to
#'   retrieve from biomaRt. Default:
#'   \code{c("ensembl_transcript_id", "external_gene_name", "description",
#'   "transcript_biotype")}.
#'
#' @param filters Character string. Column of the input data frame to match
#'   with the target mart. Default: \code{"ensembl_transcript_id"}.
#'
#' @param organism Character string or \code{NULL}. Organism shorthand for
#'   built-in mart presets. Currently available:
#'   \describe{
#'     \item{\code{"athaliana"}}{\emph{Arabidopsis thaliana} (Ensembl Plants)}
#'     \item{\code{"hsapiens"}}{\emph{Homo sapiens} (Ensembl)}
#'     \item{\code{"mmusculus"}}{\emph{Mus musculus} (Ensembl)}
#'     \item{\code{"scerevisiae"}}{\emph{Saccharomyces cerevisiae}
#'       (Ensembl Fungi)}
#'   }
#'   Mutually exclusive with \code{mart_to_use}. Default: \code{NULL}.
#'
#' @param mart_to_use Mart object or \code{NULL}. A mart object created with
#'   \code{\link[biomaRt]{useMart}} or \code{\link[biomaRt]{useEnsembl}}.
#'   Mutually exclusive with \code{organism}. Default: \code{NULL}.
#'
#' @return A data frame with the original ninetails output data joined with
#'   the retrieved annotation attributes via left join on transcript IDs.
#'
#' @seealso
#' \code{\link{merge_nonA_tables}} for preparing the input,
#' \code{\link[biomaRt]{getBM}} for the underlying biomaRt query,
#' \code{\link[biomaRt]{useMart}} for creating custom mart objects
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # With built-in organism preset
#' annot <- ninetails::annotate_with_biomart(
#'   input_data = merged_nonA_tables,
#'   organism = "mmusculus"
#' )
#'
#' # With custom mart object
#' mart <- biomaRt::useMart(
#'   biomart = "ENSEMBL_MART_ENSEMBL",
#'   dataset = "hsapiens_gene_ensembl"
#' )
#' annot <- ninetails::annotate_with_biomart(
#'   input_data = merged_nonA_tables,
#'   mart_to_use = mart
#' )
#'
#' }
#'
annotate_with_biomart <- function(input_data,
                                  attributes_to_get = c('ensembl_transcript_id',
                                                        'external_gene_name',
                                                        'description',
                                                        'transcript_biotype'),
                                  filters = 'ensembl_transcript_id',
                                  organism = NULL,
                                  mart_to_use = NULL) {

  #assertions
  if (missing(input_data)) {
    stop(
      "Please provide a dataframe with ninetails results as an input.",
      call. = FALSE
    )
  }

  if (!is.data.frame(input_data) || nrow(input_data) == 0) {
    stop(
      "Empty data frame provided as an input (input_data). Please provide valid input"
    )
  }

  assert_condition(
    length(attributes) > 0,
    "please provide attributes to get from biomart"
  )

  # Assert that at least one of organism or mart_to_use is provided
  if (!is.null(organism) && !is.null(mart_to_use)) {
    stop("Only one of 'organism' or 'mart_to_use' can be provided.")
  }

  # Check that at least one of organism and mart_to_use is provided
  if (is.null(organism) && is.null(mart_to_use)) {
    stop("Either 'organism' or 'mart_to_use' must be provided.")
  }

  # Check that mart_to_use is a mart object if it is provided
  if (!is.null(mart_to_use) && !inherits(mart_to_use, "Mart")) {
    stop("The 'mart_to_use' must be a mart. Please provide valid object.")
  }

  # Check the version of the biomaRt package
  if (utils::packageVersion("biomaRt") < "2.40.0") {
    stop("This function requires a biomaRt package version greater than 2.40.")
  }

  ensembl_ids <- unique(input_data$ensembl_transcript_id_short)
  ensembl_ids <- ensembl_ids[!is.na(ensembl_ids)]

  input_data <- input_data %>%
    dplyr::rename(ensembl_transcript_id = ensembl_transcript_id_short)

  if (is.null(organism)) {
    annotation_data <- biomaRt::getBM(
      attributes = attributes_to_get,
      filters = filters,
      values = ensembl_ids,
      mart = mart_to_use
    )
  } else if (organism == "mmusculus") {
    mart_to_use <- biomaRt::useMart(
      biomart = 'ENSEMBL_MART_ENSEMBL',
      dataset = 'mmusculus_gene_ensembl'
    )
  } else if (organism == "hsapiens") {
    mart_to_use <- biomaRt::useMart(
      biomart = 'ENSEMBL_MART_ENSEMBL',
      dataset = 'hsapiens_gene_ensembl'
    )
  } else if (organism == "athaliana") {
    mart_to_use <- biomaRt::useEnsemblGenomes(
      biomart = "plants_mart",
      dataset = "athaliana_eg_gene"
    )
  } else if (organism == "scerevisiae") {
    mart_to_use <- biomaRt::useEnsemblGenomes(
      biomart = "fungi_mart",
      dataset = "scerevisiae_eg_gene"
    )
  } else {
    stop(
      "Please provide valid organism name (mmusculus, hsapiens, athaliana) or mart object (mart_to_use)",
      call. = FALSE
    )
  }

  annotation_data <- biomaRt::getBM(attributes = attributes_to_get,
                                    filters = filters,
                                    values = ensembl_ids,
                                    mart = mart_to_use)

  if (!is.data.frame(annotation_data) || nrow(annotation_data) == 0) {
    stop(
      "Could not retrieve annotation data for given transcript IDs. Check if the organism or mart_to_use was defined correctly"
    )
  }

  input_data_annotated <- input_data %>%
    dplyr::left_join(annotation_data) %>%
    dplyr::rename(ensembl_transcript_id_short = ensembl_transcript_id)

  return(input_data_annotated)
}
