#' Annotates the ninetails output data
#'
#' This function is a wrapper for data annotation using the biomaRt package.
#'
#' For proper execution, this function requires the biomaRt package version
#' greater than 2.40.
#'
#' One of the arguments 'organism' or 'mart_to_use' must be defined.
#' The arguments 'organism' and 'mart_to_use' are mutually exclusive.
#' If both are declared, the function will throw an error.
#'
#' Function based on PK (smaegol) NanoTail annotate_with_biomart()
#' https://github.com/LRB-IIMCB/nanotail/
#' Many thanks to the NanoTail developer for help and kind advice!
#'
#' @param input_data [dataframe] tabular output of ninetails pipeline.
#'
#' @param attributes_to_get [character] annotation attributes to be retrieved.
#' By default: 'external_gene_name','description','transcript_biotype'.
#'
#' @param filters column of the input dataframe to be matched
#' with the target mart (e.g. "ensembl_transcript_id_short")
#'
#' @param organism [character] character string. Currently available::\itemize{
#' \item 'athaliana' - Arabidopsis thaliana
#' \item 'hsapiens' - Homo sapiens
#' \item 'mmusculus' - Mus musculus
#' \item 'scerevisiae' - Saccharomyces cerevisiae
#' }
#' This argument is optional & mutually exclusive with 'mart_to_use'.
#'
#' @param mart_to_use optional: mart object created with \link[biomaRt]{useMart}
#' or \link[biomaRt]{useEnsembl} function. This argument is optional
#' & mutually exclusive with 'organism'.
#'
#' @return A dataframe with annotated ninetails output data.
#'
#' @export
#'
#' @examples
#'  \dontrun{
#' # with provided organism
#' annot <- ninetails::annotate_with_biomart(input_data=merged_nonA_tables,
#'                                           organism="mmusculus")
#'
#' # with provided mart (where 'mart' is the output of \link[biomaRt]{useMart} function)
#' annot <- ninetails::annotate_with_biomart(input_data=merged_nonA_tables,
#'                                           mart_to_use=mart)
#' }
#'
annotate_with_biomart <- function(input_data,
                                  attributes_to_get=c('ensembl_transcript_id',
                                                      'external_gene_name',
                                                      'description',
                                                      'transcript_biotype'),
                                  filters='ensembl_transcript_id',
                                  organism=NULL,
                                  mart_to_use=NULL){

  # variable binding
  ensembl_transcript_id_short <- ensembl_transcript_id <- NULL

  #assertions
  if (missing(input_data)) {
    stop("Please provide a dataframe with ninetails results as an input.",
         call. = FALSE)
  }


  if (!is.data.frame(input_data) || nrow(input_data) == 0) {
    stop("Empty data frame provided as an input (input_data). Please provide valid input")}


  assertthat::assert_that(length(attributes)>0,
                          msg="please provide attributes to get from biomart")

  # Assert that at least one of organism or mart_to_use is provided
  if (!is.null(organism) && !is.null(mart_to_use)) {
    stop("Only one of 'organism' or 'mart_to_use' can be provided.")
  }
  # Check that at least one of organism and mart_to_use is provided
  if (is.null(organism) && is.null(mart_to_use)) {
    stop("Either 'organism' or 'mart_to_use' must be provided.")
  }

  # Check that mart_to_use is a mart object if it is provided
  if (!is.null(mart_to_use) && !class(mart_to_use)=="Mart") {
    stop("The 'mart_to_use' must be a mart. Please provide valid object.")
  }

  # Check the version of the biomaRt package
  if (utils::packageVersion("biomaRt") < "2.40.0") {
    stop("This function requires a biomaRt package version greater than 2.40.")
  }

  ensembl_ids <- unique(input_data$ensembl_transcript_id_short)
  ensembl_ids <- ensembl_ids[!is.na(ensembl_ids)]

  input_data <- input_data %>% dplyr::rename(ensembl_transcript_id = ensembl_transcript_id_short)


  if (is.null(organism)){

    annotation_data <- biomaRt::getBM(attributes=attributes_to_get,
                                      filters = filters,
                                      values = ensembl_ids,
                                      mart = mart_to_use)

  } else if (organism=="mmusculus") {

    mart_to_use <- biomaRt::useMart(biomart='ENSEMBL_MART_ENSEMBL',
                                    dataset = 'mmusculus_gene_ensembl')


  } else if(organism=="hsapiens") {

    mart_to_use <- biomaRt::useMart(biomart='ENSEMBL_MART_ENSEMBL',
                                    dataset = 'hsapiens_gene_ensembl')

  } else if (organism=="athaliana") {

    mart_to_use <- biomaRt::useEnsemblGenomes(biomart = "plants_mart",
                                              dataset = "athaliana_eg_gene")

  } else if (organism=="scerevisiae"){

    mart_to_use <- biomaRt::useEnsemblGenomes(biomart = "fungi_mart",
                                              dataset = "scerevisiae_eg_gene")

  } else {
    stop("Please provide valid organism name (mmusculus, hsapiens, athaliana) or mart object (mart_to_use)",
         call. = FALSE)
  }


  annotation_data <- biomaRt::getBM(attributes=attributes_to_get,
                                    filters = filters,
                                    values = ensembl_ids,
                                    mart = mart_to_use)


  if (!is.data.frame(annotation_data) || nrow(annotation_data) == 0) {
    stop("Could not retrieve annotation data for given transcript IDs. Check if the organism or mart_to_use was defined correctly")}


  input_data_annotated <-  input_data %>% dplyr::left_join(annotation_data) %>%
    dplyr::rename(ensembl_transcript_id_short = ensembl_transcript_id)

  return(input_data_annotated)
}
