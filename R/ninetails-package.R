#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom utils adist
## usethis namespace: end


# NSE variable declarations for R CMD check
# These variables are used in dplyr/tidyr/ggplot2 NSE contexts throughout the package
# and would otherwise generate "no visible binding for global variable" warnings

utils::globalVariables(c(
  # Read and sequence identifiers
  "readname", "read_id", "chunkname", "filename",

  # Poly(A) tail characteristics
  "polya_length", "polya_start", "transcript_start", "qc_tag",
  "tail_start", "tail_end", "tail_length",  # tailfindr compatibility
  "poly_tail_start", "poly_tail_length",     # dorado-specific

  # Non-A position and prediction data
  "est_nonA_pos", "prediction", "class", "comments",
  "prediction_C", "prediction_G", "prediction_U",  # residue-specific predictions

  # Alignment information
  "alignment_direction", "alignment_mapq", "alignment_genome",

  # Sample and file paths
  "sample_name", "class_path", "residue_path",

  # Data contents from file reading
  "class_contents", "residue_contents",

  # Summary data
  "nonA_residues", "n", "prop", "total",
  "sum_nonA", "counts_nonA", "counts_blank", "counts_total", "count_nonA",

  # Quality control and correction variables
  "qc_pos", "n_resid", "no_qc_pos_N",
  "corr_class", "corr_comments",
  "mode_pos", "mode_len",
  "seg_err_quart", "pos_err_quart",

  # BiomaRt annotation
  "ensembl_transcript_id", "ensembl_transcript_id_short",

  # Species-specific whitelists
  "mouse_whitelist", "human_whitelist", "saccer_whitelist",
  "celegans_whitelist", "arabidopsis_whitelist", "trypa_whitelist",

  # Signal processing variables
  "position", "time", "signal", "pA", "segment", "moves",
  "sampling_rate", "offset", "range", "digitisation",
  "adapter_start", "leader_start",  # fast5 plotting

  # Residue types and counts
  "C", "G", "U",           # nucleotide residue columns
  "pC", "pG", "pU",        # proportions
  "instances", "count",     # counting variables
  "two", "more",           # abundance categories
  "psingle", "ptwo", "pmore",  # abundance proportions

  # Grouping and categorical variables
  "group", "contig",

  # Generic plotting variables
  "Var1", "Var2", "value", "nam",
  "..density..", "..ndensity..",  # ggplot2 computed aesthetics

  # Iterator variables from loops and functional programming
  "i", "x", ".",  # foreach, purrr::map, magrittr pipe

  # Statistical variables (note: some clash with base R functions)
  "padj", "data", "stats"  # data and stats also exist as functions/packages
))
