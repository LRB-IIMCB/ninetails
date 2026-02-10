#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom utils adist
## usethis namespace: end


# NSE variable declarations for R CMD check
# These variables are used in dplyr/tidyr/ggplot2/foreach NSE contexts
# throughout the package and would otherwise generate
# "no visible binding for global variable" warnings.
#
# This is the SINGLE, CENTRALISED place for all such declarations.
# Do NOT use `variable <- NULL` inside individual functions for this purpose.

utils::globalVariables(c(

  # -- Read and sequence identifiers --
  "readname", "read_id", "chunkname", "filename",

  # -- Poly(A) tail characteristics --
  "polya_length", "polya_start", "transcript_start", "qc_tag",
  "tail_start", "tail_end", "tail_length",  # tailfindr compatibility
  "poly_tail_start", "poly_tail_length",     # dorado-specific

  # -- Non-A position and prediction data --
  "est_nonA_pos", "est_nonA_pos_2", "prediction", "class", "comments",
  "prediction_C", "prediction_G", "prediction_U",

  # -- Alignment information --
  "alignment_direction", "alignment_mapq", "alignment_genome",

  # -- Sample and file paths --
  "sample_name", "class_path", "residue_path",

  # -- Data contents from file reading --
  "class_contents", "residue_contents",

  # -- Summary and aggregation variables --
  "nonA_residues", "n", "prop", "total",
  "sum_nonA", "counts_nonA", "counts_blank", "counts_total", "count_nonA",

  # -- Quality control and correction variables --
  "qc_pos", "n_resid", "no_qc_pos_N",
  "corr_class", "corr_comments",
  "mode_pos", "mode_len",
  "seg_err_quart", "pos_err_quart",

  # -- BiomaRt annotation --
  "ensembl_transcript_id", "ensembl_transcript_id_short",

  # -- Species-specific whitelists --
  "mouse_whitelist", "human_whitelist", "saccer_whitelist",
  "celegans_whitelist", "arabidopsis_whitelist", "trypa_whitelist",

  # -- Signal processing variables --
  "position", "time", "signal", "pA", "segment", "moves",
  "sampling_rate", "offset", "range", "digitisation",
  "adapter_start", "leader_start",

  # -- Residue types and counts --
  "C", "G", "U",
  "pC", "pG", "pU",
  "instances", "count", "counts",
  "counts_C", "counts_G", "counts_U",
  "two", "more", "single",
  "psingle", "ptwo", "pmore",

  # -- Grouping and categorical variables --
  "group", "contig",

  # -- Plotting panel: distribution / density / position --
  "ygreki", "type", "label", "cat",
  "binned_lengths", "binned_positions",
  "polya_median", "polya_mean",
  "median_value", "mean_value", "mode_value",
  "counts_of_reads_equal_or_longer_than_est_position",
  "counts_of_reads_with_nonA_in_given_position",
  "center_value", "xintercept", "variable_to_plot",

  # -- Generic plotting variables --
  "Var1", "Var2", "value", "nam",
  "..density..", "..ndensity..",

  # -- Statistics variables --
  "transcript_id", "significance",
  "padj", "data", "stats",

  # -- Iterator / pipe placeholders --
  "i", "x", "."
))
