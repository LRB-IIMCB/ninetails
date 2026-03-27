################################################################################
# Testing core functions Dorado cDNA
################################################################################

# Helpers
################################################################################

#' Create dummy cli_log function for testing
#' @keywords internal
dummy_cli_log_cdna <- function(msg, level = "INFO", ...) {
  return(invisible(NULL))
}

#' Write a minimal sequence TSV file for detect_orientation tests
#' @keywords internal
write_sequence_tsv <- function(sequences, path) {
  df <- data.frame(
    read_id = paste0("read-", seq_along(sequences)),
    sequence = sequences,
    stringsAsFactors = FALSE
  )
  vroom::vroom_write(df, path, delim = "\t")
  return(invisible(path))
}

# SSP primer used by detect_orientation_single (forward strand indicator)
SSP_PRIMER <- "TTTCTGTTGGTGCTGATATTGCTTT"
# VNP primer (reverse strand indicator)
VNP_PRIMER <- "ACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAGTTTT"


################################################################################
# count_trailing_chars
################################################################################

test_that("count_trailing_chars returns correct count for trailing characters", {
  expect_equal(count_trailing_chars("ACGTTTTT", "T"), 5L)
  expect_equal(count_trailing_chars("ACGT", "A"), 0L)
  expect_equal(count_trailing_chars("TTTT", "T"), 4L)
  expect_equal(count_trailing_chars("A", "A"), 1L)
})

test_that("count_trailing_chars returns 0 when no trailing match", {
  expect_equal(count_trailing_chars("TTTAAA", "T"), 0L)
})

test_that("count_trailing_chars validates input types", {
  expect_error(count_trailing_chars(123, "T"))
  expect_error(count_trailing_chars("ACGT", "TT"))
  expect_error(count_trailing_chars("ACGT", 1))
})


################################################################################
# reverse_complement
################################################################################

test_that("reverse_complement produces correct output", {
  expect_equal(reverse_complement("ATCG"), "CGAT")
  expect_equal(reverse_complement("AAATTT"), "AAATTT")
  expect_equal(reverse_complement("AAAA"), "TTTT")
  expect_equal(reverse_complement("GCGC"), "GCGC")
  expect_equal(reverse_complement("A"), "T")
})

test_that("reverse_complement handles N bases", {
  expect_equal(reverse_complement("ATCGN"), "NCGAT")
})

test_that("reverse_complement handles lowercase by uppercasing first", {
  expect_equal(reverse_complement("atcg"), "CGAT")
})

test_that("reverse_complement validates input", {
  expect_error(reverse_complement(123))
  expect_error(reverse_complement(c("ATCG", "GCTA")))
})


################################################################################
# edit_distance_hw
################################################################################

test_that("edit_distance_hw returns 0 for exact substring match", {
  # Query is a substring of target — sliding window should find distance 0
  expect_equal(edit_distance_hw("ATCG", "XXATCGXX"), 0L)
})

test_that("edit_distance_hw returns 0 for identical strings", {
  expect_equal(edit_distance_hw("ATCG", "ATCG"), 0L)
})

test_that("edit_distance_hw returns positive distance for mismatches", {
  dist <- edit_distance_hw("AAAA", "TTTT")
  expect_true(dist > 0)
})

test_that("edit_distance_hw handles query longer than target", {
  # When query_len > target_len, sliding loop doesn't run;
  # full adist is used and returns a positive value
  dist <- edit_distance_hw("ATCGATCG", "ATCG")
  expect_true(dist >= 0)
})

test_that("edit_distance_hw rejects empty query (is_string requires nchar >= 1)", {
  # is_string enforces nchar >= 1, so "" is not a valid input
  expect_error(edit_distance_hw("", "ATCG"), "query must be a character string")
})

test_that("edit_distance_hw validates input types", {
  expect_error(edit_distance_hw(123, "ATCG"))
  expect_error(edit_distance_hw("ATCG", 123))
})


################################################################################
# detect_orientation_single
################################################################################

test_that("detect_orientation_single returns unknown for short sequences", {
  expect_equal(detect_orientation_single("ACGT"), "unknown")
  expect_equal(detect_orientation_single(""), "unknown")
  expect_equal(detect_orientation_single(NA_character_), "unknown")
})

test_that("detect_orientation_single returns unknown for NA input", {
  expect_equal(detect_orientation_single(NA), "unknown")
})

test_that("detect_orientation_single returns one of three valid values", {
  # Use a 50-character sequence — may be too ambiguous for confident
  # classification, but must return one of the three valid codes.
  result <- detect_orientation_single(paste(rep("A", 100), collapse = ""))
  expect_true(result %in% c("A", "T", "unknown"))
})

test_that("detect_orientation_single classifies forward polyA orientation", {
  # Construct a sequence with the SSP primer at the 5' end and reverse
  # complement of VNP (minus trailing Ts) at the 3' end — the canonical
  # forward (polyA) configuration.
  vnp_trimmed <- sub("T+$", "", VNP_PRIMER)
  vnp_rc <- reverse_complement(vnp_trimmed)
  # Build a sequence: SSP + 200 As + VNP_RC
  seq_fwd <- paste0(SSP_PRIMER,
                    paste(rep("A", 200), collapse = ""),
                    vnp_rc)
  result <- detect_orientation_single(seq_fwd)
  expect_equal(result, "A")
})

test_that("detect_orientation_single classifies reverse polyT orientation", {
  # Reverse configuration: VNP (trimmed) at 5' end and RC of SSP at 3' end
  vnp_trimmed <- sub("T+$", "", VNP_PRIMER)
  ssp_rc <- reverse_complement(SSP_PRIMER)
  seq_rev <- paste0(vnp_trimmed,
                    paste(rep("T", 200), collapse = ""),
                    ssp_rc)
  result <- detect_orientation_single(seq_rev)
  expect_equal(result, "T")
})


################################################################################
# detect_orientation_multiple
################################################################################

test_that("detect_orientation_multiple errors when sequence_files is missing", {
  expect_error(detect_orientation_multiple(),
               "Sequence files are missing")
})

test_that("detect_orientation_multiple errors on non-character sequence_files", {
  expect_error(detect_orientation_multiple(sequence_files = 123),
               "must be a character vector")
})

test_that("detect_orientation_multiple errors when files do not exist", {
  expect_error(
    detect_orientation_multiple(sequence_files = "/nonexistent/file.tsv"),
    "All sequence files must exist"
  )
})

test_that("detect_orientation_multiple errors on non-numeric num_cores", {
  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp), add = TRUE)
  write_sequence_tsv(paste(rep("A", 100), collapse = ""), tmp)

  expect_error(
    detect_orientation_multiple(sequence_files = tmp, num_cores = "one"),
    "positive numeric"
  )
})

test_that("detect_orientation_multiple returns tibble with tail_type column", {
  skip_if_not_installed("vroom")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tibble")

  vnp_trimmed <- sub("T+$", "", VNP_PRIMER)
  vnp_rc  <- reverse_complement(vnp_trimmed)
  ssp_rc  <- reverse_complement(SSP_PRIMER)

  seqs <- c(
    paste0(SSP_PRIMER, paste(rep("A", 200), collapse = ""), vnp_rc),   # polyA
    paste0(vnp_trimmed, paste(rep("T", 200), collapse = ""), ssp_rc),  # polyT
    paste(rep("N", 100), collapse = "")                                 # unidentified
  )

  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp), add = TRUE)
  write_sequence_tsv(seqs, tmp)

  result <- detect_orientation_multiple(sequence_files = tmp,
                                        num_cores = 1,
                                        cli_log = dummy_cli_log_cdna)

  expect_s3_class(result, "tbl_df")
  expect_true("tail_type" %in% colnames(result))
  expect_true(all(result$tail_type %in% c("polyA", "polyT", "unidentified")))
})

test_that("detect_orientation_multiple classifies polyA and polyT reads correctly", {
  skip_if_not_installed("vroom")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tibble")

  vnp_trimmed <- sub("T+$", "", VNP_PRIMER)
  vnp_rc <- reverse_complement(vnp_trimmed)
  ssp_rc <- reverse_complement(SSP_PRIMER)

  seq_polya <- paste0(SSP_PRIMER, paste(rep("A", 200), collapse = ""), vnp_rc)
  seq_polyt <- paste0(vnp_trimmed, paste(rep("T", 200), collapse = ""), ssp_rc)

  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp), add = TRUE)
  write_sequence_tsv(c(seq_polya, seq_polyt), tmp)

  result <- detect_orientation_multiple(sequence_files = tmp,
                                        num_cores = 1,
                                        cli_log = dummy_cli_log_cdna)

  # unname() is needed because sapply in the source preserves the sequence
  # string as the name of each result element
  expect_equal(unname(result$tail_type[1]), "polyA")
  expect_equal(unname(result$tail_type[2]), "polyT")
})

test_that("detect_orientation_multiple handles missing sequence column gracefully", {
  skip_if_not_installed("vroom")
  skip_if_not_installed("tibble")

  # File exists but lacks a 'sequence' column — should skip and return empty tibble
  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp), add = TRUE)
  df <- data.frame(read_id = c("r1", "r2"), other_col = c("x", "y"),
                   stringsAsFactors = FALSE)
  vroom::vroom_write(df, tmp, delim = "\t")

  result <- detect_orientation_multiple(sequence_files = tmp,
                                        num_cores = 1,
                                        cli_log = dummy_cli_log_cdna)

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 0L)
})

test_that("detect_orientation_multiple processes multiple files and combines results", {
  skip_if_not_installed("vroom")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tibble")

  vnp_trimmed <- sub("T+$", "", VNP_PRIMER)
  vnp_rc <- reverse_complement(vnp_trimmed)
  seq_polya <- paste0(SSP_PRIMER, paste(rep("A", 200), collapse = ""), vnp_rc)

  tmp1 <- tempfile(fileext = ".tsv")
  tmp2 <- tempfile(fileext = ".tsv")
  on.exit({ unlink(tmp1); unlink(tmp2) }, add = TRUE)

  write_sequence_tsv(rep(seq_polya, 2), tmp1)
  write_sequence_tsv(rep(seq_polya, 3), tmp2)

  result <- detect_orientation_multiple(sequence_files = c(tmp1, tmp2),
                                        num_cores = 1,
                                        cli_log = dummy_cli_log_cdna)

  expect_equal(nrow(result), 5L)
})


################################################################################
# split_bam_file_cdna — guards testable without Rsamtools or a real BAM
################################################################################

test_that("split_bam_file_cdna errors when Rsamtools is absent", {
  skip_if(requireNamespace("Rsamtools", quietly = TRUE),
          "Rsamtools is installed — skipping unavailability test")

  expect_error(
    split_bam_file_cdna(bam_file = "file.bam",
                        dorado_summary = "summary.txt",
                        save_dir = tempdir(),
                        part_number = 1),
    "Rsamtools.*required"
  )
})

test_that("split_bam_file_cdna errors on missing bam_file argument", {
  skip_if_not_installed("Rsamtools")
  skip_if_not_installed("S4Vectors")

  expect_error(
    split_bam_file_cdna(dorado_summary = "summary.txt",
                        save_dir = tempdir(),
                        part_number = 1),
    "BAM file is missing"
  )
})

test_that("split_bam_file_cdna errors on missing dorado_summary argument", {
  skip_if_not_installed("Rsamtools")
  skip_if_not_installed("S4Vectors")

  expect_error(
    split_bam_file_cdna(bam_file = "file.bam",
                        save_dir = tempdir(),
                        part_number = 1),
    "Summary file is missing"
  )
})

test_that("split_bam_file_cdna errors on missing save_dir argument", {
  skip_if_not_installed("Rsamtools")
  skip_if_not_installed("S4Vectors")

  expect_error(
    split_bam_file_cdna(bam_file = "file.bam",
                        dorado_summary = "summary.txt",
                        part_number = 1),
    "Output directory is missing"
  )
})

test_that("split_bam_file_cdna errors on missing part_number argument", {
  skip_if_not_installed("Rsamtools")
  skip_if_not_installed("S4Vectors")

  expect_error(
    split_bam_file_cdna(bam_file = "file.bam",
                        dorado_summary = "summary.txt",
                        save_dir = tempdir()),
    "Part number is missing"
  )
})

test_that("split_bam_file_cdna errors on non-existent bam_file", {
  skip_if_not_installed("Rsamtools")
  skip_if_not_installed("S4Vectors")

  expect_error(
    split_bam_file_cdna(bam_file = "/nonexistent/file.bam",
                        dorado_summary = "summary.txt",
                        save_dir = tempdir(),
                        part_number = 1),
    "BAM"
  )
})


################################################################################
# extract_data_from_bam — guards testable without Rsamtools or a real BAM
################################################################################

test_that("extract_data_from_bam errors when Rsamtools is absent", {
  skip_if(requireNamespace("Rsamtools", quietly = TRUE),
          "Rsamtools is installed — skipping unavailability test")

  expect_error(
    extract_data_from_bam(bam_file = "file.bam",
                          summary_file = "summary.txt"),
    "Rsamtools.*required"
  )
})

test_that("extract_data_from_bam errors on non-existent bam_file", {
  skip_if_not_installed("Rsamtools")
  skip_if_not_installed("S4Vectors")

  expect_error(
    extract_data_from_bam(bam_file = "/nonexistent/file.bam",
                          summary_file = "summary.txt"),
    "BAM file does not exist"
  )
})

test_that("extract_data_from_bam errors on non-existent summary_file", {
  skip_if_not_installed("Rsamtools")
  skip_if_not_installed("S4Vectors")

  # Create a real (empty) temp file so bam_file check passes, then
  # provide a non-existent summary_file to reach its check.
  tmp_bam <- tempfile(fileext = ".bam")
  file.create(tmp_bam)
  on.exit(unlink(tmp_bam), add = TRUE)

  expect_error(
    extract_data_from_bam(bam_file = tmp_bam,
                          summary_file = "/nonexistent/summary.txt"),
    "Summary file does not exist"
  )
})

test_that("extract_data_from_bam errors when seq_only is not logical", {
  skip_if_not_installed("Rsamtools")
  skip_if_not_installed("S4Vectors")

  tmp_bam <- tempfile(fileext = ".bam")
  tmp_sum <- tempfile(fileext = ".txt")
  file.create(tmp_bam)
  file.create(tmp_sum)
  on.exit({ unlink(tmp_bam); unlink(tmp_sum) }, add = TRUE)

  expect_error(
    extract_data_from_bam(bam_file = tmp_bam,
                          summary_file = tmp_sum,
                          seq_only = "yes"),
    "seq_only.*must be logical"
  )
})


################################################################################
# preprocess_inputs_cdna — validation guards only
################################################################################

test_that("preprocess_inputs_cdna errors on non-numeric num_cores", {
  expect_error(
    preprocess_inputs_cdna(bam_file = "file.bam",
                           dorado_summary = "summary.txt",
                           pod5_dir = tempdir(),
                           num_cores = "one",
                           qc = TRUE,
                           save_dir = tempdir(),
                           prefix = "",
                           part_size = 100,
                           cli_log = dummy_cli_log_cdna),
    "Number of cores must be a positive numeric"
  )
})

test_that("preprocess_inputs_cdna errors on non-character pod5_dir", {
  expect_error(
    preprocess_inputs_cdna(bam_file = "file.bam",
                           dorado_summary = "summary.txt",
                           pod5_dir = 12345,
                           num_cores = 1,
                           qc = TRUE,
                           save_dir = tempdir(),
                           prefix = "",
                           part_size = 100,
                           cli_log = dummy_cli_log_cdna),
    "Pod5 directory path must be a character"
  )
})

test_that("preprocess_inputs_cdna errors on non-existent pod5_dir", {
  expect_error(
    preprocess_inputs_cdna(bam_file = "file.bam",
                           dorado_summary = "summary.txt",
                           pod5_dir = "/nonexistent/pod5",
                           num_cores = 1,
                           qc = TRUE,
                           save_dir = tempdir(),
                           prefix = "",
                           part_size = 100,
                           cli_log = dummy_cli_log_cdna),
    "Pod5.*directory does not exist"
  )
})

test_that("preprocess_inputs_cdna errors on non-logical qc", {
  expect_error(
    preprocess_inputs_cdna(bam_file = "file.bam",
                           dorado_summary = "summary.txt",
                           pod5_dir = tempdir(),
                           num_cores = 1,
                           qc = "yes",
                           save_dir = tempdir(),
                           prefix = "",
                           part_size = 100,
                           cli_log = dummy_cli_log_cdna),
    "qc must be logical"
  )
})

test_that("preprocess_inputs_cdna errors when bam_file is not a string", {
  expect_error(
    preprocess_inputs_cdna(bam_file = 12345,
                           dorado_summary = "summary.txt",
                           pod5_dir = tempdir(),
                           num_cores = 1,
                           qc = TRUE,
                           save_dir = tempdir(),
                           prefix = "",
                           part_size = 100,
                           cli_log = dummy_cli_log_cdna),
    "BAM file must be provided as a character string"
  )
})

test_that("preprocess_inputs_cdna errors when bam_file does not exist", {
  expect_error(
    preprocess_inputs_cdna(bam_file = "/nonexistent/file.bam",
                           dorado_summary = "summary.txt",
                           pod5_dir = tempdir(),
                           num_cores = 1,
                           qc = TRUE,
                           save_dir = tempdir(),
                           prefix = "",
                           part_size = 100,
                           cli_log = dummy_cli_log_cdna),
    "BAM"
  )
})

test_that("preprocess_inputs_cdna errors when dorado_summary is empty data frame", {
  tmp_bam <- tempfile(fileext = ".bam")
  file.create(tmp_bam)
  on.exit(unlink(tmp_bam), add = TRUE)

  expect_error(
    preprocess_inputs_cdna(bam_file = tmp_bam,
                           dorado_summary = data.frame(),
                           pod5_dir = tempdir(),
                           num_cores = 1,
                           qc = TRUE,
                           save_dir = tempdir(),
                           prefix = "",
                           part_size = 100,
                           cli_log = dummy_cli_log_cdna),
    "non-empty data frame"
  )
})


################################################################################
# process_polya_reads_cdna — validation guards only
################################################################################

test_that("process_polya_reads_cdna errors on missing polya_sequences", {
  expect_error(
    process_polya_reads_cdna(signal_files = character(0),
                             save_dir = tempdir()),
    "PolyA sequences are missing"
  )
})

test_that("process_polya_reads_cdna errors on missing signal_files", {
  expect_error(
    process_polya_reads_cdna(polya_sequences = data.frame(read_id = "r1"),
                             save_dir = tempdir()),
    "Signal files are missing"
  )
})

test_that("process_polya_reads_cdna errors on missing save_dir", {
  expect_error(
    process_polya_reads_cdna(polya_sequences = data.frame(read_id = "r1"),
                             signal_files = character(0)),
    "Output directory is missing"
  )
})

test_that("process_polya_reads_cdna errors when polya_sequences is not a data frame", {
  expect_error(
    process_polya_reads_cdna(polya_sequences = "not_a_df",
                             signal_files = character(0),
                             save_dir = tempdir()),
    "polya_sequences must be a data frame"
  )
})

test_that("process_polya_reads_cdna errors when polya_sequences lacks read_id column", {
  expect_error(
    process_polya_reads_cdna(polya_sequences = data.frame(wrong = "r1"),
                             signal_files = character(0),
                             save_dir = tempdir()),
    "read_id.*column"
  )
})

test_that("process_polya_reads_cdna errors when signal_files do not exist", {
  expect_error(
    process_polya_reads_cdna(polya_sequences = data.frame(read_id = "r1"),
                             signal_files = "/nonexistent/file.rds",
                             save_dir = tempdir()),
    "All signal files must exist"
  )
})

test_that("process_polya_reads_cdna returns empty result when no matching signals found", {
  skip_if_not_installed("vroom")

  tmp_rds <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp_rds), add = TRUE)

  # Signal list has a different read ID than polya_sequences
  saveRDS(list("other-read-001" = as.integer(rnorm(200, 680, 15))), tmp_rds)

  result <- process_polya_reads_cdna(
    polya_sequences = data.frame(read_id = "no-match-read", stringsAsFactors = FALSE),
    signal_files = tmp_rds,
    save_dir = tempdir(),
    cli_log = dummy_cli_log_cdna
  )

  expect_type(result, "list")
  expect_true("read_classes" %in% names(result))
  expect_true("nonadenosine_residues" %in% names(result))
  expect_equal(nrow(result$read_classes), 0L)
  expect_equal(result$tail_type, "polyA")
})


################################################################################
# process_polyt_reads_cdna — validation guards only
################################################################################

test_that("process_polyt_reads_cdna errors on missing polyt_sequences", {
  expect_error(
    process_polyt_reads_cdna(signal_files = character(0),
                             save_dir = tempdir()),
    "PolyT sequences are missing"
  )
})

test_that("process_polyt_reads_cdna errors on missing signal_files", {
  expect_error(
    process_polyt_reads_cdna(polyt_sequences = data.frame(read_id = "r1"),
                             save_dir = tempdir()),
    "Signal files are missing"
  )
})

test_that("process_polyt_reads_cdna errors on missing save_dir", {
  expect_error(
    process_polyt_reads_cdna(polyt_sequences = data.frame(read_id = "r1"),
                             signal_files = character(0)),
    "Output directory is missing"
  )
})

test_that("process_polyt_reads_cdna returns empty result when no matching signals found", {
  skip_if_not_installed("vroom")

  tmp_rds <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp_rds), add = TRUE)

  saveRDS(list("other-read-001" = as.integer(rnorm(200, 680, 15))), tmp_rds)

  result <- process_polyt_reads_cdna(
    polyt_sequences = data.frame(read_id = "no-match-read", stringsAsFactors = FALSE),
    signal_files = tmp_rds,
    save_dir = tempdir(),
    cli_log = dummy_cli_log_cdna
  )

  expect_type(result, "list")
  expect_true("read_classes" %in% names(result))
  expect_equal(nrow(result$read_classes), 0L)
  expect_equal(result$tail_type, "polyT")
})


################################################################################
# merge_cdna_results
################################################################################

test_that("merge_cdna_results errors on missing save_dir", {
  expect_error(merge_cdna_results(), "Output directory is missing")
})

test_that("merge_cdna_results errors when save_dir is not a character string", {
  expect_error(merge_cdna_results(save_dir = 123), "save_dir must be a character")
})

test_that("merge_cdna_results returns empty output when both polya and polyt are NULL", {
  result <- merge_cdna_results(save_dir = tempdir(),
                               cli_log = dummy_cli_log_cdna)
  expect_type(result, "list")
  expect_true("read_classes" %in% names(result))
  expect_true("nonadenosine_residues" %in% names(result))
  expect_equal(nrow(result$read_classes), 0L)
  expect_equal(nrow(result$nonadenosine_residues), 0L)
})

test_that("merge_cdna_results merges polyA and polyT read_classes correctly", {
  polya_classes <- data.frame(readname = c("rA-1", "rA-2"),
                              class = c("decorated", "blank"),
                              tail_type = "polyA",
                              stringsAsFactors = FALSE)
  polyt_classes <- data.frame(readname = c("rT-1", "rT-2"),
                              class = c("blank", "unclassified"),
                              tail_type = "polyT",
                              stringsAsFactors = FALSE)

  polya_res <- list(read_classes = polya_classes,
                    nonadenosine_residues = data.frame())
  polyt_res <- list(read_classes = polyt_classes,
                    nonadenosine_residues = data.frame())

  result <- merge_cdna_results(polya_results = polya_res,
                               polyt_results = polyt_res,
                               save_dir = tempdir(),
                               cli_log = dummy_cli_log_cdna)

  expect_equal(nrow(result$read_classes), 4L)
  expect_true(all(c("polyA", "polyT") %in% result$read_classes$tail_type))
})

test_that("merge_cdna_results adds tail_type column when absent from input", {
  # read_classes without tail_type — the function should add it
  polya_classes <- data.frame(readname = "rA-1",
                              class = "decorated",
                              stringsAsFactors = FALSE)
  polya_res <- list(read_classes = polya_classes,
                    nonadenosine_residues = data.frame())

  result <- merge_cdna_results(polya_results = polya_res,
                               save_dir = tempdir(),
                               cli_log = dummy_cli_log_cdna)

  expect_true("tail_type" %in% colnames(result$read_classes))
  expect_equal(result$read_classes$tail_type, "polyA")
})

test_that("merge_cdna_results merges nonadenosine_residues from both sources", {
  polya_mods <- data.frame(readname = "rA-1",
                           prediction = "C",
                           est_nonA_pos = 20L,
                           tail_type = "polyA",
                           stringsAsFactors = FALSE)
  polyt_mods <- data.frame(readname = "rT-1",
                           prediction = "G",
                           est_nonA_pos = 30L,
                           tail_type = "polyT",
                           stringsAsFactors = FALSE)

  polya_res <- list(read_classes = data.frame(),
                    nonadenosine_residues = polya_mods)
  polyt_res <- list(read_classes = data.frame(),
                    nonadenosine_residues = polyt_mods)

  result <- merge_cdna_results(polya_results = polya_res,
                               polyt_results = polyt_res,
                               save_dir = tempdir(),
                               cli_log = dummy_cli_log_cdna)

  expect_equal(nrow(result$nonadenosine_residues), 2L)
  expect_true(all(c("polyA", "polyT") %in% result$nonadenosine_residues$tail_type))
})


################################################################################
# save_cdna_outputs
################################################################################

test_that("save_cdna_outputs errors on missing outputs", {
  expect_error(save_cdna_outputs(save_dir = tempdir()),
               "Outputs are missing")
})

test_that("save_cdna_outputs errors on missing save_dir", {
  expect_error(save_cdna_outputs(outputs = list()),
               "Output directory is missing")
})

test_that("save_cdna_outputs errors when outputs is not a list", {
  expect_error(save_cdna_outputs(outputs = "not_a_list",
                                 save_dir = tempdir()),
               "outputs must be a list")
})

test_that("save_cdna_outputs returns standard ninetails format from complete input", {
  outputs <- list(
    read_classes = data.frame(readname = "r1", class = "decorated",
                              stringsAsFactors = FALSE),
    nonadenosine_residues = data.frame(readname = "r1", prediction = "C",
                                       stringsAsFactors = FALSE)
  )

  result <- save_cdna_outputs(outputs = outputs, save_dir = tempdir())

  expect_type(result, "list")
  expect_true("read_classes" %in% names(result))
  expect_true("nonadenosine_residues" %in% names(result))
  expect_equal(nrow(result$read_classes), 1L)
  expect_equal(nrow(result$nonadenosine_residues), 1L)
})

test_that("save_cdna_outputs returns empty data frames when outputs have NULL components", {
  result <- save_cdna_outputs(outputs = list(read_classes = NULL,
                                             nonadenosine_residues = NULL),
                              save_dir = tempdir())

  expect_true(is.data.frame(result$read_classes))
  expect_true(is.data.frame(result$nonadenosine_residues))
  expect_equal(nrow(result$read_classes), 0L)
  expect_equal(nrow(result$nonadenosine_residues), 0L)
})
