################################################################################
# Testing ninetails_statistic functions
################################################################################


# Helper: Create test data suitable for Fisher's test
# Requires: grouping_factor with 2 levels, counts for blank/decorated reads
create_fisher_test_data <- function() {
  data.frame(
    readname = paste0("read_", 1:40),
    contig = rep(c("tx_A", "tx_B"), each = 20),
    polya_length = sample(50:150, 40, replace = TRUE),
    qc_tag = rep("PASS", 40),
    class = c(rep(c("decorated", "blank"), 10), rep(c("decorated", "blank"), 10)),
    comments = c(rep(c("YAY", "MPU"), 10), rep(c("YAY", "MPU"), 10)),
    group = factor(rep(c("WT", "KO"), each = 20)),
    sample_name = factor(rep(c("WT", "KO"), each = 20)),
    prediction_C = c(sample(0:2, 20, replace = TRUE), sample(0:3, 20, replace = TRUE)),
    prediction_G = c(sample(0:1, 20, replace = TRUE), sample(0:2, 20, replace = TRUE)),
    prediction_U = c(sample(0:1, 20, replace = TRUE), sample(0:1, 20, replace = TRUE)),
    nonA_residues = c(
      rep(c("G15", NA, "C20:G30", NA), 5),
      rep(c("G10:C25", NA, "U40", NA), 5)
    ),
    ensembl_transcript_id_short = rep(c("ENST001", "ENST002"), each = 20),
    stringsAsFactors = FALSE
  )
}


# nonA_fisher
################################################################################

test_that("nonA_fisher validates missing ninetails_data", {
  expect_error(nonA_fisher(grouping_factor = "group",
                           base = "C",
                           transcript_id_column = "contig"),"Ninetails data are missing")
})


test_that("nonA_fisher validates missing base argument", {
  test_df <- create_fisher_test_data()

  expect_error(
    nonA_fisher(ninetails_data = test_df,
                grouping_factor = "group",
                transcript_id_column = "contig"),"Base is missing")
})


test_that("nonA_fisher validates missing transcript_id_column", {
  test_df <- create_fisher_test_data()

  expect_error(
    nonA_fisher(ninetails_data = test_df,
                grouping_factor = "group",
                base = "C"),"Transcript_id_column is missing")
})


test_that("nonA_fisher validates empty data frame", {
  empty_df <- data.frame(group = factor(levels = c("WT", "KO")))

  expect_error(
    nonA_fisher(ninetails_data = empty_df,
                grouping_factor = "group",
                base = "C",
                transcript_id_column = "contig"),"Empty data frame")
})


test_that("nonA_fisher validates non-numeric min_reads", {
  test_df <- create_fisher_test_data()

  expect_error(
    nonA_fisher(
      ninetails_data = test_df,
      grouping_factor = "group",
      base = "C",
      min_reads = "ten",
      transcript_id_column = "contig"),"Non-numeric parameter.*min_reads")
})


test_that("nonA_fisher validates grouping_factor exists in data", {
  test_df <- create_fisher_test_data()

  expect_error(
    nonA_fisher(
      ninetails_data = test_df,
      grouping_factor = "nonexistent_column",
      base = "C",
      transcript_id_column = "contig"),"is not a column")
})


test_that("nonA_fisher validates grouping_factor has more than 1 level", {
  test_df <- create_fisher_test_data()
  # Force single level
  test_df$group <- factor(rep("WT", nrow(test_df)))

  expect_error(
    nonA_fisher(
      ninetails_data = test_df,
      grouping_factor = "group",
      base = "C",
      transcript_id_column = "contig"),"Only 1 level")
})


test_that("nonA_fisher validates base argument values", {
  test_df <- create_fisher_test_data()

  expect_error(
    nonA_fisher(
      ninetails_data = test_df,
      grouping_factor = "group",
      base = "X",  # Invalid base
      transcript_id_column = "contig"),"Wrong non-A nucleotide")
})


test_that("nonA_fisher accepts valid base values (C, G, U, all)", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tibble")

  test_df <- create_fisher_test_data()

  # Test each valid base
  for (base_val in c("C", "G", "U", "all")) {
    result <- nonA_fisher(ninetails_data = test_df,
                          grouping_factor = "group",
                          base = base_val,
                          transcript_id_column = "contig")

    expect_s3_class(result, "tbl_df")
    expect_true("p.value" %in% colnames(result))
    expect_true("stats_code" %in% colnames(result))
  }
})


test_that("nonA_fisher returns tibble with correct structure", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tibble")

  test_df <- create_fisher_test_data()

  result <- nonA_fisher(ninetails_data = test_df,
                        grouping_factor = "group",
                        base = "C",
                        transcript_id_column = "contig")

  expect_s3_class(result, "tbl_df")
  expect_equal(ncol(result), 2)
  expect_true("p.value" %in% colnames(result))
  expect_true("stats_code" %in% colnames(result))
})


test_that("nonA_fisher returns appropriate stats_code for low counts", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tibble")

  test_df <- create_fisher_test_data()

  # Set high min_reads to trigger low count codes
  result <- nonA_fisher(
    ninetails_data = test_df,
    grouping_factor = "group",
    base = "C",
    min_reads = 1000,  # Higher than available reads
    transcript_id_column = "contig")

  expect_true(result$stats_code %in% c("B_LC", "G_LC", "G_NA", "B_NA"))
})


################################################################################
# calculate_fisher - Input Validation
################################################################################

test_that("calculate_fisher validates missing ninetails_data", {
  expect_error(
    calculate_fisher(
      transcript_id_column = "contig",
      min_reads = 10,
      min_nonA_reads = 5,
      base = "C"),"Ninetails data are missing")
})


test_that("calculate_fisher validates missing transcript_id_column", {
  test_df <- create_fisher_test_data()

  expect_error(
    calculate_fisher(ninetails_data = test_df,
                     min_reads = 10,
                     min_nonA_reads = 5,
                     base = "C"),"Transcript_id_column is missing")
})


test_that("calculate_fisher validates missing min_reads", {
  test_df <- create_fisher_test_data()

  expect_error(
    calculate_fisher(
      ninetails_data = test_df,
      transcript_id_column = "contig",
      min_nonA_reads = 5,
      base = "C"),"Min_reads are missing")
})


test_that("calculate_fisher validates missing min_nonA_reads", {
  test_df <- create_fisher_test_data()

  expect_error(
    calculate_fisher(
      ninetails_data = test_df,
      transcript_id_column = "contig",
      min_reads = 10,
      base = "C"),"Min_nonA_reads.*are missing")
})


test_that("calculate_fisher validates missing base", {
  test_df <- create_fisher_test_data()

  expect_error(
    calculate_fisher(
      ninetails_data = test_df,
      transcript_id_column = "contig",
      min_reads = 10,
      min_nonA_reads = 5),"Base definition is missing" )
})


test_that("calculate_fisher validates empty data frame", {
  empty_df <- data.frame(group = factor(levels = c("WT", "KO")))

  expect_error(
    calculate_fisher(
      ninetails_data = empty_df,
      transcript_id_column = "contig",
      min_reads = 10,
      min_nonA_reads = 5,
      base = "C"),"Empty data frame")
})


test_that("calculate_fisher validates numeric parameters", {
  test_df <- create_fisher_test_data()

  # Non-numeric min_reads
  expect_error(
    calculate_fisher(ninetails_data = test_df,
                     transcript_id_column = "contig",
                     min_reads = "ten",
                     min_nonA_reads = 5,
                     base = "C"),
    "Min_reads must be numeric")

  # Non-numeric min_nonA_reads
  expect_error(calculate_fisher(ninetails_data = test_df,
                                transcript_id_column = "contig",
                                min_reads = 10,
                                min_nonA_reads = "five",
                                base = "C"),
               "Min_nonA_reads must be numeric")

  # Non-numeric alpha
  expect_error(calculate_fisher(ninetails_data = test_df,
                                transcript_id_column = "contig",
                                min_reads = 10,
                                min_nonA_reads = 5,
                                base = "C",
                                alpha = "Putin chuj"),
               "Alpha must be numeric")
})


test_that("calculate_fisher validates grouping_factor has exactly 1 level fails", {
  test_df <- create_fisher_test_data()
  test_df$group <- factor(rep("WT", nrow(test_df)))

  expect_error(
    calculate_fisher(
      ninetails_data = test_df,
      transcript_id_column = "contig",
      min_reads = 0,
      min_nonA_reads = 0,
      grouping_factor = "group",
      base = "C"),"Only 1 level")
})


test_that("calculate_fisher requires conditions for >2 level grouping_factor", {
  test_df <- create_fisher_test_data()
  # Add third level
  extra_rows <- test_df[1:10, ]
  extra_rows$group <- factor("HET")
  test_df <- rbind(test_df, extra_rows)
  test_df$group <- factor(test_df$group, levels = c("WT", "KO", "HET"))

  expect_error(
    calculate_fisher(
      ninetails_data = test_df,
      transcript_id_column = "contig",
      min_reads = 0,
      min_nonA_reads = 0,
      grouping_factor = "group",
      base = "C"),"more than 2 levels")
})


test_that("calculate_fisher validates condition1 and condition2 exist", {
  test_df <- create_fisher_test_data()
  # Add third level
  extra_rows <- test_df[1:10, ]
  extra_rows$group <- factor("HET")
  test_df <- rbind(test_df, extra_rows)
  test_df$group <- factor(test_df$group, levels = c("WT", "KO", "HET"))

  expect_error(
    calculate_fisher(
      ninetails_data = test_df,
      transcript_id_column = "contig",
      min_reads = 0,
      min_nonA_reads = 0,
      grouping_factor = "group",
      condition1 = "WT",
      condition2 = "NICNIMOZNACZYNIO",  # Invalid level
      base = "C"), "is not a level")
})


test_that("calculate_fisher validates condition1 != condition2", {
  test_df <- create_fisher_test_data()
  # Add third level
  extra_rows <- test_df[1:10, ]
  extra_rows$group <- factor("HET")
  test_df <- rbind(test_df, extra_rows)
  test_df$group <- factor(test_df$group, levels = c("WT", "KO", "HET"))

  expect_error(
    calculate_fisher(
      ninetails_data = test_df,
      transcript_id_column = "contig",
      min_reads = 0,
      min_nonA_reads = 0,
      grouping_factor = "group",
      condition1 = "WT",
      condition2 = "WT",
      base = "C"),"should be different") # same conditions
})



test_that("calculate_fisher returns tibble with expected columns", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("purrr")
  skip_if_not_installed("stringr")
  skip_if_not_installed("tibble")

  test_df <- create_fisher_test_data()

  result <- suppressMessages(
    calculate_fisher(
      ninetails_data = test_df,
      transcript_id_column = "ensembl_transcript_id_short",
      min_reads = 0,
      min_nonA_reads = 0,
      grouping_factor = "group",
      base = "C",
      alpha = 0.05))

  expect_s3_class(result, "tbl_df")
  expect_true("p.value" %in% colnames(result))
  expect_true("stats_code" %in% colnames(result))
  expect_true("padj" %in% colnames(result))
  expect_true("significance" %in% colnames(result))
  expect_true("ensembl_transcript_id_short" %in% colnames(result))
  expect_true(all(result$significance %in% c("FDR<0.05", "NotSig")))
  expect_s3_class(result, "tbl_df")
  expect_true(nrow(result) > 0)
})



test_that("calculate_fisher filters by min_nonA_reads", {

  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("purrr")
  skip_if_not_installed("stringr")
  skip_if_not_installed("tibble")

  test_df <- create_fisher_test_data()

  # With very high min_nonA_reads, should filter out all transcripts
  result <- suppressMessages(suppressWarnings(
    calculate_fisher(
      ninetails_data = test_df,
      transcript_id_column = "ensembl_transcript_id_short",
      min_reads = 0,
      min_nonA_reads = 1000,  # Very high
      grouping_factor = "group",
      base = "C",
      alpha = 0.05)))

  # Result should have 0 rows (all filtered out)
  # Note: When all transcripts are filtered, the result may be an empty tibble
  # which triggers warnings about uninitialized columns - this is expected behavior
  expect_equal(nrow(result), 0)
})





