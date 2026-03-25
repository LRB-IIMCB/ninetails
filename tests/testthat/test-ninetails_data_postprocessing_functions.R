################################################################################
# Testing data postprocessing functions
################################################################################


# read_class_single
################################################################################

test_that("read_class_single errors on missing class_path", {
  expect_error(read_class_single(),
               "class_path.*is missing")
})

test_that("read_class_single errors on NA class_path", {
  expect_error(read_class_single(class_path = NA),
               "missing or invalid")
})

test_that("read_class_single errors on empty string class_path", {
  expect_error(read_class_single(class_path = ""),
               "missing or invalid")
})

test_that("read_class_single errors on nonexistent file", {
  expect_error(read_class_single(class_path = "/nonexistent/path/file.tsv"),
               "does not exist")
})

test_that("read_class_single errors on empty file", {
  # Create empty temp file
  tmp_file <- tempfile(fileext = ".tsv")
  file.create(tmp_file)
  on.exit(unlink(tmp_file))

  expect_error(read_class_single(class_path = tmp_file),
               "empty")
})





# read_class_multiple
################################################################################

test_that("read_class_multiple errors on missing samples_table", {
  expect_error(read_class_multiple(),
               "Samples table.*is missing")
})

test_that("read_class_multiple errors on empty samples_table", {
  expect_error(read_class_multiple(samples_table = data.frame()),
               "Empty data frame")
})

test_that("read_class_multiple errors on missing class_path column", {
  bad_table <- data.frame(sample_name = "WT_1")
  expect_error(read_class_multiple(samples_table = bad_table),
               "class_path.*sample_name")
})

test_that("read_class_multiple errors on missing sample_name column", {
  bad_table <- data.frame(class_path = "/some/path.tsv")
  expect_error(read_class_multiple(samples_table = bad_table),
               "class_path.*sample_name")
})


# count_class
################################################################################

test_that("count_class validates empty data frame", {
  expect_error(count_class(class_data = data.frame()),
               "Empty data frame")
})


test_that("count_class validates grouping_factor exists", {
  skip_if(!exists("class_data_single"), "class_data_single not loaded")

  expect_error(
    count_class(class_data = class_data_single,
                grouping_factor = "nonexistent_column"),
    "is not a column")
})


test_that("count_class returns detailed counts (comments) without grouping", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("forcats")
  skip_if(!exists("class_data_single"), "class_data_single not loaded")

  result <- count_class(class_data = class_data_single,
                        detailed = TRUE)

  expect_s3_class(result, "tbl_df")
  expect_true("comments" %in% colnames(result))
  expect_true("n" %in% colnames(result))
  expect_equal(ncol(result), 2)

})


test_that("count_class returns crude counts (class) without grouping", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("forcats")
  skip_if(!exists("class_data_single"), "class_data_single not loaded")

  result <- count_class(class_data = class_data_single,
                        detailed = FALSE)

  expect_s3_class(result, "tbl_df")
  expect_true("class" %in% colnames(result))
  expect_true("n" %in% colnames(result))
  expect_equal(ncol(result), 2)

})


test_that("count_class returns detailed counts with grouping", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("forcats")
  skip_if(!exists("class_data_grouped"), "class_data_grouped not loaded")

  result <- count_class(class_data = class_data_grouped,
                        grouping_factor = "sample_name",
                        detailed = TRUE)

  expect_s3_class(result, "tbl_df")
  expect_true("sample_name" %in% colnames(result))
  expect_true("comments" %in% colnames(result))
  expect_true("n" %in% colnames(result))
  expect_equal(ncol(result), 3)

})


test_that("count_class returns crude counts with grouping", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("forcats")
  skip_if(!exists("class_data_grouped"), "class_data_grouped not loaded")

  result <- count_class(class_data = class_data_grouped,
                        grouping_factor = "group",
                        detailed = FALSE)

  expect_s3_class(result, "tbl_df")
  expect_true("group" %in% colnames(result))
  expect_true("class" %in% colnames(result))
  expect_true("n" %in% colnames(result))

})


test_that("count_class reorders class factor with decorated last", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("forcats")
  skip_if(!exists("class_data_single"), "class_data_single not loaded")

  result <- count_class(class_data = class_data_single, detailed = FALSE)

  # decorated should be last level
  levels_order <- levels(result$class)
  expect_equal(levels_order[length(levels_order)], "decorated")
})

# read_residue_single
################################################################################

test_that("read_residue_single errors on NA residue_path", {
  expect_error(read_residue_single(residue_path = NA),
               "missing or an empty string")
})

test_that("read_residue_single errors on empty string residue_path", {
  expect_error(read_residue_single(residue_path = ""),
               "missing or an empty string")
})

test_that("read_residue_single errors on nonexistent file", {
  expect_error(read_residue_single(residue_path = "/nonexistent/path/file.tsv"),
               "missing or an empty file")
})

test_that("read_residue_single errors on empty file", {
  tmp_file <- tempfile(fileext = ".tsv")
  file.create(tmp_file)
  on.exit(unlink(tmp_file))

  expect_error(read_residue_single(residue_path = tmp_file),
               "missing or an empty file")
})

test_that("read_residue_single reads valid file", {
  tmp_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp_file))

  test_df <- data.frame(
    readname = c("read_1", "read_2"),
    contig = c("ENST00000123|a|b|c|d|GENE1|e", "ENST00000456|a|b|c|d|GENE2|e"),
    prediction = c("G", "C"),
    est_nonA_pos = c(50, 30),
    polya_length = c(100, 80),
    qc_tag = c("PASS", "PASS"),
    stringsAsFactors = FALSE
  )
  utils::write.table(test_df, tmp_file, sep = "\t", row.names = FALSE, quote = FALSE)

  result <- suppressMessages(read_residue_single(residue_path = tmp_file))

  expect_s3_class(result, "tbl_df")
  expect_true("readname" %in% colnames(result))
  expect_true("transcript" %in% colnames(result))
  expect_true("ensembl_transcript_id_full" %in% colnames(result))
  expect_true("ensembl_transcript_id_short" %in% colnames(result))
})

test_that("read_residue_single adds sample_name when provided", {
  tmp_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp_file))

  test_df <- data.frame(
    readname = c("read_1"),
    contig = c("gene1"),
    prediction = c("G"),
    est_nonA_pos = c(50),
    polya_length = c(100),
    qc_tag = c("PASS"),
    stringsAsFactors = FALSE
  )
  utils::write.table(test_df, tmp_file, sep = "\t", row.names = FALSE, quote = FALSE)

  result <- suppressMessages(suppressWarnings(
    read_residue_single(residue_path = tmp_file, sample_name = "KO_1")
  ))

  expect_true("sample_name" %in% colnames(result))
  expect_equal(as.character(result$sample_name[1]), "KO_1")
  expect_true(is.factor(result$sample_name))
})


# read_residue_multiple
################################################################################

test_that("read_residue_multiple errors on missing samples_table", {
  expect_error(read_residue_multiple(),
               "samples_table.*is missing")
})

test_that("read_residue_multiple errors on empty samples_table", {
  expect_error(read_residue_multiple(samples_table = data.frame()),
               "Empty data frame")
})

test_that("read_residue_multiple errors on missing residue_path column", {
  bad_table <- data.frame(sample_name = "WT_1")
  expect_error(read_residue_multiple(samples_table = bad_table),
               "residue_path.*sample_name")
})

test_that("read_residue_multiple errors on missing sample_name column", {
  bad_table <- data.frame(residue_path = "/some/path.tsv")
  expect_error(read_residue_multiple(samples_table = bad_table),
               "residue_path.*sample_name")
})
# count_residues
################################################################################

test_that("count_residues validates empty data frame", {
  expect_error(count_residues(residue_data = data.frame()),"Empty data frame")
})


test_that("count_residues validates grouping_factor exists", {
  skip_if(!exists("residue_data_single"), "residue_data_single not loaded")

  expect_error(count_residues(
    residue_data = residue_data_single,
    grouping_factor = "nonexistent_column"),
    "is not a column")
})


test_that("count_residues returns counts without grouping", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("forcats")
  skip_if(!exists("residue_data_single"), "residue_data_single not loaded")

  result <- count_residues(residue_data = residue_data_single)

  expect_s3_class(result, "tbl_df")
  expect_true("prediction" %in% colnames(result))
  expect_true("n" %in% colnames(result))
  expect_equal(ncol(result), 2)

  # Should have counts for C, G, U
  expect_true(all(c("C", "G", "U") %in% result$prediction))

})


test_that("count_residues returns counts with grouping", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("forcats")
  skip_if(!exists("residue_data_grouped"), "residue_data_grouped not loaded")

  result <- count_residues(residue_data = residue_data_grouped,
                           grouping_factor = "group")

  expect_s3_class(result, "tbl_df")
  expect_true("group" %in% colnames(result))
  expect_true("prediction" %in% colnames(result))
  expect_true("n" %in% colnames(result))
  expect_equal(ncol(result), 3)
})


test_that("count_residues works with sample_name grouping", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("forcats")
  skip_if(!exists("residue_data_grouped"), "residue_data_grouped not loaded")

  result <- count_residues(residue_data = residue_data_grouped,
                           grouping_factor = "sample_name")

  expect_s3_class(result, "tbl_df")
  expect_true("sample_name" %in% colnames(result))
  expect_equal(ncol(result), 3)
})



# spread_nonA_residues
################################################################################

test_that("spread_nonA_residues validates required arguments", {
  expect_error(
    spread_nonA_residues(),
    "is missing"
  )

  expect_error(spread_nonA_residues(residue_data = data.frame()), "Empty data frame")
})


test_that("spread_nonA_residues validates required columns", {
  bad_df <- data.frame(a = 1, b = 2)

  expect_error(spread_nonA_residues(residue_data = bad_df),"Required columns missing")
})


test_that("spread_nonA_residues returns wide format with prediction counts", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  skip_if(!exists("residue_data_grouped"), "residue_data_grouped not loaded")

  result <- spread_nonA_residues(residue_data = residue_data_grouped)

  expect_s3_class(result, "tbl_df")
  expect_true("readname" %in% colnames(result))
  expect_true("nonA_residues" %in% colnames(result))

  # Should have prediction_C, prediction_G, prediction_U columns
  prediction_cols <- grep("^prediction_", colnames(result), value = TRUE)
  expect_true(length(prediction_cols) > 0)
})


test_that("spread_nonA_residues creates correct nonA_residues string", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  skip_if(!exists("residue_data_grouped"), "residue_data_grouped not loaded")

  result <- spread_nonA_residues(residue_data = residue_data_grouped)

  # nonA_residues should be in format like "G38:C83" or "U33"
  nona_strings <- result$nonA_residues[!is.na(result$nonA_residues)]
  expect_true(length(nona_strings) > 0)

  # Each string should contain letters and numbers
  expect_true(all(grepl("[CGU]\\d+", nona_strings)))
})




# merge_nonA_tables
################################################################################

test_that("merge_nonA_tables validates required arguments", {
  expect_error( merge_nonA_tables(residue_data = data.frame(a = 1)),"class data is missing")

  expect_error(merge_nonA_tables(class_data = data.frame(a = 1)),"residue data is missing")

  expect_error(merge_nonA_tables(class_data = data.frame(),residue_data = data.frame(a = 1)),
               "Empty data frame.*class_data" )

  expect_error(merge_nonA_tables(class_data = data.frame(a = 1),residue_data = data.frame()),
               "Empty data frame.*residue_data")
})


test_that("merge_nonA_tables validates pass_only is logical", {
  skip_if(!exists("class_data_grouped"), "class_data_grouped not loaded")
  skip_if(!exists("residue_data_grouped"), "residue_data_grouped not loaded")

  expect_error(merge_nonA_tables(class_data = class_data_grouped,
                                 residue_data = residue_data_grouped,
                                 pass_only = "yes"), "TRUE/FALSE" )
})


test_that("merge_nonA_tables validates qc_tag column exists", {

  bad_class <- data.frame(
    readname = "read_1",
    class = "decorated",
    comments = "YAY",
    stringsAsFactors = FALSE)

  bad_residue <- data.frame(
    readname = "read_1",
    prediction = "G",
    est_nonA_pos = 10,
    group = "WT",
    stringsAsFactors = FALSE)

  expect_error(merge_nonA_tables(class_data = bad_class,
                                 residue_data = bad_residue), "qc_tag.*missing")
})


test_that("merge_nonA_tables works with numeric qc_tag (Dorado)", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  skip_if(!exists("class_data_grouped"), "class_data_grouped not loaded")
  skip_if(!exists("residue_data_grouped"), "residue_data_grouped not loaded")

  # class_data_grouped and residue_data_grouped have numeric qc_tag
  result <- suppressMessages(suppressWarnings(
    merge_nonA_tables(class_data = class_data_grouped,
                      residue_data = residue_data_grouped,
                      pass_only = TRUE)))

  expect_s3_class(result, "data.frame")
  expect_true("readname" %in% colnames(result))
  expect_true("class" %in% colnames(result))
  expect_true("nonA_residues" %in% colnames(result))

  # Should have prediction count columns
  expect_true(any(grepl("prediction_", colnames(result))))
})


test_that("merge_nonA_tables works with character qc_tag (Guppy)", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  skip_if(!exists("merged_nonA_char"), "merged_nonA_char not loaded")

  # Create class_data and residue_data with character qc_tag from merged_nonA_char
  # This is a workaround since we need character qc_tag data
  class_data_char <- data.frame(
    readname = c("read_1", "read_2", "read_3"),
    contig = c("beta_gal", "beta_gal", "beta_gal"),
    polya_length = c(100, 80, 60),
    qc_tag = c("PASS", "PASS", "SUFFCLIP"),
    class = c("decorated", "blank", "decorated"),
    comments = c("YAY", "MPU", "YAY"),
    group = c("WT", "WT", "KO"),
    stringsAsFactors = FALSE)

  residue_data_char <- data.frame(
    readname = c("read_1", "read_3"),
    contig = c("beta_gal", "beta_gal"),
    prediction = c("G", "C"),
    est_nonA_pos = c(50, 30),
    polya_length = c(100, 60),
    qc_tag = c("PASS", "SUFFCLIP"),
    group = c("WT", "KO"),
    stringsAsFactors = FALSE)

  # Test pass_only = TRUE
  result_pass <- suppressMessages(suppressWarnings(
    merge_nonA_tables(
      class_data = class_data_char,
      residue_data = residue_data_char,
      pass_only = TRUE)))

  expect_s3_class(result_pass, "data.frame")

  # Test pass_only = FALSE
  result_all <- suppressMessages(suppressWarnings(
    merge_nonA_tables(class_data = class_data_char,
                      residue_data = residue_data_char,
                      pass_only = FALSE)))

  expect_s3_class(result_all, "data.frame")
})


# summarize_nonA
################################################################################

test_that("summarize_nonA validates required arguments", {
  expect_error(
    summarize_nonA(),
    "merged_nonA_tables.*is missing"
  )

  expect_error(summarize_nonA(merged_nonA_tables = data.frame()),"Empty data frame")

})


test_that("summarize_nonA validates summary_factors type", {
  skip_if(!exists("merged_nonA_char"), "merged_nonA_char not loaded")

  expect_error(summarize_nonA(merged_nonA_tables = merged_nonA_char,
                              summary_factors = 123 ),"Non-character")

})


test_that("summarize_nonA validates summary_factors exist in data", {
  skip_if(!exists("merged_nonA_char"), "merged_nonA_char not loaded")

  expect_error(summarize_nonA(merged_nonA_tables = merged_nonA_char,
                              summary_factors = "nonexistent_column"),"Non-existent column")
})


test_that("summarize_nonA returns transcript-level summary", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("stringr")
  skip_if(!exists("merged_nonA_char"), "merged_nonA_char not loaded")

  # Add required transcript column for testing
  test_data <- merged_nonA_char
  test_data$ensembl_transcript_id_short <- test_data$contig

  result <- summarize_nonA(merged_nonA_tables = test_data,
                           summary_factors = "group",
                           transcript_id_column = "ensembl_transcript_id_short")

  expect_s3_class(result, "tbl_df")
  expect_true("ensembl_transcript_id_short" %in% colnames(result))
  expect_true("group" %in% colnames(result))
  expect_true("polya_median" %in% colnames(result))
  expect_true("polya_mean" %in% colnames(result))
  expect_true("counts_total" %in% colnames(result))
  expect_true("counts_blank" %in% colnames(result))
})




# nanopolish_qc
################################################################################

test_that("nanopolish_qc validates empty data frame", {
  expect_error( nanopolish_qc(class_data = data.frame()),"Empty data frame")
})


test_that("nanopolish_qc validates grouping_factor exists", {
  skip_if(!exists("class_data_single"), "class_data_single not loaded")

  expect_error( nanopolish_qc(class_data = class_data_single,
                              grouping_factor = "nonexistent_column"),"is not a column")

})


test_that("nanopolish_qc returns qc counts without grouping", {
  skip_if_not_installed("dplyr")
  skip_if(!exists("class_data_single"), "class_data_single not loaded")

  # Create class_data with character qc_tag for proper testing
  class_data_qc <- data.frame(
    readname = paste0("read_", 1:10),
    contig = rep("beta_gal", 10),
    polya_length = sample(50:150, 10),
    qc_tag = c(rep("PASS", 6), rep("ADAPTER", 2), rep("SUFFCLIP", 2)),
    class = c(rep("decorated", 5), rep("blank", 3), rep("unclassified", 2)),
    comments = c(rep("YAY", 5), rep("MPU", 3), rep("IRL", 2)),
    stringsAsFactors = FALSE)

  result <- nanopolish_qc(class_data = class_data_qc)

  expect_s3_class(result, "tbl_df")
  expect_true("qc_tag" %in% colnames(result))
  expect_true("n" %in% colnames(result))
  expect_equal(ncol(result), 2)
})


test_that("nanopolish_qc returns qc counts with grouping", {
  skip_if_not_installed("dplyr")

  # Create grouped class_data with character qc_tag
  class_data_qc <- data.frame(
    readname = paste0("read_", 1:20),
    contig = rep("beta_gal", 20),
    polya_length = sample(50:150, 20),
    qc_tag = rep(c("PASS", "PASS", "ADAPTER", "SUFFCLIP"), 5),
    class = rep(c("decorated", "blank", "unclassified", "decorated"), 5),
    comments = rep(c("YAY", "MPU", "IRL", "YAY"), 5),
    sample_name = rep(c("WT_1", "KO_1"), each = 10),
    stringsAsFactors = FALSE)

  result <- nanopolish_qc(class_data = class_data_qc,grouping_factor = "sample_name")

  expect_s3_class(result, "tbl_df")
  expect_true("sample_name" %in% colnames(result))
  expect_true("qc_tag" %in% colnames(result))
  expect_true("n" %in% colnames(result))
  expect_equal(ncol(result), 3)
})



# correct_residue_data
################################################################################

test_that("correct_residue_data validates required arguments", {
  expect_error(correct_residue_data(residue_data = data.frame(a = 1),
                                    transcript_column = "contig" ), "class_data argument is missing")

  expect_error(correct_residue_data(class_data = data.frame(a = 1),
                                    transcript_column = "contig"),"residue_data argument is missing")

  expect_error(correct_residue_data(class_data = data.frame(a = 1),
                                    residue_data = data.frame(b = 2)),"transcript_column argument is missing")

})


test_that("correct_residue_data validates empty data frames", {
  skip_if(!exists("residue_data_single"), "residue_data_single not loaded")

  expect_error(correct_residue_data(class_data = data.frame(),
                                    residue_data = residue_data_single,
                                    transcript_column = "contig"),"Empty data frame.*class_data")

  skip_if(!exists("class_data_single"), "class_data_single not loaded")

  expect_error(correct_residue_data(class_data = class_data_single,
                                    residue_data = data.frame(),
                                    transcript_column = "contig"),"Empty data frame.*residue_data")
})



# correct_class_data
################################################################################

test_that("correct_class_data validates required arguments", {
  expect_error(correct_class_data(class_data = data.frame(a = 1)),
               "residue_data_edited argument is missing")

  expect_error(correct_class_data(residue_data_edited = data.frame(a = 1)),
               "class_data argument is missing")

})



# reclassify_ninetails_data
################################################################################

test_that("reclassify_ninetails_data validates required arguments", {

#  Missing residue_data (checked first)
expect_error(reclassify_ninetails_data(class_data = data.frame(a = 1),
                                       transcript_column = "contig"),
             "The residue_data argument is missing")

#  Missing class_data (residue_data must be present to reach this check)
expect_error(reclassify_ninetails_data(residue_data = data.frame(a = 1),
                                       transcript_column = "contig"),
             "The class_data argument is missing")

#  Missing transcript_column (both data inputs must be present)
expect_error(reclassify_ninetails_data(residue_data = data.frame(a = 1),
                                       class_data = data.frame(b = 2)),
             "The transcript_column argument is missing")

})


# count_nonA_abundance
################################################################################

test_that("count_nonA_abundance validates empty data frame", {
  expect_error(
    count_nonA_abundance(residue_data = data.frame()), "Empty data frame")
})


test_that("count_nonA_abundance validates grouping_factor exists", {
  skip_if(!exists("residue_data_single"), "residue_data_single not loaded")

  expect_error(
    count_nonA_abundance(residue_data = residue_data_single,
                         grouping_factor = "nonexistent_column"),"is not a column")

})


test_that("count_nonA_abundance returns abundance counts", {
  skip_if_not_installed("dplyr")
  skip_if(!exists("residue_data_grouped"), "residue_data_grouped not loaded")

  result <- count_nonA_abundance(residue_data = residue_data_grouped,
                                 grouping_factor = "group" )

  expect_s3_class(result, "tbl_df")
  expect_true("group" %in% colnames(result))
  expect_true("instances" %in% colnames(result))
  expect_true("count" %in% colnames(result))

  # instances should be "single", "two", or "more"
  expect_true(all(result$instances %in% c("single", "two", "more")))
})


test_that("count_nonA_abundance works with sample_name grouping", {
  skip_if_not_installed("dplyr")
  skip_if(!exists("residue_data_grouped"), "residue_data_grouped not loaded")

  result <- count_nonA_abundance(residue_data = residue_data_grouped,
                                 grouping_factor = "sample_name")

  expect_s3_class(result, "tbl_df")
  expect_true("sample_name" %in% colnames(result))
})

