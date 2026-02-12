# ################################################################################
# # Testing ninetails_statistic
# ################################################################################
#
#
# # nonA_fisher
# ################################################################################
#
# test_that("nonA_fisher errors on empty data frame", {
#   empty_df <- data.frame(group = factor(levels = c("WT", "KO")))
#
#   expect_error(nonA_fisher(ninetails_data = empty_df,
#                            grouping_factor = "group",
#                            base = "C",
#                            transcript_id_column = "contig"), "Empty data frame")
# })
#
# test_that("nonA_fisher errors on non-numeric min_reads", {
#   dummy_df <- data.frame(group = factor(c("WT", "KO")),
#                          contig = c("tx1", "tx1"),
#                          stringsAsFactors = FALSE)
#
#   expect_error(nonA_fisher(ninetails_data = dummy_df,
#                            grouping_factor = "group",
#                            base = "C",
#                            min_reads = "ten",
#                            transcript_id_column = "contig"),"Non-numeric parameter.*min_reads")
# })
#
# test_that("nonA_fisher errors when grouping factor has only 1 level", {
#   dummy_df <- data.frame(group = factor(c("WT", "KO", "HET"), levels = c("WT", "KO", "HET")),
#                          contig = c("tx1", "tx1", "tx1"),
#                          stringsAsFactors = FALSE)
#
#   expect_error(nonA_fisher(ninetails_data = dummy_df,
#                            grouping_factor = "group",
#                            base = "C",
#                            transcript_id_column = "contig"),"only 1 level")
# })
