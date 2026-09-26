# Counts peak and valley pseudomove runs per read in a training-set feature list (lookup for the cDNA training data).

In the Guppy/DRS training routine the polarity of the signal deviation
caused by a given residue was known (G produces a peak, C and U produce
valleys). The polarity for cDNA (DNA chemistry, two orientations, and
complementary bases in polyT reads) has to be established empirically
before
[`filter_nonA_chunks_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/filter_nonA_chunks_trainingset.md)
is called with the correct `value`. This function tabulates, for every
read, how many qualifying peak (+1) and valley (-1) runs the pseudomove
vector contains, so that the dominant polarity of a labelled dataset can
be read off the column sums.

## Usage

``` r
count_pseudomove_runs_trainingset_cdna(tail_feature_list, min_run_length = 4)
```

## Arguments

- tail_feature_list:

  List object produced by
  [`create_tail_feature_list_trainingset_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list_trainingset_cdna.md).

- min_run_length:

  Numeric `[4]`. Minimum length of a non-zero pseudomove run to be
  counted. The default matches the chunking criterion of
  [`split_tail_centered_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/split_tail_centered_trainingset.md).

## Value

A data frame with one row per read and columns:

- readname:

  Character. Read ID.

- signal_length:

  Integer. Length of the downsampled tail signal.

- peak_runs:

  Integer. Number of +1 runs of length \>= `min_run_length`.

- valley_runs:

  Integer. Number of -1 runs of length \>= `min_run_length`.

## See also

[`create_tail_feature_list_trainingset_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list_trainingset_cdna.md)
for the input,
[`plot_tail_features_trainingset_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/plot_tail_features_trainingset_cdna.md)
for visual inspection of single reads,
[`filter_nonA_chunks_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/filter_nonA_chunks_trainingset.md)
where the polarity is used.

## Examples

``` r
if (FALSE) { # \dontrun{

run_counts <- ninetails::count_pseudomove_runs_trainingset_cdna(
  tail_feature_list = tfl)
colSums(run_counts[, c("peak_runs", "valley_runs")])

} # }
```
