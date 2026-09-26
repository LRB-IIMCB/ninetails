# Draws the downsampled tail signal of a single training-set read with pseudomove runs highlighted.

Visual lookup for the cDNA training data. The shaded regions mark the
non-zero pseudomove runs (length \>= 4) which will become the centres of
the extracted chunks, labelled as peak (+1) or valley (-1). Inspect a
handful of reads per labelled dataset and orientation to confirm that
the deviations look like genuine residue-induced distortions and to
determine their polarity.

## Usage

``` r
plot_tail_features_trainingset_cdna(readname, tail_feature_list)
```

## Arguments

- readname:

  Character string. Name (UUID) of the given read within the feature
  list.

- tail_feature_list:

  List object produced by
  [`create_tail_feature_list_trainingset_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list_trainingset_cdna.md).

## Value

ggplot2 object with the tail signal and highlighted pseudomove runs.

## See also

[`count_pseudomove_runs_trainingset_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/count_pseudomove_runs_trainingset_cdna.md)
for the tabular lookup,
[`plot_tail_chunk`](https://LRB-IIMCB.github.io/ninetails/reference/plot_tail_chunk.md)
for plotting extracted chunks,
[`plot_gaf`](https://LRB-IIMCB.github.io/ninetails/reference/plot_gaf.md)
for plotting the resulting GAFs.

## Examples

``` r
if (FALSE) { # \dontrun{

example <- ninetails::plot_tail_features_trainingset_cdna(
  readname = "5c2386e6-32e9-4e15-a5c7-2831f4750b2b",
  tail_feature_list = tfl)

print(example)

} # }
```
