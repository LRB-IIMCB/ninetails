# Creates the training-set feature list (signal + pseudomoves) from Dorado cDNA tail signals of a single orientation.

This is the cDNA counterpart of
[`create_tail_feature_list_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list_trainingset.md)
(for C, G, U) and
[`create_tail_feature_list_A`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list_A.md)
(for A). It computes pseudomoves for every tail signal in parallel with
[`filter_signal_by_threshold_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/filter_signal_by_threshold_trainingset.md)
and then applies the nucleotide-specific retention criterion:

- `"C"`, `"G"`, `"U"`: keep reads whose pseudomove vector contains at
  least one non-zero run of length \>= 4 (potential modification
  present).

- `"A"`: keep reads whose pseudomove vector contains *no* non-zero run
  of length \>= 4 (pure homopolymer tail).

## Usage

``` r
create_tail_feature_list_trainingset_cdna(signal_list, num_cores, nucleotide)
```

## Arguments

- signal_list:

  Named list of numeric vectors. Tail signals of a single orientation,
  named by read ID (as produced by
  [`extract_tail_signals_trainingset_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/extract_tail_signals_trainingset_cdna.md)).

- num_cores:

  Numeric `[1]`. Number of physical cores to use. Do not exceed 1 less
  than the number of cores at your disposal.

- nucleotide:

  Character. One of `"A"`, `"C"`, `"G"` or `"U"`. The residue inserted
  into the tails of the analysed construct (known from the experimental
  design). Selects the retention criterion described above.

## Value

A named list with two elements:

- tail_feature_list:

  Named list of per-read feature lists with four slots: `pod5_filename`
  (`NA`), `tail_signal`, `tail_moves` (`NA`) and `tail_pseudomoves`.

- discarded_readnames:

  Character vector. Read IDs discarded by the nucleotide-specific
  criterion (no qualifying pseudomove run for C/G/U; a qualifying run
  present for A).

Always assign this returned list to a variable; printing the full list
to the console may crash the R session.

## Details

Dorado does not provide a move table, so the per-read feature sublists
carry `NA` in the `pod5_filename` and `tail_moves` slots. The four-slot
layout (`pod5_filename`, `tail_signal`, `tail_moves`,
`tail_pseudomoves`) is preserved on purpose, so that the Guppy
training-set chunkers
([`create_tail_chunk_list_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_chunk_list_trainingset.md),
[`create_tail_chunk_list_A`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_chunk_list_A.md))
and everything downstream can be reused without modification.
Consequently, no zero-moved read category exists in this variant.

Provide signals of a single orientation only (either `polya_signals` or
`polyt_signals` from
[`extract_tail_signals_trainingset_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/extract_tail_signals_trainingset_cdna.md));
mixing orientations in one feature list would produce a mixed training
set.

## See also

[`extract_tail_signals_trainingset_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/extract_tail_signals_trainingset_cdna.md)
for the preceding pipeline step,
[`count_pseudomove_runs_trainingset_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/count_pseudomove_runs_trainingset_cdna.md)
and
[`plot_tail_features_trainingset_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/plot_tail_features_trainingset_cdna.md)
for data inspection,
[`create_tail_chunk_list_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_chunk_list_trainingset.md)
and
[`create_tail_chunk_list_A`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_chunk_list_A.md)
for the next pipeline step,
[`prepare_trainingset_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/prepare_trainingset_cdna.md)
for the top-level wrapper.

## Examples

``` r
if (FALSE) { # \dontrun{

tfl <- ninetails::create_tail_feature_list_trainingset_cdna(
  signal_list = extracted$polya_signals,
  num_cores = 10,
  nucleotide = "C")

} # }
```
