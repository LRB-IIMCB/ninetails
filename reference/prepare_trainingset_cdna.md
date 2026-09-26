# Produces GAF training data of a given nucleotide and read orientation from Dorado cDNA data.

Top-level convenience wrapper mirroring
[`prepare_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/prepare_trainingset.md)
for the cDNA pipeline. It chains
[`extract_tail_signals_trainingset_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/extract_tail_signals_trainingset_cdna.md),
[`create_tail_feature_list_trainingset_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list_trainingset_cdna.md)
and the reused Guppy training-set chunking, filtering and GAF functions
for a single nucleotide and a single orientation.

## Usage

``` r
prepare_trainingset_cdna(
  nucleotide,
  tail_type,
  dorado_summary,
  bam_file,
  pod5_dir,
  num_cores = 1,
  contig = NA,
  value = NA
)
```

## Arguments

- nucleotide:

  Character. One of `"A"`, `"C"`, `"G"` or `"U"`. The residue inserted
  into the tails of the analysed construct (known from the experimental
  design).

- tail_type:

  Character. Either `"polyA"` or `"polyT"`. Read orientation for which
  the training data are produced.

- dorado_summary:

  Character string or data frame. Dorado summary (see
  [`extract_tail_signals_trainingset_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/extract_tail_signals_trainingset_cdna.md)).

- bam_file:

  Character string. Full path of the aligned BAM file.

- pod5_dir:

  Character string. Full path of the directory containing POD5 files.

- num_cores:

  Numeric `[1]`. Number of physical cores to use. Do not exceed 1 less
  than the number of cores at your disposal.

- contig:

  Character vector `[NA]`. Reference contigs to retain (see
  [`extract_tail_signals_trainingset_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/extract_tail_signals_trainingset_cdna.md)).

- value:

  Numeric `[NA]`. Pseudomove polarity of the residue of interest in the
  given orientation: `1` to retain chunks with peaks, `-1` to retain
  chunks with valleys. Required for `"C"`, `"G"` and `"U"`; ignored for
  `"A"`.

## Value

A named list of GAF arrays (100, 100, 2) organised by
`<read_ID>_<index>`. Always assign this returned list to a variable;
printing the full list to the console may crash the R session.

## Details

The internal pipeline differs by nucleotide:

- `"A"`:

  [`create_tail_chunk_list_A`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_chunk_list_A.md)
  \\\rightarrow\\
  [`create_gaf_list_A`](https://LRB-IIMCB.github.io/ninetails/reference/create_gaf_list_A.md)
  (overlapping windows, data augmentation).

- `"C"`, `"G"`, `"U"`:

  [`create_tail_chunk_list_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_chunk_list_trainingset.md)
  \\\rightarrow\\
  [`filter_nonA_chunks_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/filter_nonA_chunks_trainingset.md)
  (with the supplied `value`) \\\rightarrow\\
  [`create_gaf_list`](https://LRB-IIMCB.github.io/ninetails/reference/create_gaf_list.md).

Unlike the Guppy/DRS routine, the pseudomove polarity of each residue is
not hardcoded, because it has not been established for the DNA chemistry
and for both orientations. Determine it first with
[`count_pseudomove_runs_trainingset_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/count_pseudomove_runs_trainingset_cdna.md)
and
[`plot_tail_features_trainingset_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/plot_tail_features_trainingset_cdna.md),
then pass it as `value`. Run the wrapper once per nucleotide and per
orientation; the polyA and polyT sets train two separate models.

## See also

[`prepare_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/prepare_trainingset.md)
for the Guppy counterpart,
[`extract_tail_signals_trainingset_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/extract_tail_signals_trainingset_cdna.md),
[`create_tail_feature_list_trainingset_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list_trainingset_cdna.md),
[`count_pseudomove_runs_trainingset_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/count_pseudomove_runs_trainingset_cdna.md),
[`filter_nonA_chunks_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/filter_nonA_chunks_trainingset.md),
[`create_gaf_list`](https://LRB-IIMCB.github.io/ninetails/reference/create_gaf_list.md),
[`create_gaf_list_A`](https://LRB-IIMCB.github.io/ninetails/reference/create_gaf_list_A.md).

## Examples

``` r
if (FALSE) { # \dontrun{

# C-containing construct, polyT orientation, valleys established beforehand
C_polyt_gafs <- ninetails::prepare_trainingset_cdna(
  nucleotide = "C",
  tail_type = "polyT",
  dorado_summary = '/path/to/dorado_summary.txt',
  bam_file = '/path/to/aligned.bam',
  pod5_dir = '/path/to/pod5_dir/',
  num_cores = 10,
  contig = "RlucB",
  value = -1)

} # }
```
