# Extracts fragments of poly(A) tail signal (Dorado mode) containing potential modifications along with their delimitation (positional indices; coordinates) within the tail.

This function processes raw poly(A) tail signal and Dorado-derived
pseudomoves to identify and extract signal segments (chunks) potentially
corresponding to modified positions (e.g., non-A residues). Each
extracted chunk spans 100 signal points, centered on the midpoint of a
pseudomove run.

## Usage

``` r
split_tail_centered_dorado(readname, tail_feature_list)
```

## Arguments

- readname:

  character string. Name of the given read within the analyzed dataset.

- tail_feature_list:

  list object produced by
  [`create_tail_features_list_dorado`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_features_list_dorado.md)
  Dorado-tail feature extraction function. Must contain `$tail_signal`
  (numeric vector) and `$tail_pseudomoves` (integer vector) for each
  read.

## Value

a nested list where each element corresponds to a signal fragment. Each
fragment is itself a list with three entries:

- `chunk_sequence`: numeric vector of raw signal values

- `chunk_start_pos`: integer, start index of the chunk

- `chunk_end_pos`: integer, end index of the chunk

## Details

In the Dorado pipeline, moves are not used: \* retrieving them from BAM
files is computationally expensive \* processing is non-intuitive

Instead, only pseudomoves are considered. As a safeguard against
Doradoâ€™s tendency to extend poly(A) boundaries into the transcript
body, the last 3 pseudomove values are forced to 0. This prevents
misclassification of transcript nucleotides as part of the tail.

Candidate modification regions are detected by: \* run-length encoding
(RLE) of the pseudomove vector \* filtering runs of pseudomoves with
length â‰¥ 5

Extracted fragments are padded/imputed if they extend beyond signal
boundaries: \* upstream/downstream missing values (NAs) are replaced \*
imputation is based on random draws from the 5 most frequent signal
values

The function returns a list object (nested), where each element
represents one candidate modification region, containing: \*
\`chunk_sequence\`: the raw signal subsequence (length = 100, imputed if
needed) \* \`chunk_start_pos\`: starting index of the subsequence \*
\`chunk_end_pos\`: ending index of the subsequence

## Examples

``` r
if (FALSE) { # \dontrun{
ninetails::split_tail_centered_dorado(
  readname = "1234-anexample-r3adn4m3",
  tail_feature_list = tail_feature_list
)
} # }
```
