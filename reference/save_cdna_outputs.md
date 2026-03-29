# Save cDNA pipeline outputs in standard ninetails format

This function takes the merged cDNA results and formats them as the
standard ninetails output format (same as create_outputs_dorado) with
two tables: read_classes and nonadenosine_residues. The key difference
from DRS is the additional tail_type column indicating polyA or polyT.

## Usage

``` r
save_cdna_outputs(outputs, save_dir, prefix = "")
```

## Arguments

- outputs:

  List. Merged results from merge_cdna_results().

- save_dir:

  Character string. Directory where outputs will be saved.

- prefix:

  Character string. Optional prefix for output file names.

## Value

List containing standard ninetails output format:

- read_classes:

  Data frame with read classifications including tail_type

- nonadenosine_residues:

  Data frame with modifications including tail_type

## Examples

``` r
if (FALSE) { # \dontrun{
final_results <- save_cdna_outputs(
  outputs = merged_results,
  save_dir = "path/to/output/",
  prefix = "experiment1"
)
} # }
```
