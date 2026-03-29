# Creates a visual representation of multiple gramian angular fields based on provided gaf_list (plots all gafs from the given list).

The function saves the plots to files with predefined size 100x100
pixels.

## Usage

``` r
plot_multiple_gaf(gaf_list, num_cores)
```

## Arguments

- gaf_list:

  A list of gaf matrices organized by the read ID_index.

- num_cores:

  numeric \[1\]. Number of physical cores to use in processing the data.
  Do not exceed 1 less than the number of cores at your disposal.

## Value

multiple png files containing gramian angular fields representing given
fragment of nanopore read signal.

## Details

A feature useful for creating custom data sets for network training. It
is recommended to use this feature with caution. Producing multiple
graphs from extensive data sets may cause the system to crash.

This function plots one-channel (100,100,1) as well as multi-channel
gafs (e.g. 100,100,2).

IMPORTANT NOTE! In current version, this function plots multi-channel
matrices in collapsed manner. If one wants to separate color spaces to
GASF/GADF channels or to R,G,B space, this function would not be
suitable. The data would require additional processing steps!

## Examples

``` r
if (FALSE) { # \dontrun{

ninetails::plot_multiple_gaf(gaf_list = gl, num_cores = 10)

} # }
```
