# Draws a portion of poly(A) tail squiggle (chunk) for given read.

This function allows to visualise a single fragment of the poly(A) tail
area defined by the fragment name (chunk_name). Intended for use on the
output of create_tail_chunk_list_moved function.

## Usage

``` r
plot_tail_chunk(chunk_name, tail_chunk_list)
```

## Arguments

- chunk_name:

  character string. Name of the given read segment (chunk) within the
  given tail_chunk_list object.

- tail_chunk_list:

  character string. The list object produced by create_chunk_list
  function.

## Value

ggplot2 object with squiggle plot depicting given fragment of analyzed
poly(A) tail.

## Details

Currently, draws only the raw signal, without the option to scale to
picoamperes \[pA\].

## Examples

``` r
if (FALSE) { # \dontrun{

example <- ninetails::plot_tail_chunk(
 chunk_name = "5c2386e6-32e9-4e15-a5c7-2831f4750b2b_1",
 tail_chunk_list = tcl)

print(example)

} # }
```
