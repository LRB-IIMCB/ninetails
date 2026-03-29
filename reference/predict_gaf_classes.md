# Classify Gramian Angular Field matrices with a pretrained CNN

Assigns GAF representations of signal chunks to one of four nucleotide
categories (A, C, G, U) using a pretrained convolutional neural network.
The model is loaded via
[`load_keras_model`](https://LRB-IIMCB.github.io/ninetails/reference/load_keras_model.md)
and inference is performed through the TensorFlow/Keras backend.

## Usage

``` r
predict_gaf_classes(gaf_list)
```

## Arguments

- gaf_list:

  List of GAF arrays (100, 100, 2), as produced by
  [`create_gaf_list`](https://LRB-IIMCB.github.io/ninetails/reference/create_gaf_list.md).

## Value

A named list with two elements:

- chunkname:

  Character vector. Names of the classified signal chunks (format:
  `<readname>_<index>`)

- prediction:

  Integer vector. Predicted nucleotide class for each chunk (0 = A, 1 =
  C, 2 = G, 3 = U)

## Details

The function reshapes the input GAF list into a 4-D tensor (n_chunks x
100 x 100 x 2) and applies the pretrained CNN to predict the nucleotide
class for each chunk. Predictions are returned as integer codes: 0 = A,
1 = C, 2 = G, 3 = U.

## See also

[`create_gaf_list`](https://LRB-IIMCB.github.io/ninetails/reference/create_gaf_list.md)
for preparing the input,
[`load_keras_model`](https://LRB-IIMCB.github.io/ninetails/reference/load_keras_model.md)
for model loading,
[`create_outputs`](https://LRB-IIMCB.github.io/ninetails/reference/create_outputs.md)
for integrating predictions into final output

## Examples

``` r
if (FALSE) { # \dontrun{

pl <- ninetails::predict_gaf_classes(gl)

} # }
```
