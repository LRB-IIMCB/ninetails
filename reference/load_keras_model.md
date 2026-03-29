# Load Keras model for multiclass signal prediction

Loads a pretrained convolutional neural network (CNN) model for
classifying Gramian Angular Field representations of poly(A) tail signal
chunks into nucleotide categories (A, C, G, U). When called without
arguments, the default model bundled with ninetails is loaded. A custom
model path can be provided for user-trained models.

## Usage

``` r
load_keras_model(keras_model_path)
```

## Arguments

- keras_model_path:

  Character string or missing. Full path of the `.h5` file with a Keras
  model. If missing (default), the pretrained model shipped with the
  package is loaded. Otherwise, the specified custom model is loaded.

## Value

A Keras model object ready for prediction.

## See also

[`predict_gaf_classes`](https://LRB-IIMCB.github.io/ninetails/reference/predict_gaf_classes.md)
where this function is called,
[`combine_gafs`](https://LRB-IIMCB.github.io/ninetails/reference/combine_gafs.md)
for preparing GAF input arrays

## Examples

``` r
if (FALSE) { # \dontrun{

# Load default pretrained model
model <- load_keras_model()

# Load custom model
model <- load_keras_model(
  keras_model_path = "/path/to/the/model/in/hdf5_format"
)

} # }
```
