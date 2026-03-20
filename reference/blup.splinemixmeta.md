# Obtain predictions (BLUPs) from a `splinemixmeta` object

Obtain predictions (BLUPs) from a `splinemixmeta` object

## Usage

``` r
blup.splinemixmeta(
  object,
  se = FALSE,
  pi = FALSE,
  vcov = FALSE,
  pi.level = 0.95,
  type = "outcome",
  level,
  format,
  aggregate = "stat",
  ...
)
```

## Arguments

- object:

  An object of class `splinemixmeta`, returned by
  [`splinemixmeta()`](https://fawda123.github.io/splinemixmeta/reference/splinemixmeta.md).

- se:

  Logical indicating whether to return standard errors of the
  predictions.

- pi:

  Logical indicating whether to return prediction intervals.

- vcov:

  Logical indicating whether to return the variance-covariance matrix of
  the predictions.

- pi.level:

  Numeric value between 0 and 1 indicating the confidence level for the
  prediction intervals. Default is 0.95.

- type:

  Character string specifying the type of prediction: "outcome" for
  predicted outcomes or "residual" for predicted residuals.

- level:

  Integer indicating the random effects level for which to obtain
  predictions. Default is the highest level.

- format:

  Character string specifying the format of the output: "matrix" or
  "list". Default depends on the number of outcomes and whether `vcov`
  is requested.

- aggregate:

  Character string specifying how to aggregate the output when multiple
  outcomes are present

- ...:

  Additional arguments (currently unused).

## Value

A matrix or list of predicted values (BLUPs), optionally including
standard errors, prediction intervals, and variance-covariance matrices.

## Details

This function is modified from
[`mixmeta::blup.mixmeta`](https://rdrr.io/pkg/mixmeta/man/blup.mixmeta.html)
with acknowledgement of the original authors. It is modified to handle
intermediate levels of random effects more carefully. This is currently
a temporary solution that may be replaced in future versions.
