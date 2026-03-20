# Make predictions from a fitted `splinemixmeta` model

Make predictions from a fitted `splinemixmeta` model

## Usage

``` r
# S3 method for class 'splinemixmeta'
predict(
  object,
  include_smooths = TRUE,
  include_REs = FALSE,
  include_residuals = FALSE,
  type = "outcome",
  ...
)
```

## Arguments

- object:

  A fitted `mixmeta` object returned from `splinemixmeta`

- include_smooths:

  `TRUE` to include the smooth (spline) terms in predictions. Typically
  one wants these.

- include_REs:

  `TRUE` to include any additional random effects (beyond the smooths)
  that were provided in the `random` argument to `splinemixmeta`. Omit
  these if you want to see just the spline predictions.

- include_residuals:

  TRUE to include the random effects (one for each datum) for residual
  variation not accounted for in the measurement variation (S). Only
  matters if `residual_re = TRUE` when calling `splinemixmeta`, which
  should typically be the case. Typically one does not want these in
  predictions.

- type:

  Type of predictions. This can be "outcome" or "residual" and will be
  passed to the `type` argument of
  [`mixmeta::blup()`](https://rdrr.io/pkg/mixmeta/man/blup.html).

- ...:

  Additional arguments (currently unused)

## Value

A matrix with columns "blup" for the predicted values, "se" for the
standard errors of the predictions, and "vcov" for the variance of the
predictions. These are returned from
[`mixmeta::blup()`](https://rdrr.io/pkg/mixmeta/man/blup.html) with
`vcov=TRUE` and `se=TRUE`.

## Details

This is a convenience function that calls
[`mixmeta::blup`](https://rdrr.io/pkg/mixmeta/man/blup.html) without
requiring you to know which random-effects "levels" of the fitted
mixmeta object correspond to which parts of the model. For more
fine-grained control (such as including one spline term but not
another), one can use
[`mixmeta::blup()`](https://rdrr.io/pkg/mixmeta/man/blup.html) directly.
