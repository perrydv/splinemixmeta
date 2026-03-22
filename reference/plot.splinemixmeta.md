# Plot results from a univariate `splinemixmeta` model

Plot results from a univariate `splinemixmeta` model

## Usage

``` r
# S3 method for class 'splinemixmeta'
plot(
  x,
  xvar,
  title,
  xlab,
  ylab,
  ylim,
  linecolor = "blue",
  fillcolor = "blue",
  ...
)
```

## Arguments

- x:

  object of class `mixmeta` returned from `splinemixmeta`

- xvar:

  name of the variable to plot on the horizontal axis. If missing, this
  will will be taken from the right side of the fixed effects formula.

- title:

  to add to the plot

- xlab:

  label for the horizontal axis

- ylab:

  label for the vertical axis

- ylim:

  limits for the vertical axis

- linecolor:

  color for the prediction line

- fillcolor:

  color for the prediction confidence band

- ...:

  additional arguments passed to
  [`predict.splinemixmeta()`](https://perrydv.github.io/splinemixmeta/reference/predict.splinemixmeta.md)

## Value

`ggplot2` object

## Details

This is not a very general plotting function. It is intended to provide
a basic feature for visualizing a univariate `splinemixmeta` fit in a
way that:

- includes fixed effects and spline terms in the predicted values, with
  95% confidence bands

- Shows the data points with 95% confidence intervals obtained as +/- 2
  times the standard errors (`se` or `diag(S)`).

- returns a `ggplot2` object that can be further updated.

If the `x` object comes from a call to
[`splinemixmeta()`](https://perrydv.github.io/splinemixmeta/reference/splinemixmeta.md)
that had a simple (fixed effects) `formula`, such as `y ~ w`, then
`xvar` does not need to be provided, since it is easily determined to be
`w`. However, if the formula was more complicated, or if `w` was not
provided in `formula` because `manual_fixed = FALSE` (the default), then
`xvar = "w"` must be provided so that this plotting function knows what
to put on the horizontal axis.
