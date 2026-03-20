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

  xvar name of variable to plot on the horizontal axis.

- title:

  title to add to the plot

- xlab:

  xlab label for the horizontal axis

- ylab:

  ylab label for the vertical axis

- ylim:

  ylim limits for the vertical axis

- linecolor:

  color for the prediction line

- fillcolor:

  color for the prediction confidence band

- ...:

  additional arguments passed to
  [`predict.splinemixmeta()`](https://fawda123.github.io/splinemixmeta/reference/predict.splinemixmeta.md)

## Value

ggplot2 object

## Details

This is not a very general plotting function. It is intended to provide
a basic feature for visualizing a univariate `splinemixmeta` fit in a
way that:

- includes fixed effects and spline terms in the predicted values, with
  95% confidence bands

- Shows the data points with 95% confidence intervals obtains as +/- 2
  times the standard errors (`se` or `diag(S)`).

- returns a `ggplot2` object that can be further updated.
