# Fit a mixed-effects meta-analysis model with inclusion of one or more spline terms

Fit a mixed-effects meta-analysis model with inclusion of one or more
spline terms

## Usage

``` r
splinemixmeta(
  smooth = NULL,
  formula,
  se,
  S = se^2,
  manual_fixed = FALSE,
  residual_re = TRUE,
  data,
  random = list(),
  method = "reml",
  bscov = "unstr",
  ...
)
```

## Arguments

- smooth:

  A smoothing term created by
  [`mgcv::s()`](https://rdrr.io/pkg/mgcv/man/s.html), or a list of such
  terms.

- formula:

  A formula for the fixed effects part of the model.

- se:

  A vector of standard errors for the response variables.

- S:

  As an alternative to `se`, `S` can be provided in several formats to
  give variance-covariance information for the response variables.

- manual_fixed:

  If `TRUE`, the fixed-effect component (if any) of any `smooth` terms
  is being manually included in the `formula` and hence should not be
  obtained as the unpenalized component(s) of `smooth`. Normally the
  "fixed-effect component" is the linear component. Hence one should
  either provide `smooth = s(x)` with `x` *omitted* from `formula` (e.g.
  `formula = y` or `formula = y ~ 1`) and thus `manual_fixed = FALSE`,
  or provide `smooth = s(x)` with `x` *included* in `formula` (e.g.
  `formula = y ~ x`) and `manual_fixed = TRUE`. The model fits should be
  identical but the coefficient for `x` will differ because `x` will be
  scaled differently if it was automatically obtained from the spline
  basis functions (i.e. with `manual_fixed = FALSE`). In either case,
  the unpenalized dimensions of the smooth term and not include in the
  spline random effects.

- residual_re:

  If `TRUE`, a datum-level random effect for residual variation (beyond
  the measurement error specified by `se` or `S`) is automatically
  included (similar to the default behavior of
  [`mixmeta::mixmeta()`](https://rdrr.io/pkg/mixmeta/man/mixmeta.html)
  when `random` is not specified). Normally this should be `TRUE` unless
  there is a clear reason to set it `FALSE`.

- data:

  A data frame containing the variables in the model. If not provided,
  variables are sought from where the function was called.

- random:

  See
  [`mixmeta::mixmeta()`](https://rdrr.io/pkg/mixmeta/man/mixmeta.html).
  This is a list of one-sided formulas specifying additional random
  effects beyond those that will be created from the `smooth` argument
  and `residual_re`.

- method:

  This *must* be `reml`. It is provided as an argument to make clear
  that only `reml` is supported for estimating models where spline
  formulations are set up as random effects. This simplifies catching
  cases where a user might try to pass a different `method` value to
  [`mixmeta::mixmeta()`](https://rdrr.io/pkg/mixmeta/man/mixmeta.html)
  via `...`.

- bscov:

  See
  [`mixmeta::mixmeta()`](https://rdrr.io/pkg/mixmeta/man/mixmeta.html).
  This is relevant only if `random` is provided.

- ...:

  Additional arguments passed to
  [`mixmeta::mixmeta()`](https://rdrr.io/pkg/mixmeta/man/mixmeta.html).'

## Value

An object of class `splinemixmeta`, unless there are no `smooth` terms,
in which case an object of class `mixmeta` is returned.

## Details

This function combines capabilities of `mgcv` and `mixmeta` in order to
provide spline meta-regression, which means a meta-regression model
where the shape of the relationship is unspecified and estimated from
the data with smoothing splines. Spline components built from `mgcv` can
be represented as random effects (along with fixed effects, which are
unpenalized, typically for linear terms). `mixmeta` supports fairly
flexible specification of fixed and random effects in (univariate or
multivariate) meta-analysis regression (meta-regression) models.
`splinemixmeta` takes `mgcv`-style specifications of smooth (spline)
terms, sets them up for `mixmeta`, and the calls `mixmeta` to fit the
model by REML.

Only a limited set of `s` arguments and options are supported. Argument
`k` should work. For `bs`, supported basis functions include "cr", "cs",
and "cc". Note that the default choice `bs = "tp"` does not work well
and so results in a warning. (The supported basis functions are those
for which
[`mgcv::smoothCon`](https://rdrr.io/pkg/mgcv/man/smoothCon.html) can
produce diagonal penalty matrices are supported ("cr", "cs", "cc"), with
the exception of "tp". Arguments `fx`, `m`, `by`, `id`, and `sp` are not
supported should not be provided. Argument `xt` should work but is
untested. Argument `pc` is untested.

Note that `bs="cc"` (cyclic cubic regression spline) does not have
unpenalized components, so if this is used, `manual_fixed` is not
relevant.

`splinemixmeta` is not particularly optimized for large data sets.

## See also

- [`predict.splinemixmeta()`](https://perrydv.github.io/splinemixmeta/reference/predict.splinemixmeta.md)
  for predictions based on BLUPs (best linear unbiased predictors) from
  fitted `splinemixmeta` models. This is an S3 method that will be
  called from `predict(x)` where `x` is a `splinemixmeta` object.

- [`plot.splinemixmeta()`](https://perrydv.github.io/splinemixmeta/reference/plot.splinemixmeta.md)
  for plotting fitted spline meta-regression models. This is an S3
  method that will be called from `plot(x)` where `x` is a
  `splinemixmeta` object.

- [`blup.splinemixmeta()`](https://perrydv.github.io/splinemixmeta/reference/blup.splinemixmeta.md)
  for obtaining BLUPs (best linear unbiased predictors) from fitted
  `splinemixmeta` models. This is used by `predict.splinemixmeta`, which
  is typically easier to call directly. This is an S3 method that will
  be called from `blup(x)` where `x` is a `splinemixmeta` object.

- [`make_smm_smooth()`](https://perrydv.github.io/splinemixmeta/reference/make_smm_smooth.md)
  for the internal function that sets up spline terms for use in
  `splinemixmeta`.
