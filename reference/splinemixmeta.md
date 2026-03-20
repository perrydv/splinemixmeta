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
  extracted from the `smooth`. Normally the "fixed-effect component" is
  the linear component. Hence one should either provide `smooth = s(x)`
  with `x` *omitted* from `formula` (e.g. `formula = y` or
  `formula = y ~ 1`) and thus `manual_fixed = FALSE`, or provide
  `smooth = s(x)` with `x` *included* in `formula` (e.g.
  `formula = y ~ x`) and `manual_fixed = TRUE`. The model fits should be
  identical but the coefficient for `x` will differ because `x` will be
  scaled differently if it was automatically extracted from the spline
  basis functions (i.e. with `manual_fixed = FALSE`).

- residual_re:

  If `TRUE`, a datum-level random effect for residual variation (beyond
  the measurement error specified by `se` or `S`) is automatically
  included (similar to the default behavior of
  [`mixmeta::mixmeta()`](https://rdrr.io/pkg/mixmeta/man/mixmeta.html)
  when `random` is not specified). Normally this should be `TRUE`.

- data:

  A data frame containing the variables in the model. If not provided,
  variables are sought from the where the function was called.

- random:

  See
  [`mixmeta::mixmeta()`](https://rdrr.io/pkg/mixmeta/man/mixmeta.html).
  This is a list of one-sided formulas specifying additional random
  effects beyond those that will be created from the `smooth` argument.

- method:

  This *must* be `reml`. It is provided as an argument to make clear
  that only `reml` is supported for estimating models where spline
  formulations are set up as random effects. This simplifies catching
  cases where a user might have tried to pass a different `method` value
  to
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

Only a limited set of `s` options are supported. `k` should work. For
`bs`, only basis functions for which `smoothCon` can produce diagonal
penalty matrices are supported (e.g. "tp", "cr", "cs", "cc"). `fx`, `m`,
`by`, `id`, and `sp` are not supported should not be provided. `xt`
should work but is untested. `pc` is untested.

Note that `bs="cc"` (cyclic cubic regression spline) does not have
unpenalized components, so if this is used, `manual_fixed` is not
relevant.

`splinemixmeta` is not particularly optimized for large data sets.

## See also

- [`predict.splinemixmeta()`](https://fawda123.github.io/splinemixmeta/reference/predict.splinemixmeta.md)
  for predictions based on BLUPs (best linear unbiased predictors) from
  [`mixmeta::blup.mixmeta()`](https://rdrr.io/pkg/mixmeta/man/blup.mixmeta.html)
  from fitted `splinemixmeta` models

- [`plot.splinemixmeta()`](https://fawda123.github.io/splinemixmeta/reference/plot.splinemixmeta.md)
  for plotting fitted spline meta-regression models

- [`make_smm_smooth()`](https://fawda123.github.io/splinemixmeta/reference/make_smm_smooth.md)
  for the internal function that sets up spline terms for use in
  `splinemixmeta`.
