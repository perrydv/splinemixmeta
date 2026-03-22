# Create random effects and fixed effects specifications from an mgcv smooth term for use in `splinemixmeta`.

This function is for internal use by
[`splinemixmeta()`](https://perrydv.github.io/splinemixmeta/reference/splinemixmeta.md).

## Usage

``` r
make_smm_smooth(
  smooth,
  data,
  vnames = character(),
  manual_fixed = FALSE,
  envir = parent.frame()
)
```

## Arguments

- smooth:

  A term created by [`mgcv::s()`](https://rdrr.io/pkg/mgcv/man/s.html)

- data:

  Data frame containing the variables used in the smooth. See `envir`.

- vnames:

  "a vector of names to avoid as dummy variable names in the random
  effects form", per
  [`mgcv::smooth2random()`](https://rdrr.io/pkg/mgcv/man/smooth2random.html),
  to which `vnames` is passed.

- manual_fixed:

  If `FALSE`, the unpenalized dimensions (typically a linear term) of
  the smooth are used as fixed effects. If `TRUE`, the user should
  provide any desired fixed effects directly. In either case,
  unpenalized dimensions of the smooth term are not included in the
  spline. See details.

- envir:

  The environment in which to evaluate variable names if `data` is not
  provided.

## Value

A list with elements `basisFxns` (the basis functions for the random
effects) and `x_fixed` (the fixed effects design matrix, or `NULL` if
`manual_fixed` is `TRUE`).

## Details

This function uses
[`mgcv::smoothCon()`](https://rdrr.io/pkg/mgcv/man/smoothCon.html) and
[`mgcv::smooth2random()`](https://rdrr.io/pkg/mgcv/man/smooth2random.html)
to obtain basis functions, penalty matrix, and (optionally) fixed
effects terms from the `smooth` specification.

The fixed effects (if `manual_fixed` is `FALSE`) represent unpenalized
directions of the smooth. Typically, this means the fixed effects will
include linear terms, because splines usually penalize curvature, so any
parameters that give a line are unpenalized. If one is not particularly
interested in the linear terms (more generally, unpenalized terms), then
the default of `manual_fixed = FALSE` is a good option. However, if one
is interested in coefficients for the linear terms, it is important to
note that when they are extracted from the basis function setup, they
may be (typically will be) also re-scaled, and it is not particularly
easy to determine the scaling factor. Hence, one may prefer to set
`manual_fixed=FALSE` and provide the linear term directly in the
`formula` argument to
[`splinemixmeta()`](https://perrydv.github.io/splinemixmeta/reference/splinemixmeta.md).
