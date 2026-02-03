#' Create random effects and fixed effects specifications from an mgcv smooth term for use in `splinemixmeta`.
#'
#' This function is for internal use by `splinemixmeta`.
#'
#' @param smooth A term created by `mgcv::s()`
#' @param data Data frame containing the variables used in the smooth
#' @param vnames "a vector of names to avoid as dummy variable names in the random effects form",
#'  per `help(mgcv::smooth2random)`, to which `vnames` is passed.
#' @param manual_fixed If `FALSE`, the unpenalized dimensions of the smooth are extracted for use
#'  as fixed effects. If `TRUE`, the user should provide any desired fixed effects directly. See details.
#' @param envir The environment in which to evaluate variable names if `data` is not provided.
#'
#' @returns A list with elements `basisFxns` (the basis functions for the random effects)
#'         and `x_fixed` (the fixed effects design matrix, or `NULL` if
#'         `manual_fixed` is `TRUE`).
#'
#' @details
#'
#' This function uses `mgcv::smoothCon()` and `mgcv::smooth2random()` to obtain basis functions, penalty matrix, and
#'  (optionally) fixed effects terms from the `smooth` specification.
#'
#'  The fixed effects (if `manual_fixed` is `FALSE`) represent unpenalized directions of the smooth. Typically, this
#'  means the fixed effects will include linear terms, because splines usually penalize curvature, so any parameters
#'  that give a line are unpenalized. If one is not particularly interested in the linear terms (more generally, unpenalized terms),
#'  then the default of `manual_fixed = FALSE` is a good option. However, if one is interested in coefficients for the
#'  linear terms, it is important to note that when they are extracted from the basis function setup, they may be (typically will be)
#'  also re-scaled (and it is not particularly easy to determine the scaling factor). Hence, one may prefer to set `manual_fixed=FALSE`
#'  and provide the linear term directly in the `formula` argument to `splinemixmeta()`.
#'
#' @export
#' @importFrom mgcv smoothCon smooth2random
make_smm_smooth <- function(smooth, data, vnames = character(), manual_fixed = FALSE, envir=parent.frame()) {
  if(missing(data)) {
    data <- smooth$term |> lapply(get, envir = envir) |>
      setNames(smooth$term) |> as.data.frame()
  } else {
    # provided data can be list or data.frame,
    # but if it is a list, smoothCon needs and n element
    # and if there are multiple terms in the smooth, only a data.frame is supported
    if(length(smooth$term)==1) {
      data_for_smooths <- data
      if(is.list(data_for_smooths))
        data_for_smooths$n <- length(data[[smooth$term]])
    } else {
      if(is.list(data))
        if(!is.null(data$n))
          stop("When a smooth term involves multiple variables, do not include n in the data.")
      okay <- try(data_for_smooths <- as.data.frame(data))
      if(inherits(okay, "try-error")) {
        stop("When a smooth term involves multiple variables, the 'data' argument must be a data frame, or something that `as.data.frame` will convert to a data frame.")
      }
    }
  }

  sCdiag <- mgcv::smoothCon(smooth, data = data,
                            absorb.cons = TRUE,
                            diagonal.penalty = TRUE)
  s2r_list <- sCdiag |> lapply(mgcv::smooth2random, vnames = vnames, type = 2)
  make_result <- \(s2r) {
    basisFxns <- s2r$rand$Xr
    x_fixed <- NULL
    if(!manual_fixed) {
      x_fixed <- s2r$Xf
    }
    list(basisFxns = basisFxns,
         x_fixed = x_fixed)
  }
  results <- lapply(s2r_list, make_result)
  results
}
