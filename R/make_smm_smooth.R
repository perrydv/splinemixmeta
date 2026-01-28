#' Create random effects and fixed effects specifications from an mgcv smooth term
#'
#' This function is for internal use.
#'
#' @param smooth A term created by `mgcv::s()`
#' @param data Data frame containing the variables used in the smooth
#' @param vnames "a vector of names to avoid as dummy variable names in the random effects form", per `help(mgcv::smooth2random)`, to which `vnames` is passed.
#' @param auto_fixed_effects If TRUE, the unpenalized dimensions of the smooth are extracted for use as fixed effects.
#'
#' @returns A list with elements `basisFxns` (the basis functions for the random effects)
#'         and `x_fixed` (the fixed effects design matrix, or NULL if
#'         `auto_fixed_effects` is FALSE)
#'
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
