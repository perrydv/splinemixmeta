#' Make predictions from a fitted `splinemixmeta` model
#'
#' @param object A fitted `mixmeta` object returned from `splinemixmeta`
#' @param include_smooths `TRUE` to include the smooth (spline) terms in predictions. Typically one wants these.
#' @param include_REs `TRUE` to include any additional random effects (beyond the smooths) that
#'   were provided in the `random` argument to `splinemixmeta`. Omit these if you want to see just the spline predictions.
#' @param include_residuals TRUE to include the random effects (one for each datum) for residual
#'   variation not accounted for in the measurement variation (S). Only matters if `residual_re = TRUE`
#'   when calling `splinemixmeta`, which should typically be the case. Typically one does not want these in predictions.
#' @param type Type of predictions. This can be "outcome" or "residual" and will be passed to the `type` argument of `mixmeta::blup()`.
#' @param ... Additional arguments (currently unused)
#' 
#' @returns A matrix with columns "blup" for the predicted values, "se" for the standard errors of the predictions,
#'   and "vcov" for the variance of the predictions. These are returned from `mixmeta::blup()` with `vcov=TRUE` and `se=TRUE`.
#'
#' @details This is a convenience function that calls `mixmeta::blup` without requiring you to know which random-effects "levels"
#'  of the fitted mixmeta object correspond to which parts of the model. For more fine-grained control (such as including one spline
#'  term but not another), one can use `mixmeta::blup()` directly.
#'
#' @method predict splinemixmeta
#' 
#' @export
predict.splinemixmeta <- function(object, include_smooths = TRUE,
                        include_REs = FALSE, include_residuals = FALSE, type = "outcome", ...) {
  if(!inherits(object, "splinemixmeta"))
    stop("object must be of class 'splinemixmeta'")
  level <- 0
  if(include_smooths) {
    level <- level + object$num_smooth
  }
  if(include_REs) {
    level <- level + object$num_random
  }
  if(include_residuals) {
    level <- level + object$num_residual
  }
  pred <- mixmeta::blup(object, vcov = TRUE, se = TRUE, level = level, type = type)
  pred
}
