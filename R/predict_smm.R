
#' Title
#'
#' @param object A fitted `mixmeta` object returned from `splinemixmeta`
#' @param include_smooths TRUE to include the smooth (spline) terms in predictions
#' @param include_REs TRUE to include any additional random effects (beyond the smooths) that
#'   were provided in the `random` argument to `splinemixmeta`
#' @param include_residuals TRUE to include the random effects (one for each datum) for residual
#'   variation not accounted for in the measurement variation (S). Only matters if `residual_re = TRUE`
#'   when calling `splinemixmeta`, which should typically be the case.
#'
#' @returns A matrix with columns "blup" for the predicted values, "se" for the standard errors of the predictions,
#'   and "vcov" for the variance of the predictions. These are returns from `mixmeta::blup()`.
#' @export
predict_smm <- function(object, include_smooths = TRUE,
                        include_REs = FALSE, include_residuals = FALSE) {
  if(!inherits(object, "mixmeta"))
    stop("object must be of class 'mixmeta'")
  level <- 0
  if(include_smooths) {
    level <- level + 1
  }
  if(include_REs) {
    level <- level + 1
  }
  if(include_residuals) {
    level <- level + 1
  }
  pred <- mixmeta::blup(object, vcov = TRUE, se = TRUE, level = level, type = "outcome")
  pred
}
