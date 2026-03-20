#' Plot results from a univariate `splinemixmeta` model
#'
#' @param x object of class `mixmeta` returned from `splinemixmeta`
#' @param xvar name of the variable to plot on the horizontal axis. If missing, this will
#' will be taken from the right side of the fixed effects formula.
#' @param title to add to the plot
#' @param xlab label for the horizontal axis
#' @param ylab label for the vertical axis
#' @param ylim limits for the vertical axis
#' @param linecolor color for the prediction line
#' @param fillcolor color for the prediction confidence band
#' @param ... additional arguments passed to `predict.splinemixmeta()`
#'
#' @details
#' This is not a very general plotting function. It is intended to provide a basic feature
#' for visualizing a univariate `splinemixmeta` fit in a way that:
#'
#' - includes fixed effects and spline terms in the predicted values, with 95% confidence bands
#' - Shows the data points with 95% confidence intervals obtained as +/- 2 times the standard errors (`se` or `diag(S)`).
#' - returns a `ggplot2` object that can be further updated.
#'
#' If the `x` object comes from a call to [splinemixmeta()] that had a simple (fixed effects) `formula`, such as `y ~ w`, then
#' `xvar` does not need to be provided, since it is easily determined to be `w`. However, if the
#' formula was more complicated, or if `w` was not provided in `formula` because `manual_fixed = FALSE` (the default),
#' then `xvar = "w"` must be provided so that this plotting function knows what to put on the horizontal axis.
#'
#' @method plot splinemixmeta
#'
#' @returns `ggplot2` object
#' @export
plot.splinemixmeta <- function(x, xvar, title, xlab, ylab, ylim, linecolor = "blue", fillcolor = "blue", ...) {

  if(!inherits(x, "splinemixmeta"))
    stop("x must be of class 'splinemixmeta'")
  if(missing(xvar)) {
    xvar_expr <- x$formula[[3]]
    if(!is.name(xvar_expr))
      stop("Please provide xvar. A default xvar could not be chosen.")
  } else {
    xvar_expr <- substitute(xvar)
  }
  xvar <- as.character(xvar_expr)
  if(missing(xlab)) xlab <- xvar
  pred <- predict.splinemixmeta(x, ...)
  xdata <- x$model[[xvar]]
  order <- order(xdata)
  yvar_expr <- x$formula[[2]]
  if(missing(ylab)) ylab <- as.character(yvar_expr)
  ydata <- x$model[[as.character(yvar_expr)]]
  xdata <- xdata[order]
  ydata <- ydata[order]
  pred <- pred[order, , drop=FALSE]
  se <- (if(is.matrix(x$S)) diag(x$S) else x$S) |> sqrt()
  se <- se[order]
  df <- data.frame(x_plot_ = xdata,
                   y_plot_ = ydata,
                   se_plot_ = se,
                   mmpred_ = pred[,'blup'],
                   mmse_ = pred[,'se'])
  df <- cbind(df, x$model[order, , drop=FALSE])
  aesmm <- ggplot2::aes(x = x_plot_, y = mmpred_)
  fig <- ggplot2::ggplot(df, ggplot2::aes(x = x_plot_, y = y_plot_)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = y_plot_ - 2*se_plot_, ymax = y_plot_ + 2*se_plot_)) +
    ggplot2::geom_point() +
    ggplot2::geom_ribbon(ggplot2::aes(x = x_plot_,
                    ymin = mmpred_ - 2*mmse_,
                    ymax = mmpred_ + 2*mmse_),
                data = df,
                inherit.aes = FALSE,
                fill = fillcolor,alpha = 0.25) +
    ggplot2::geom_line(mapping = aesmm, data = df, inherit.aes=FALSE, color = linecolor,linewidth=1)

  fig <- fig + ggplot2::xlab(xlab)
  fig <- fig + ggplot2::ylab(ylab)
  if(!missing(title)) fig <- fig + ggplot2::ggtitle(title)
  if(!missing(ylim)) fig <- fig + ggplot2::coord_cartesian(ylim = ylim)
  fig
}
