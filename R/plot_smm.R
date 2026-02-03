#' Plot results from a univariate `splinemixmeta` model
#'
#' @param object object of class `mixmeta` returned from `splinemixmeta`
#' @param xvar xvar name of variable to plot on the horizontal axis.
#' @param title title to add to the plot
#' @param xlab xlab label for the horizontal axis
#' @param ylab ylab label for the vertical axis
#' @param ylim ylim limits for the vertical axis
#' @param linecolor color for the prediction line
#' @param fillcolor color for the prediction confidence band
#'
#' @details
#' This is not a very general plotting function. It is intended to provide a basic feature
#' for visualizing a univariate `splinemixmeta` fit in a way that:
#'
#' - includes fixed effects and spline terms in the predicted values, with 95% confidence bands
#' - Shows the data points with 95% confidence intervals obtains as +/- 2 times the standard errors (`se` or `diag(S)`).
#' - returns a `ggplot2` object that can be further updated.
#'
#' @returns ggplot2 object
#' @export
plot.splinemixmeta <- function(object, xvar, title, xlab, ylab, ylim, linecolor = "blue", fillcolor = "blue", ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop(
      "Package 'ggplot2' is required to use `plot.splinemixmeta\n",
      call. = FALSE
    )
  if(!inherits(object, "splinemixmeta"))
    stop("object must be of class 'splinemixmeta'")
  if(missing(xvar)) {
    xvar_expr <- object$formula[[3]]
    if(!is.name(xvar_expr))
      stop("Please provide xvar. A default xvar could not be chosen.")
  } else {
    xvar_expr <- substitute(xvar)
  }
  xvar <- as.character(xvar_expr)
  if(missing(xlab)) xlab <- xvar
  pred <- predict.splinemixmeta(object, ...)
  xdata <- object$model[[xvar]]
  order <- order(xdata)
  yvar_expr <- object$formula[[2]]
  if(missing(ylab)) ylab <- as.character(yvar_expr)
  ydata <- object$model[[as.character(yvar_expr)]]
  xdata <- xdata[order]
  ydata <- ydata[order]
  pred <- pred[order, , drop=FALSE]
  se <- (if(is.matrix(object$S)) diag(object$S) else object$S) |> sqrt()
  se <- se[order]
  df <- data.frame(x_plot_ = xdata,
                   y_plot_ = ydata,
                   se_plot_ = se,
                   mmpred_ = pred[,'blup'],
                   mmse_ = pred[,'se'])
  df <- cbind(df, object$model[order, , drop=FALSE])
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
