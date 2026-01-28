#' Title
#'
#' @param object mmfit
#' @param xvar xvar
#' @param title title
#' @param xlab xlab
#' @param ylab ylab
#' @param ylim ylim
#' @param color color
#'
#' @returns ggplot2 object
#' @export
smm_plot <- function(object, xvar, title, xlab, ylab, ylim, color = "blue") {
  browser()
  if(!inherits(object, "mixmeta"))
    stop("object must be of class 'mixmeta'")
  if(missing(xvar)) {
    xvar_expr <- object$formula[[3]]
    if(!is.name(xvar_expr))
      stop("Please provide xvar. A default xvar could not be chosen.")
  } else {
    xvar_expr <- substitute(xvar)
  }
  xvar <- as.character(xvar_expr)
  if(missing(xlab)) xlab <- xvar
  pred <- predict_smm(object, include_smooths = TRUE,
                      include_REs = FALSE, include_residuals = FALSE)
  xdata <- object$model[[xvar]]
  order <- order(xdata)
  yvar_expr <- mmfit$formula[[2]]
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
                fill = "blue",alpha = 0.25) +
    ggplot2::geom_line(mapping = aesmm, data = df, inherit.aes=FALSE, color = color,linewidth=1)

  fig <- fig + ggplot2::xlab(xlab)
  fig <- fig + ggplot2::ylab(ylab)
  if(!missing(title)) fig <- fig + ggplot2::ggtitle(title)
  if(!missing(ylim)) fig <- fig + ggplot2::coord_cartesian(ylim = ylim)
  fig
}
