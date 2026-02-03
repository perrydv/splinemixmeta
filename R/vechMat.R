#' Copyright (C) 2025 Antonio Gasparrini, Francesco Sera
#' Copied from mixmeta version 1.2.2
#' Licensed under GPL (>= 3)
#' 
#' @keywords internal
#' @noRd
vechMat <-
function(mat, diag=TRUE) {
#
################################################################################
#
  if(!is.matrix(mat)) mat <- as.matrix(mat)
  if(diff(dim(mat))!=0) stop("Non-square matrix")
#
  mat[lower.tri(mat,diag=diag)]
}