#' Copyright (C) 2025 Antonio Gasparrini, Francesco Sera
#' Copied from mixmeta version 1.2.2
#' Licensed under GPL (>= 3)
#' 
#' @keywords internal
#' @noRd
xpndMat <-
function(vech) {
#
################################################################################
#
  dim <- (-1+sqrt(1+8*length(vech)))/2
  if(dim%%1!=0L) stop("dimension of 'vech' not consistent")
  mat <- matrix(nrow=dim,ncol=dim)
  mat[lower.tri(mat,diag=TRUE)] <- as.matrix(vech)
  mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
#
  mat
}