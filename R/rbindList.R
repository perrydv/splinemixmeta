#' Copyright (C) 2025 Antonio Gasparrini, Francesco Sera
#' Copied from mixmeta version 1.2.2
#' Licensed under GPL (>= 3)
#' 
#' @keywords internal
#' @noRd
rbindList <-
function(list, ncol) {
#
################################################################################
# FUNCTION TO EXTRACT THE PARAMETERS DEFINING THE RANDOM PART
#
  # DEFINE IF MATRIX AND COLUMNS
  ismat <- !is.null(dim(list[[1]]))
#  
  # IF A MATRIX, TRANSPOSE AND GET THE ELEMENTS
  if(ismat) list <- lapply(list, function(x) c(t(x)))
#
  # DEFINE THE MATRIX
  matrix(unlist(list), ncol=ncol, byrow=T)
}