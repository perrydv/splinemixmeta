#' Copyright (C) 2025 Antonio Gasparrini, Francesco Sera
#' Copied from mixmeta version 1.2.2
#' Licensed under GPL (>= 3)
#' 
#' @keywords internal
#' @noRd
dropList <-
function(object) {
#
################################################################################
# DROP THE LIST STRUCTURE IF THE LIST HAS ONLY ONE COMPONENT
#
  if(is.list(object) && length(object)==1L) object[[1L]] else object
}