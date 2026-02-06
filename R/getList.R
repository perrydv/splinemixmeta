#' Copyright (C) 2025 Antonio Gasparrini, Francesco Sera
#' Copied from mixmeta version 1.2.2
#' Licensed under GPL (>= 3)
#' 
#' @keywords internal
#' @noRd
getList <-
function(object) {
#
################################################################################
# TRANFORM THE OBJECT IN A LIST
#
  if(is.list(object)) object else list(object)
}