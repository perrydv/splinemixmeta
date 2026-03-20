#' Copyright (C) 2025 Antonio Gasparrini, Francesco Sera
#' Copied from mixmeta version 1.2.2
#' Licensed under GPL (>= 3)
#' 
#' @keywords internal
#' @noRd
sumList <-
function(list) {
#
################################################################################
#
  res <- 0
  for(i in seq(list)) res <- res + list[[i]]
#
  res
}