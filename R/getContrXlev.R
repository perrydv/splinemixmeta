#' Copyright (C) 2025 Antonio Gasparrini, Francesco Sera
#' Copied from mixmeta version 1.2.2
#' Licensed under GPL (>= 3)
#' 
#' @keywords internal
#' @noRd
getContrXlev <-
function(terms, list)  {
#
################################################################################
# FUNCTION TO EXTRACT THE CONTRASTS/LEVELS RELATED TO A GIVEN FORMULA
#
  # IF NO CONSTRASTS, RETURN NULL
  if(is.null(list)) return(NULL)
#
  # RETURN CONSTRASTS RELATED TO TERMS IN THE FORMULA
  vars <- vapply(attr(terms, "variables"), deparse, "")[-1L]
  ind <- names(list) %in% vars
  if(any(ind)) list[ind] else NULL
}