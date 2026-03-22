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

#' @keywords internal
#' @noRd
getGroups <-
function(random, data)  {
#
################################################################################
# FUNCTION TO DEFINE THE GROUPING STRUCTURE FOR THE RANDOM PART
#
  # IF random IS NULL, RETURN JUST SEQUENCE OF ROWS (ONE OBS PER STUDY)
  if(is.null(random)) return(matrix(seq(nrow(data))))
#
  # EXTRACT A LIST WITH GROUPING VARIABLES
  random <- getList(random)
  groups <- lapply(random, function(form) {
    form[[2]] <- form[[2]][[3]]
    model.frame(form,data)[[1]]
  })
#
  # DEFINE GROUPING THROUGH FACTORS, ACCOUNTING FOR INTERNAL NESTING
  # NB: LEVELS FROM 2 ON ALWAYS NESTED IN 1, AS THE LATTER DEFINES THE LISTS
  groups[[1]] <- factor(groups[[1]])
  if((len <- length(groups))>1L) for(i in 2:len)
    groups[[i]] <- factor(paste(groups[[i-1]],groups[[i]],sep="-"))
#
  # TRANFORM IN MATRIX
  groups <- do.call(cbind,lapply(groups,unclass))
#
  groups
}

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

#' @keywords internal
#' @noRd
getSlist <-
function(S, nay, groups, m, k, addSlist=NULL, checkPD=NULL) {
#
################################################################################
# TRANSFORM S IN A LIST OF MATRICES
#
  # IF S IS PROVIDED
  if(!is.null(S)) {
#
    # CHECK THAT S AND addSlist ARE NOT BOTH PROVIDED
    if(!is.null(addSlist)) stop("'addSlist' only allowed without 'S'")
#
    # CREATE THE LIST
    Slist <- lapply(seq(m),function(i)
      mixmeta::bdiagMat(lapply(which(groups[,1]%in%i),function(j)
        mixmeta::xpndMat(S[j,])[!nay[j,],!nay[j,],drop=FALSE])))
    if(any(is.na(unlist(Slist))))
      stop("missing pattern in 'y' and S' is not consistent")
#
    # RETURN
    return(Slist)
  }
#
  # IF NOT INSTEAD, MUST BE PROVIDED
  # CHECKS (ASSUMED CONSISTENT WITH ORDER OF GROUPING AND DIMENSIONS)
  # NOT NULL, A LIST, NO MISSING, RIGHT LENGTH, RIGHT DIMENSIONS, POS-DEF
  if(is.null(addSlist))
    stop("within-unit errors must be provided either through 'S' or 'control'")
  if(!is.list(addSlist)||length(addSlist)!=m)
    stop(paste("'addSlist' not consistent with required format and",
      "grouping length. See help(mixmeta.control)"))
  if(any(is.na(unlist(addSlist)))) stop("no missing allowed in 'addSlist'")
  ind <- sapply(seq(m), function(i)
    dim(as.matrix(addSlist[[i]]))!=sum(!nay[groups[,1]%in%i,]))
  if(any(ind)) stop("wrong dimensions in 'addSlist'. See help(mixmetaControl)")
  # CHECK POSITIVE-DEFINITENESS (BY DEFAULT)
  if(is.null(checkPD) || checkPD)
    addSlist <- checkPD(addSlist,error=TRUE,label="Slist")
#
  # RETURN
  return(addSlist)
}

#' @keywords internal
#' @noRd
getZ <-
function(random, data, contrasts=NULL)  {
#
################################################################################
# FUNCTION TO DEFINE THE DESIGN MATRICES FOR THE RANDOM PART
#
  # IF random IS NULL, JUST RETURN NULL
  if(is.null(random)) return(NULL)
#
  # OTHERWISE, GENERATE THE LIST
  random <- getList(random)
  Z <- lapply(random, function(form) {
    # REMOVE THE GROUPING FACTOR FROM THE FORMULA
    form[[2]] <- form[[2]][[2]]
    # EXTRACT THE CONTRASTS
    contr <- getContrXlev(form, contrasts)
    # DERIVE THE MODEL MATRIX
    model.matrix(form, data, contr)
  })
#
  # RETURN A SINGLE MATRIX IF SINGLE LEVEL
  dropList(Z)
}

#' @keywords internal
#' @noRd
getZlist <-
function(Z, nay, groups, m, k, q)  {
#
################################################################################
# FUNCTION TO DEFINE THE LIST OF DESIGN MATRICES FOR THE RANDOM PART
#
  # IF NULL, RETURN SO, OTHERWISE TRANSFORM IN LIST
  if(is.null(Z)) return(NULL)
  Z <- getList(Z)
#
  # OTHERWISE, GENERATE THE LIST ACCOUNTING FOR MULTIPLE LEVELS
  # FOR EACH GROUP, A q-LENGTH LIST OF LIST OF MATRICES
  lapply(seq(m),function(i) {
    Zi <- lapply(Z, function(x) x[groups[,1]%in%i,,drop=FALSE])
    gi <- groups[groups[,1]%in%i,,drop=FALSE]
    nayi <- nay[groups[,1]%in%i,,drop=FALSE]
    Zij <- lapply(seq_along(q),function(j)
      lapply(unique(gi[,j]), function(rep) {
        ind <- gi[,j]%in%rep
        Zind <- Zi[[j]][ind,,drop=FALSE]%x%diag(k)
        Zind[!c(t(nayi[ind,])),,drop=FALSE]
      }))
  })
}

#' @keywords internal
#' @noRd
getZPZlist <-
function(Zlist, nalist, Psi)  {
#
################################################################################
# FUNCTION TO DEFINE THE LIST OF RANDOM-EFFECTS (CO)VARIANCE MATRICES
#
  # IF NO Psi, SIMPLY RETURN NULL
  if(is.null(Psi)) return(NULL)
#
  # IF NO Zlist, USE INFO IN nalist TO EXCLUDE MISSING
  if(is.null(Zlist)) return(lapply(nalist, function(na) Psi[!na,!na,drop=FALSE]))
#
  # OTHERWISE, IN EACH LEVEL MULTIPLY BY Z, BLOCK-DIAG, ADD BY LEVEL
  Psi <- getList(Psi)
  lapply(seq_along(Zlist), function(i)
    sumList(lapply(seq_along(Psi),function(j)
      mixmeta::bdiagMat(lapply(Zlist[[i]][[j]],function(x)
        x%*%Psi[[j]]%*%t(x))))))
}

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
