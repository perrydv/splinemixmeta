#' Fit a mixed-effects meta-analysis model with inclusion of a spline term
#'
#' @param smooth A smoothing term created by `mgcv::s()`
#' @param ... Additional arguments to be passed to `mixmeta::mixmeta()`
#'
#' @returns An object of class `mixmeta`
#' @export
#'
#' @importFrom mixmeta mixmeta
splinemixmeta <- function(smooth = NULL, formula, se, S=se^2, manual_fixed=FALSE, residual_re=TRUE, data, random = list(),
                          method="reml", bscov="unstr", ...) {
  input_call <- match.call(expand.dots = FALSE)
  formula_expr <- input_call$formula
  if(is.numeric(formula_expr)) {
    stop("formula must be a name or a formula. It can't be a numeric object itself because splinemixmeta must work with a formula.")
  }
  if(is.name(formula_expr))
    formula <- as.formula(paste0(formula_expr, " ~ 1"))
  ## Next steps
  ## To do:
  ##  - Arrange for a user to easily determine the order of random effects terms from splines and residuals
  ##  - Tell a user to change the order if they want to extract a different sequence of results.
  ##. - Explain manual_fixed
  ##. - Compare to blups from prototype version
  ##  - Write more tests including for bivariate splines
  ##  - fill in roxygen and make a vignette
  ##  - Test when the formula was made in a different environment
  ##. - mixmeta gets an error from the blup invocation within vcov
  eval_env <- new.env(parent = parent.frame())
  mixmeta_call <- input_call
  if(is.null(input_call$S)) {
    if(is.null(input_call$se)) stop("One of se or S must be provided.")
    mixmeta_call$S <- substitute((SE)^2, list(SE = input_call$se))
  }
  mixmeta_call$smooth <- NULL
  mixmeta_call$se <- NULL
  mixmeta_call$manual_fixed <- NULL
  mixmeta_call$residual_re <- NULL
  mixmeta_call[[1]] <- quote(mixmeta::mixmeta)
  if(is.null(smooth)) {
    res <- eval(mixmeta_call, envir = eval_env)
    return(res)
  }

  # The result of mgcv::s is a list, but we want a list of
  # such results. We'll assume that if an input list has a
  # class attribute, it was returned from mgcv::s.
  if(!is.null(attr(smooth, "class"))) {
    smooth <- list(smooth)
  }

  if(method != "reml")
    stop("Only 'reml' method is supported. The argument is included in this function in order to match mixmeta()'s arguments after 'smooth'.")

  if(!is.null(random)) {
    if(!is.list(random)) {
      random <- list(random)
    }
  } else {
    random <- list()
  }
  # A single bscov is recycled in mixmeta,
  # but we expand it manually so that we can add new
  # entries for the smooth terms
  if(length(bscov) == 1) {
    bscov <- rep(bscov, length(random))
  } else {
    if(length(bscov) != length(random))
      stop("Length of 'bscov' must be 1 or equal to the length of 'random'.")
  }
  smooth_var_names <- smooth |> lapply(\(s) s$term) |> unlist() |> unique()
  random_var_names <- random |> lapply(all.vars) |> unlist() |> unique()
  var_names <- smooth_var_names |> c(all.vars(formula), random_var_names) |> unique()

  # Create basis functions for each smooth term
  pf <- parent.frame()
  smooth_components <- smooth |>
    lapply(make_smm_smooth, data = data, vnames = var_names, envir=pf, manual_fixed = manual_fixed) |>
    unlist(recursive = FALSE)

  missing_data <- missing(data)
  names_data <- if(missing_data) "" else names(data)

  # browser()
  # Determine number of rows in data
  # Assume the formula is a Y name or has a Y name on the LHS

  if(is.name(formula_expr)) Yname <- formula_expr
  else {
    if(!is.name(formula[[2]]))
      stop("The LHS of the formula must be a single variable name.")
    Yname <- formula[[2]]
  }

  if(!missing(data)) {
    Yobj <- eval(Yname, envir = data, enclos = pf)
  } else {
    Yobj <- eval(Yname, envir = pf)
  }
  if(!is.numeric(Yobj))
    stop("Could not find the data for the response variable in the formula.")

  nrow_data <- NROW(Yobj)
  rm(Yobj)

  random_name <- function(len)
    sample(letters, len, replace=TRUE) |> paste0(collapse="")
  find_available_name <- function(proposed_names) {
    for(name in proposed_names) {
      if(name %in% var_names) next;
      if(name %in% names_data) next;
      return(name)
    }
    stop("Could not generate a new name for a grouping factor. Something is likely wrong.")
  }

  allname <- find_available_name(c("all", "All", "ALL", "all_", "All_", "ALL_",
                                           "smooth_grouping",
                                           paste0("all_", random_name(4))))
  eval_env[[allname]] <- as.factor(rep(allname, nrow_data))

  # Versions of fits needed for later potential blup and prediction requests:
  # mmfit model

  smooth_random <- list()

  for(i in seq_along(smooth_components)) {
    bFname <- paste0("basisFxns_", LETTERS[i])
    xFname <- paste0("s_fixed_", LETTERS[i])
    eval_env[[bFname]] <- smooth_components[[i]]$basisFxns
    new_random <- substitute(~ BFN - 1 | ALLNAME,
                             list(BFN = as.name(bFname),
                                  ALLNAME = as.name(allname))) |>
      eval() |> list()
    environment(new_random[[1]]) <- eval_env
    smooth_random <- c(smooth_random, new_random)
    bscov <- c(bscov, "id")
    if(!manual_fixed) {
      if(!is.null(smooth_components[[i]]$x_fixed)) {
        eval_env[[xFname]] <- smooth_components[[i]]$x_fixed
        # Add fixed effects design matrix to dotsList
        formula <- substitute(update(formula, ~ . + XFN),
                              list(XFN = as.name(xFname))) |>
          eval()
      }
    }
  }
  random <- c(smooth_random, random)

  residual_re <- isTRUE(residual_re)
  if(residual_re) {
    IDname <- find_available_name(c("ID", "id", "POINT_ID", "point_id",
                                    "obs_id", "OBS_ID",
                                    "row_id", "ROW_ID",
                                    paste0("ID_", random_name(4))))
    new_random <- substitute(~ 1 | IDNAME,
                             list(IDNAME = as.name(IDname))) |>
      eval() |> list()
    environment(new_random[[1]]) <- eval_env
    random <- c(random, new_random)
    eval_env[[IDname]] <- as.factor(1:nrow_data); bscov <- c(bscov, "unstr")
  }
  mixmeta_call$random <- random
  mixmeta_call$bscov <- bscov
  mixmeta_call$formula <- formula
  environment(mixmeta_call$formula) <- eval_env
  # browser()

  mmfit <- eval(mixmeta_call, envir = eval_env)
  mmfit
  # mixmeta::mixmeta sets data <- parent.frame() if data was missing
  # If we simply call mixmeta::mixmeta(...) from here, this function
  # will be the parent.frame, and that may not make sense.
  # Therefore, we construct an environment in which to evaluate
  # a call to mixmeta::mixmeta and set its parent frame so scoping
  # will find things in this call's parent.frame().
}
