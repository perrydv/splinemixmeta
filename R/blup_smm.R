#' Obtain predictions (BLUPs) from a `splinemixmeta` object
#'
#' @param object An object of class `splinemixmeta`, returned by `splinemixmeta()`.
#' @param se Logical indicating whether to return standard errors of the predictions.
#' @param pi Logical indicating whether to return prediction intervals.
#' @param vcov Logical indicating whether to return the variance-covariance matrix of the predictions.
#' @param pi.level Numeric value between 0 and 1 indicating the confidence level for the prediction intervals. Default is 0.95.
#' @param type Character string specifying the type of prediction: "outcome" for predicted outcomes or "residual" for predicted residuals.
#' @param level Integer indicating the random effects level for which to obtain predictions. Default is the highest level.
#' @param format Character string specifying the format of the output: "matrix" or "list". Default depends on the number of outcomes and whether `vcov` is requested.
#' @param aggregate Character string specifying how to aggregate the output when multiple outcomes are present
#' @param ... Additional arguments (currently unused).
#' @return A matrix or list of predicted values (BLUPs), optionally including standard errors, prediction intervals, and variance-covariance matrices.

#' @details This function is modified from `mixmeta::blup.mixmeta` with acknowledgement of the original authors.
#' It is modified to handle intermediate levels of random effects more carefully.
#' This is currently a temporary solution that may be replaced in future versions.
#'
#' @export
blup.splinemixmeta <- function (object, se = FALSE, pi = FALSE, vcov = FALSE, pi.level = 0.95,
    type = "outcome", level, format, aggregate = "stat", ...)
{
    if (missing(format))
        format <- ifelse(vcov && object$dim$k > 1, "list", "matrix")
    type <- match.arg(type, c("outcome", "residual"))
    format <- match.arg(format, c("matrix", "list"))
    aggregate <- match.arg(aggregate, c("stat", "outcome"))
    maxlevel <- length(object$dim$q) # added by Perry
    if (missing(level))
        level <- length(object$dim$q)
    if (object$method == "fixed")
        level <- 0L
    if (level < 0L || level > length(object$dim$q))
        stop("'level' non compatible with random levels")
    mf <- model.frame(object)
    na.action <- object$na.action
    nm <- rownames(mf)
    groups <- getGroups(object$random, mf)
    ord <- do.call(order, lapply(seq(ncol(groups)), function(i) groups[,
        i]))
    groups <- groups[ord, seq(max(1, level)), drop = FALSE]
    allGroups <- getGroups(object$random, mf)
    allGroups <- allGroups[ord, seq(maxlevel), drop = FALSE] # added by Perry
    mf <- mf[ord, , drop = FALSE]
    y <- as.matrix(model.response(mf, "numeric"))
    offset <- model.offset(mf)
    if (!is.null(offset))
        y <- y - offset
    X <- model.matrix(object)[ord, , drop = FALSE]
    Z <- if (level > 0L)
      getZ(object$random, mf, object$contrasts)
    else NULL
    allZ <- Z # added by Perry
    if (!is.null(Z) && is.list(Z))
        Z <- if (level == 1L)
            Z[[1L]]
        else Z[seq(level)]
    S <- if (!is.null(object$S))
        as.matrix(object$S)[ord, , drop = FALSE]
    else NULL
    m <- object$dim$m
    k <- object$dim$k
    q <- if (level == 0L)
        0
    else object$dim$q[seq(level)]
    allq <- object$dim$q[seq(maxlevel)] # added by Perry
    gp <- groups[, 1]
    rep <- do.call(cbind, lapply(seq(ncol(groups)), function(i) tapply(groups[,
        i], gp, function(xi) length(unique(xi)))))
    nay <- is.na(y)
    if (any(nay)) {
        yS <- mixmeta::inputna(y, S)
        y <- yS[, seq(k), drop = FALSE]
        S <- yS[, -seq(k), drop = FALSE]
        nay[nay] <- FALSE
    }
    nalist <- lapply(seq(m), function(i) rep(FALSE, k))
    ylist <- lapply(seq(m), function(i) c(t(y[gp %in% i, ])))
    Xlist <- lapply(seq(m), function(i) X[gp %in% i, , drop = FALSE] %x%
        diag(k))
    Zlist <- getZlist(Z, nay, groups, m, k, q)

    allZlist <- getZlist(allZ, nay, allGroups, m, k, allq) # added by Perry

    Slist <- getSlist(S, nay, groups, m, k, object$control$addSlist,
        object$control$checkPD)

    allSlist <- getSlist(S, nay, allGroups, m, k, object$control$addSlist,
        object$control$checkPD) # added by Perry.

    predlist <- lapply(seq(m), function(i) Xlist[[i]] %*% object$coefficients)
    reslist <- mapply(function(y, pred) y - pred, ylist, predlist,
        SIMPLIFY = FALSE)
    if (level > 0L && !is.null(Psi <- object$Psi) && is.list(Psi))
        Psi <- Psi[seq(level)]
    ZPZlist <- if (level == 0L)
        NULL
    else getZPZlist(Zlist, nalist, Psi)

    allPsi <- object$Psi
    allZPZlist <- getZPZlist(allZlist, nalist, allPsi) # added by Perry
    # if (!is.null(ZPZlist)) {
    #     Ulist <- mapply(function(ZPZ, S) chol(ZPZ + S), ZPZlist,
    #         Slist, SIMPLIFY = FALSE)
    #     invUlist <- lapply(Ulist, function(U) backsolve(U, diag(ncol(U))))
    # }
    if (!is.null(allZPZlist)) {
        Ulist <- mapply(function(ZPZ, S) chol(ZPZ + S), allZPZlist,
            allSlist, SIMPLIFY = FALSE)
        invUlist <- lapply(Ulist, function(U) backsolve(U, diag(ncol(U))))
    }
    complist <- lapply(seq(m), function(i) {
        blup <- if (type == "residual")
            array(0, dim(predlist[[i]]))
        else predlist[[i]]
        V <- if (type == "residual")
            matrix(0, nrow(Xlist[[i]]), nrow(Xlist[[i]]))
        else Xlist[[i]] %*% tcrossprod(object$vcov, Xlist[[i]])
        if (!is.null(ZPZlist)) {
            ZPZinvSigma <- ZPZlist[[i]] %*% tcrossprod(invUlist[[i]])
            blup <- blup + ZPZinvSigma %*% reslist[[i]]
            V <- V + ZPZlist[[i]] - ZPZinvSigma %*% ZPZlist[[i]]
        }
        stderr <- sqrt(diag(V))
        seqlist <- lapply(seq(length(blup)/k), function(i) c(i *
            k - k + 1, i * k))
        blup <- rbindList(lapply(seqlist, function(x) blup[x[1]:x[2],
            ]), k)
        V <- rbindList(lapply(seqlist, function(x) mixmeta::vechMat(V[x[1]:x[2],
            x[1]:x[2]])), k * (k + 1)/2)
        stderr <- rbindList(lapply(seqlist, function(x) stderr[x[1]:x[2]]),
            k)
        return(list(blup, V, stderr))
    })
    blup <- rbindList(lapply(complist, "[[", 1), k)
    if (!is.null(offset) && type == "outcome")
        blup <- blup + offset
    blup <- blup[order(ord), , drop = FALSE]
    V <- rbindList(lapply(complist, "[[", 2), k * (k + 1)/2)[order(ord),
        , drop = FALSE]
    stderr <- rbindList(lapply(complist, "[[", 3), k)[order(ord),
        , drop = FALSE]
    colnames(blup) <- colnames(stderr) <- object$lab$k
    if (!is.null(na.action)) {
        blup <- napredict(na.action, blup)
        V <- napredict(na.action, V)
        stderr <- napredict(na.action, stderr)
        if (class(na.action) %in% c("exclude", "pass")) {
            nm <- napredict(na.action, nm)
            nm[na.action] <- names(na.action)
        }
    }
    zvalpi <- qnorm((1 - pi.level)/2, lower.tail = FALSE)
    pilb <- blup - zvalpi * stderr
    piub <- blup + zvalpi * stderr
    if (format == "matrix" && k == 1) {
        rownames(blup) <- nm
        if (se)
            blup <- cbind(blup, stderr)
        if (pi)
            blup <- cbind(blup, pilb, piub)
        if (vcov)
            blup <- cbind(blup, V)
        colnames(blup) <- c("blup", "se", "pi.lb", "pi.ub", "vcov")[c(TRUE,
            se, pi, pi, vcov)]
        if (ncol(blup) == 1L)
            blup <- drop(blup)
    }
    else if (format == "matrix" && any(se, pi, vcov) && aggregate ==
        "stat") {
        blup <- list(blup = blup)
        if (se)
            blup$se <- stderr
        if (pi)
            blup[c("pi.lb", "pi.ub")] <- list(pilb, piub)
        if (vcov)
            blup$vcov <- V
        for (i in seq(blup)) rownames(blup[[i]]) <- nm
    }
    else if (format == "matrix" && any(se, pi) && !vcov && aggregate ==
        "outcome") {
        rownames(blup) <- nm
        blup <- lapply(colnames(blup), function(j) cbind(blup = blup[,
            j], se = stderr[, j], pi.lb = pilb[, j], pi.ub = piub[,
            j])[, c(TRUE, se, pi, pi)])
        names(blup) <- object$lab$k
    }
    else if (format == "list" || (vcov && aggregate == "outcome")) {
        blup <- lapply(seq(nrow(blup)), function(i) {
            temp <- list(blup = blup[i, ])
            if (se)
                temp$se <- stderr[i, ]
            if (pi)
                temp[c("pi.lb", "pi.ub")] <- list(pilb[i, ],
                  piub[i, ])
            if (vcov) {
                temp$vcov <- mixmeta::xpndMat(V[i, ])
                dimnames(temp$vcov) <- list(object$lab$k, object$lab$k)
            }
            temp <- lapply(temp, function(x) if (any(is.na(x)))
                NA
            else x)
            return(dropList(temp))
        })
        names(blup) <- nm
    }
    else rownames(blup) <- nm
    blup
}
