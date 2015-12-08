#' Spatially Local Polynomial Regression Prediction
#'
#' The second layer of the prediction of Spatial locally weighted regression. It calls 
#' the computation engine functions wrote in C and Fortran, Concretely, the "loess_dfit"
#' and "loess_ifit". Notice, it has to be called by predloess function, do not call 
#' newPredLoess directly.
#'
#' @param y
#'    the response variable in the form of vector
#' @param x
#'    a data.frame of the original design matrix 
#' @param newx
#'    a data.frame containing the new locations for prediction
#' @param s
#'    standard error from the orginial local fit. If se = FALSE, s will not be useful.
#' @param weights
#'    optional initial weights for the original data set. Get from spaloess$weights.
#' @param robust
#'    robust factor to each location of original data set, get from spaloess$pars$robust.
#' @param span
#'    The parameter alpha which controls the portion of data points used in the local fit.
#' @param degree
#'    The degree of the polynomials to be used, normally 1 or 2. (Degree 0 is also allowed, but see 
#'    the ‘Note’.)
#' @param distance 
#'    Specify the method of distance calculation. Options: "Euclid", or "Latlong" which is for great 
#'    circle distance
#' @param parametric 
#'    should any terms be fitted globally rather than locally? Terms can be specified by name, 
#'    number or as a logical vector of the same length as the number of predictors.
#' @param drop.square
#'    For fits with more than one predictor and 'degree = 2', should the quadratic term be dropped 
#'    for particular predictors?  Terms are specified in the same way as for 'parametric'.
#' @param normalize
#'    Should the predictors be normalized to a common scale if there is more than one?  The 
#'    normalization used is to set the 10% trimmed standard deviation to one. Set to false for 
#'    "Latlong" distance.
#' @param surface
#'    should the fitted surface be computed exactly or via interpolation from a kd tree? 
#'    can be either "interpolate" or "direct", default is "interpolate".
#' @param cell
#'    if interpolation is used this controls the accuracy of the approximation via the maximum 
#'    number of points in a cell in the kd tree. Cells with more than ‘floor(n*span*cell)’ points
#'    are subdivided.
#' @param family 
#'    If 'gaussian' fitting is by least-squares, and if 'symmetric' a re-descending M estimator is 
#'    used with Tukey's bi-weight function.
#' @param kd
#'    a list object which contains all information about the kd-tree built based on original data 
#'    set.
#' @param divisor
#'    the scaler for the normalization. x/divisor if normalize = TRUE.
#' @param se 
#'    should standard errors be computed? Default is FALSE
#' @param ...
#'    arguments passed to or from other methods.
#' @author 
#'    Xiaosu Tong, based on 'loess' function of B. D. Ripley, and 'cloess' package of Cleveland,
#'    Grosse and Shyu. 
#' @export

newPredLoess <- function (y, x, newx, s, weights, robust, span, degree, normalize,
    parametric, drop.square, surface, cell, family, kd, divisor, se = FALSE, distance)
{
    D <- NCOL(x)
    N <- NROW(x)
    M <- NROW(newx)
    x <- as.matrix(x)
    newx <- as.matrix(newx)
    newx <- newx/rep(divisor, rep(M, D))
    x <- x/rep(divisor, rep(N, D))
    sum.drop.sqr <- sum(drop.square)
    nonparametric <- sum(!parametric)
    order.parametric <- order(parametric)
    x <- x[, order.parametric, drop = FALSE]
    x.evaluate <- newx[, order.parametric, drop = FALSE]
    order.drop.sqr <- (2L - drop.square)[order.parametric]
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"

    if (distance == "Latlong") {
        xtdist <- 1
    } else if (distance == "Euclid") {
        xtdist <- 0
    }

    if (surface == "direct") {
        nas <- rowSums(is.na(newx)) > 0L
        fit <- rep(NA_real_, length(nas))
        x.evaluate <- x.evaluate[!nas, , drop = FALSE]
        M <- nrow(x.evaluate)
        if (se) {
            se.fit <- fit
            z <- .C("loess_dfitse", y, x, as.double(x.evaluate),
                as.double(weights * robust), as.double(robust),
                as.integer(family == "gaussian"), as.double(span),
                as.integer(degree), as.integer(nonparametric),
                as.integer(order.drop.sqr), as.integer(sum.drop.sqr),
                as.integer(D), as.integer(N), as.integer(M),
                fit = double(M), L = double(N * M), as.integer(xtdist))[c("fit",
                "L")]
            fit[!nas] <- z$fit
            ses <- (matrix(z$L^2, M, N)/rep(weights, rep(M, N))) %*%
                rep(1, N)
            se.fit[!nas] <- drop(s * sqrt(ses))
        }
        else {
            fit[!nas] <- .C("loess_dfit", y, x, as.double(x.evaluate),
                as.double(weights * robust), as.double(span),
                as.integer(degree), as.integer(nonparametric),
                as.integer(order.drop.sqr), as.integer(sum.drop.sqr),
                as.integer(D), as.integer(N), as.integer(M),
                fit = double(M), as.integer(xtdist))$fit
        }
    }
    else {
        inside <- matrix(FALSE, M, ncol = D)
        ranges <- apply(x, 2L, range)
        inside <- (x.evaluate <= rep(ranges[2L, ], rep(M, D))) &
            (x.evaluate >= rep(ranges[1L, ], rep(M, D)))
        inside <- inside %*% rep(1, D) == D
        inside[is.na(inside)] <- FALSE
        #M1 is the number of new points want to be fitted
        M1 <- sum(inside)
        fit <- rep(NA_real_, M)
        if (any(inside))
            fit[inside] <- .C("loess_ifit", as.integer(kd$parameter),
                as.integer(kd$a), as.double(kd$xi), as.double(kd$vert),
                as.double(kd$vval), as.integer(M1), as.double(x.evaluate[inside,
                  ]), fit = double(M1), as.integer(xtdist))$fit
        if (se) {
            se.fit <- rep(NA_real_, M)
            if (any(inside)) {
                L <- .C("loess_ise", y, x, as.double(x.evaluate[inside,
                  ]), as.double(weights), as.double(span), as.integer(degree),
                  as.integer(nonparametric), as.integer(order.drop.sqr),
                  as.integer(sum.drop.sqr), as.double(span *
                    cell), as.integer(D), as.integer(N), as.integer(M1),
                  double(M1), L = double(N * M1), as.integer(xtdist))$L
                tmp <- (matrix(L^2, M1, N)/rep(weights, rep(M1,
                  N))) %*% rep(1, N)
                se.fit[inside] <- drop(s * sqrt(tmp))
            }
        }
    }
    rn <- rownames(newx)
    if (se) {
        if (!is.null(rn))
            names(fit) <- names(se.fit) <- rn
        list(fit = fit, se.fit = drop(se.fit), residual.scale = s)
    }
    else {
        if (!is.null(rn))
            names(fit) <- rn
        fit
    }
}
