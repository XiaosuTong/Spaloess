

spaloess <- function (formula, data, weights, subset, na.action, model = FALSE, 
    span = 0.75, enp.target, degree = 2L, parametric = FALSE, 
    drop.square = FALSE, normalize = FALSE, family = c("gaussian", 
        "symmetric"), method = c("loess", "model.frame"), control = loess.control(...), 
    ...) 
{
    family <- match.arg(family)
    print(family)
    method <- match.arg(method)
    mf <- match.call(expand.dots = FALSE)
    mf$model <- mf$span <- mf$enp.target <- mf$degree <- mf$parametric <- mf$drop.square <- mf$normalize <- mf$family <- mf$method <- mf$control <- mf$... <- NULL
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    if (match.arg(method) == "model.frame") 
        return(mf)
    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    w <- model.weights(mf)
    if (is.null(w)) 
        w <- rep(1, length(y))
    nmx <- as.character(attr(mt, "variables"))[-(1L:2)]
    x <- mf[, nmx, drop = FALSE]
    if (any(sapply(x, is.factor))) 
        stop("predictors must all be numeric")
    x <- as.matrix(x)
    D <- ncol(x) 
    nmx <- colnames(x)
    names(nmx) <- nmx
    drop.square <- match(nmx, nmx[drop.square], 0L) > 0L 
    parametric <- match(nmx, nmx[parametric], 0L) > 0L 
    if (!match(degree, 0L:2L, 0L)) 
        stop("'degree' must be 0, 1 or 2")
    iterations <- if (family == "gaussian") 
        1
    else control$iterations
    if (!missing(enp.target)) 
        if (!missing(span)) 
            warning("both 'span' and 'enp.target' specified: 'span' will be used")
        else {
            tau <- switch(degree + 1L, 1, D + 1, (D + 1) * (D + 
                2)/2) - sum(drop.square)
            span <- 1.2 * tau/enp.target #calculate the span based on enp.target 
        }
    if (!is.list(control) || !is.character(control$surface) || 
        !is.character(control$statistics) || !is.character(control$trace.hat) || 
        !is.numeric(control$cell) || !is.numeric(iterations)) 
        stop("invalid 'control' argument")
    fit <- newsimpleLoess(y, x, w, span, degree, parametric, drop.square, 
        normalize, control$statistics, control$surface, control$cell, 
        iterations, control$trace.hat)
    fit$call <- match.call()
    fit$terms <- mt
    fit$xnames <- nmx
    fit$x <- x
    fit$y <- y
    fit$weights <- w
    if (model) 
       fit$model <- mf 
    fit$na.action <- attr(mf, "na.action")
    fit
}

newsimpleLoess <- function (y, x, weights, span = 0.75, degree = 2L, parametric = FALSE, 
    drop.square = FALSE, normalize = TRUE, statistics = "approximate", 
    surface = "interpolate", cell = 0.2, iterations = 1L, trace.hat = "exact") 
{
    D <- as.integer(NCOL(x))
    if (is.na(D)) 
        stop("invalid NCOL(X)")
    if (D > 4) 
        stop("only 1-4 predictors are allowed")
    N <- as.integer(NROW(x))
    if (is.na(N)) 
        stop("invalid NCOL(X)")
    if (!N || !D) 
        stop("invalid 'x'")
    if (!length(y)) 
        stop("invalid 'y'")
    x <- as.matrix(x)
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"
    storage.mode(weights) <- "double"
    max.kd <- max(N, 200)
    robust <- rep(1, N)
    divisor <- rep(1, D)
    if (normalize && D > 1L) {
        trim <- ceiling(0.1 * N)
        divisor <- sqrt(apply(apply(x, 2L, sort)[seq(trim + 1, 
            N - trim), , drop = FALSE], 2L, var))
        x <- x/rep(divisor, rep(N, D))
    }
    sum.drop.sqr <- sum(drop.square)
    sum.parametric <- sum(parametric)
    nonparametric <- sum(!parametric)
    order.parametric <- order(parametric)
    x <- x[, order.parametric]
    order.drop.sqr <- (2L - drop.square)[order.parametric]
    if (degree == 1L && sum.drop.sqr) 
        stop("specified the square of a factor predictor to be dropped when degree = 1")
    if (D == 1L && sum.drop.sqr) 
        stop("specified the square of a predictor to be dropped with only one numeric predictor")
    if (sum.parametric == D) 
        stop("specified parametric for all predictors")
    if((D - sum.parametric) != 2) 
        stop("dimension for distance must be 2 in spatial loess")
    if (iterations) 
        for (j in seq_len(iterations)) {
            robust <- weights * robust
            if (j > 1) 
                statistics <- "none"
            else if (surface == "interpolate" && statistics == 
                "approximate") 
                statistics <- if (trace.hat == "exact") 
                  "1.approx"
                else "2.approx"
            surf.stat <- paste(surface, statistics, sep = "/")
            if (length(span) != 1L) 
                stop("invalid argument 'span'")
            if (length(cell) != 1L) 
                stop("invalid argument 'cell'")
            if (length(degree) != 1L) 
                stop("invalid argument 'degree'")
            myv = rep(0, 5000)
            z <- .C("loess_raw", y, x, weights, robust, D, N, 
                as.double(span), as.integer(degree), as.integer(nonparametric), 
                as.integer(order.drop.sqr), as.integer(sum.drop.sqr), 
                as.double(span * cell), as.character(surf.stat), 
                fitted.values = double(N), parameter = integer(7L),
                a = integer(max.kd), xi = double(max.kd), vert = double(2L * 
                  D), vval = double((D + 1L) * max.kd), diagonal = double(N), 
                trL = double(1L), delta1 = double(1L), delta2 = double(1L), 
                as.integer(surf.stat == "interpolate/exact"), vert2 = double(D * max.kd))
            if (j == 1) {
                trace.hat.out <- z$trL
                one.delta <- z$delta1
                two.delta <- z$delta2
            }
#################################################################
#            fitted.residuals <- y - z$fitted.values
#            if (j < iterations) 
#                robust <- .Fortran(C_lowesw, fitted.residuals, 
#                  N, robust = double(N), integer(N))$robust
##################################################################
        }
    if (surface == "interpolate") {
        pars <- setNames(z$parameter, c("d", "n", "vc", "nc", 
            "nv", "liv", "lv"))
        enough <- (D + 1L) * pars["nv"]
        ## nv is the number of vertex
        fit.kd <- list(
            parameter = pars, 
            a = z$a[1L:pars[4L]], 
            xi = z$xi[1L:pars[4L]], 
            vert = z$vert, 
            vval = z$vval[1L:enough],
            vert2 = data.frame(matrix(z$vert2, ncol = D))[1:pars["nv"],] 
        )
    }
##########################################################################################
#    if (iterations > 1L) {
#        pseudovalues <- .Fortran("lowesp", N, as.double(y), as.double(z$fitted.values), 
#            as.double(weights), as.double(robust), integer(N), 
#            pseudovalues = double(N))$pseudovalues
#        zz <- .C("loess_raw", as.double(pseudovalues), x, weights, 
#            weights, D, N, as.double(span), as.integer(degree), 
#            as.integer(nonparametric), as.integer(order.drop.sqr), 
#            as.integer(sum.drop.sqr), as.double(span * cell), 
#            as.character(surf.stat), temp = double(N), parameter = integer(7L), 
#            a = integer(max.kd), xi = double(max.kd), vert = double(2L * 
#                D), vval = double((D + 1L) * max.kd), diagonal = double(N), 
#            trL = double(1L), delta1 = double(1L), delta2 = double(1L), 
#            0L)
#        pseudo.resid <- pseudovalues - zz$temp
#    }
#    sum.squares <- if (iterations <= 1L) 
#        sum(weights * fitted.residuals^2)
#    else sum(weights * pseudo.resid^2)
##########################################################################################
    enp <- one.delta + 2 * trace.hat.out - N
#########################################
#    s <- sqrt(sum.squares/one.delta)
#########################################
    pars <- list(
        robust = robust, 
        span = span, 
        degree = degree, 
        normalize = normalize, 
        parametric = parametric, 
        drop.square = drop.square, 
        surface = surface, 
        cell = cell, 
        family = if (iterations <= 1L) "gaussian" else "symmetric", 
        iterations = iterations
    )
    fit <- list(
        n = N, 
###########################################
#        fitted = z$fitted.values, 
#        residuals = fitted.residuals, 
###########################################
        enp = enp, 
#####################
#        s = s, 
#####################
        s = NULL,
        one.delta = one.delta, 
        two.delta = two.delta, 
        trace.hat = trace.hat.out, 
        divisor = divisor
    )
    fit$pars <- pars
    if (surface == "interpolate") 
        fit$kd <- fit.kd
    class(fit) <- "loess"
    fit
}
