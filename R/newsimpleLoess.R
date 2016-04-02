#' Spatially Local Polynomial Regression Fitting
#'
#' The second layer of the Spatial locally weighted regression. It calls the computation engine 
#' functions wrote in C and Fortran. Notice, it has to be called by spaloess function, do not
#' call newsimpleLoess directly.
#'
#' @param y
#'    the response variable in the form of vector
#' @param x
#'    a data.frame of the design matrix 
#' @param weights
#'    optional initial weights for each case
#' @param span
#'    The parameter alpha which controls the portion of data points used in the local fit.
#' @param degree
#'    The degree of the polynomials to be used, normally 1 or 2. (Degree 0 is also allowed, but see 
#'    the 'Note'.)
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
#' @param statistics
#'    should the statistics be computed exactly or approximately? 
#'    Exact computation can be very slow. Can be either "approximate" or "exact", default is
#'    "approximate".
#' @param surface
#'    should the fitted surface be computed exactly or via interpolation from a kd tree? 
#'    can be either "interpolate" or "direct", default is "interpolate".
#' @param cell
#'    if interpolation is used this controls the accuracy of the approximation via the maximum 
#'    number of points in a cell in the kd tree. Cells with more than 'floor(n*span*cell)' points
#'    are subdivided.
#' @param iterations
#'    the number of iterations used in robust fitting.
#' @param trace.hat
#'    should the trace of the smoother matrix be computed exactly or approximately? 
#'    It is recommended to use the approximation for more than about 1000 data points.
#' @param ...
#'    arguments passed to or from other methods.
#' @author 
#'    Xiaosu Tong, based on 'loess' function of B. D. Ripley, and 'cloess' package of Cleveland,
#'    Grosse and Shyu. 
#' @export
#' @details
#'    The second layer of the Spatial locally weigted regression fitting. In this function, all
#'    input arguments which are critical for the fit are checked, like the dimension of design
#'    matrix can only be less than or equal to 4. Or input x can only by lat and lon if distance
#'    is "Latlong".
#'    Fitting is done locally.  That is, for the fit at point x, the fit is made using points in 
#'    a neighborhood of x, weighted by their distance from x (with differences in 'parametric' 
#'    variables being ignored when computing the distance).  The size of the neighborhood is 
#'    controlled by alpha (set by 'span' or 'enp.target').  For alpha < 1, the neighborhood includes
#'    proportion alpha of the points, and these have tricubic weighting (proportional to 
#'    (1 - (dist/maxdist)^3)^3).  For alpha > 1, all points are used, with the 'maximum distance' 
#'    assumed to be alpha^(1/p) times the actual maximum distance for p explanatory variables. For 
#'    the default family, fitting is by (weighted) least squares. For 'family="symmetric"' a few 
#'    iterations of an M-estimation procedure with Tukey's bi-weight are used.  Be aware that as the
#'    initial value is the least-squares fit, this need not be a very resistant fit. It can be 
#'    important to tune the control list to achieve acceptable speed.  See 'loess.control' for details.

newsimpleLoess <- function (y, x, allx, weights, span = 0.05, degree = 2L, distance = "Latlong", parametric = FALSE, 
    drop.square = FALSE, normalize = FALSE, statistics = "approximate", 
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
    alN <- as.integer(NROW(allx))
    if (N > alN)
        stop("invalid NROW(alN)")
    if (!N || !D) 
        stop("invalid 'x'")
    if (!length(y)) 
        stop("invalid 'y'")
    x <- as.matrix(x)
    allx <- as.matrix(allx)
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"
    storage.mode(weights) <- "double"
    max.kd <- max(N, 200)
    robust <- rep(1, N)
    divisor <- rep(1, D)

    ## if normalize is TRUE, the scaler factor "divisor" is calculated based on the sd of truncated data betwen 5% - 95% quantiles
    ## set normalization with distance="Euclid" is same as we choose distance="Mahal" with M = SIGMA^(-1), which SIGMA is a diagonal
    ## matrix with variance of each predictor as diagnonal elements.
    if (normalize && D > 1L && distance != "Latlong") {
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
    if (distance == "Latlong" & (D - sum.parametric) != 2) 
        stop("dimension for great circle distance must be 2 in spatial loess")
    
    if (distance == "Latlong") {
        xtdist <- 1
    } else if (distance == "Euclid") {
        xtdist <- 0
    }

    ## The real computation engine of loess is in c code named "loess_raw"
    if (iterations) {
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

            ## the 7th argument is added by Xiaosu Tong, which is the distance calculation flag
            ## the 8th argument is added by Xiaosu Tong, which is passing all x locations to the kd-tree
            ## the 9th argument is added by Xiaosu Tong, which is the nrow of all x location
            z <- .C("loess_raw", y, x, weights, robust, D, N, as.integer(xtdist), allx, alN,
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

            fitted.residuals <- y - z$fitted.values

            ## if current iteration is less than the total number of iteration time, calculate the
            ## robust factor based on fitted residuals.
            if (j < iterations) { 
                robust <- .Fortran("lowesw", fitted.residuals, N, robust = double(N), integer(N))$robust
            }
        }
    }

    ## if the calculation is interpolate which is fitting at vertices of kd-tree
    ## and reture a list object named fit.kd which including all kd-tree information
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
    if (iterations > 1L) {                                                                
        pseudovalues <- .Fortran("lowesp", N, as.double(y), as.double(z$fitted.values),   
            as.double(weights), as.double(robust), integer(N),                            
            pseudovalues = double(N))$pseudovalues    
        ## the 7th argument is added by Xiaosu Tong, which is the distance calculation flag
        zz <- .C("loess_raw", as.double(pseudovalues), x, weights,                        
            weights, D, N, as.integer(xtdist), as.double(span), as.integer(degree),         
            as.integer(nonparametric), as.integer(order.drop.sqr),                        
            as.integer(sum.drop.sqr), as.double(span * cell),                             
            as.character(surf.stat), temp = double(N), parameter = integer(7L),           
            a = integer(max.kd), xi = double(max.kd), vert = double(2L *                  
                D), vval = double((D + 1L) * max.kd), diagonal = double(N),               
            trL = double(1L), delta1 = double(1L), delta2 = double(1L),                   
            0L)                                                                           
        pseudo.resid <- pseudovalues - zz$temp                                            
    }                                                                                     
    sum.squares <- if (iterations <= 1L)                                                  
        sum(weights * fitted.residuals^2)                                                 
    else sum(weights * pseudo.resid^2)                                                    
    enp <- one.delta + 2 * trace.hat.out - N                                              
                                                                                          
    s <- sqrt(sum.squares/one.delta)                                                      
    pars <- list(
        robust = robust, 
        span = span, 
        degree = degree, 
        normalize = normalize, 
        parametric = parametric, 
        drop.square = drop.square, 
        surface = surface, 
        cell = cell, 
        distance = distance, 
        family = if (iterations <= 1L) "gaussian" else "symmetric", 
        iterations = iterations
    )
    fit <- list(
        n = N, 
        fitted = z$fitted.values, 
        residuals = fitted.residuals,
        enp = enp,  
        s = s,      
        one.delta = one.delta, 
        two.delta = two.delta, 
        trace.hat = trace.hat.out, 
        divisor = divisor
    )
    fit$pars <- pars
    if (surface == "interpolate") 
        fit$kd <- fit.kd
    class(fit) <- "spaloess"
    fit
}
