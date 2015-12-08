#' Spatially Local Polynomial Regression Prediction
#'
#' The first layer of the prediction of Spatial locally weighted regression. Mainly used as
#' prediction function for the NAs in the original data set.
#'
#' @param object
#'     an object fitted by ‘loess’.
#' @param newdata
#'     an optional data frame in which to look for variables with which to predict, or a 
#'     matrix or vector containing exactly the variables needs for prediction.  If missing, 
#'     the original data points are used.
#' @param se
#'     should standard errors be computed? Default is FALSE
#' @param na.action
#'     function determining what should be done with missing values in data frame ‘newdata’.  
#'     The default is to predict ‘NA’.
#' @param ...
#'     arguments passed to or from other methods.
#' @details
#'     This is the first layer of prediction function of spatial locally weigted regression.
#'     In the spaloess function, NA will be removed from the fitting. By passing the spaloess
#'     object and NA observations to predloess, predction at the locations of NA is carried out.
#'
#'     When the fit was made using ‘surface = "interpolate"’ (the default), ‘predictLoessf’ 
#'     will not extrapolate - so points outside an axis-aligned hypercube enclosing the 
#'     original data will have missing (‘NA’) predictions and standard errors.
#' @author 
#'     Xiaosu Tong, based on 'loess' function of B. D. Ripley, and 'cloess' package of Cleveland,
#'     Grosse and Shyu.  
#' @export
#' @examples
#'     set.seed(66)
#'     x1 <- rnorm(100, mean=-100, sd=10)
#'     x2 <- rnorm(100, mean=38, sd=4)
#'     y <- 0.1*x1 + 1*x2 - 10 + rnorm(100, 0, 1.3); y[1:2] <- NA
#'     testdata <- data.frame(LON = x1, LAT = x2, tmax = y)
#'     cars.lo <- spaloess(tmax ~ LON + LAT, testdata, distance = "Latlong")

predloess <- function (object, newdata = NULL, se = FALSE, na.action = na.pass, ...)
{
    if (!inherits(object, "spaloess")) {
        stop("first argument must be a \"spaloess\" object")
    }

    if (is.null(newdata) && !se) {
        return(fitted(object))
    }

    newx <- if (is.null(newdata)) {
        object$x
    }else if (is.data.frame(newdata)) {
        as.matrix(model.frame(delete.response(terms(object)), newdata, na.action = na.action))
    }else { 
        as.matrix(newdata)
    }

    res <- with(object, newPredLoess(y, x, newx, s, weights, pars$robust,
        pars$span, pars$degree, pars$normalize, pars$parametric,
        pars$drop.square, pars$surface, pars$cell, pars$family,
        kd, divisor, se = se, pars$distance))

    if (!is.null(out.attrs <- attr(newdata, "out.attrs"))) {
        if (se) {
            res$fit <- array(res$fit, out.attrs$dim, out.attrs$dimnames)
            res$se.fit <- array(res$se.fit, out.attrs$dim, out.attrs$dimnames)
        } else {
            res <- array(res, out.attrs$dim, out.attrs$dimnames)
        }
    }
    if (se){
        res$df <- object$one.delta^2/object$two.delta
    }

    res
}
