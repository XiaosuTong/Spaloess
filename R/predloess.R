##############################################
##predict loess function
##############################################


#' My Loess Function Title
#'
#' My Loess details
#' @param formula FORMULA DESCRIPTION. DO THIS FOR EACH ONE
#' @param data
#' @param weights
#' @param subset
#' @param na.action
#' @param model
#' @param span
#' @param enp.target
#' @param degree
#' @param parametric
#' @param drop.square
#' @param normalize
#' @param family
#' @param method
#' @param control
#' @param distance
#' @param ...
#' @author Xiaosu Tong
#' @export
#' @examples
#'   a <- my.loess1(dist ~ speed, cars, control = loess.control(surface = "direct"))
#'   str(a)
#'   my.predict.loess(a, data.frame(speed = seq(5, 30, 1)))
predictLoessf <- function (object, newdata = NULL, se = FALSE, na.action = na.pass, ...)
{
    if (!inherits(object, "loess")) {
        stop("first argument must be a \"loess\" object")
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
