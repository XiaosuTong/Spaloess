#' Spatially Local Polynomial Regression Fitting
#'
#' The first layer of the Spatial locally weighted regression, using local fitting with different 
#' type of distance calculation.
#'
#' @useDynLib Spaloess
#'
#' @param formula 
#'     a formula specifying the numeric response and one to four numeric predictors.
#' @param data 
#'     an optional data from, list or environment containing the variables in the model. If 
#'     not found in 'data', the variables are taken from 'environment', typically the environment 
#'     from which 'loess' is called.
#' @param weights 
#'     optional weights for each case
#' @param subset 
#'     an optional specification of a subset of the data to be used
#' @param na.action 
#'     the action to be taken with missing values in the response or predictors.  The default is 
#'     given by 'getOption("na.action")'.
#' @param distance 
#'     Options: "Euclid", or "Latlong" which is for great circle distance
#' @param model 
#'     Should the model frame be returned?
#' @param span 
#'     The parameter alpha which controls the portion of data points used in the local fit.
#' @param enp.target 
#'     An alternative way to specify 'span', as the approximate equivalent number of parameters to 
#'     be used.
#' @param degree 
#'     The degree of the polynomials to be used, normally 1 or 2. (Degree 0 is also allowed, but see 
#'     the 'Note'.)
#' @param parametric 
#'     should any terms be fitted globally rather than locally? Terms can be specified by name, 
#'     number or as a logical vector of the same length as the number of predictors.
#' @param drop.square 
#'     For fits with more than one predictor and 'degree = 2', should the quadratic term be dropped 
#'     for particular predictors?  Terms are specified in the same way as for 'parametric'.
#' @param normalize 
#'     Should the predictors be normalized to a common scale if there is more than one?  The 
#'     normalization used is to set the 10% trimmed standard deviation to one. Set to false for 
#'     "Latlong" distance.
#' @param family 
#'     If 'gaussian' fitting is by least-squares, and if 'symmetric' a re-descending M estimator is 
#'     used with Tukey's bi-weight function.
#' @param method 
#'     Fit the model or just extract the model frame.
#' @param napred
#'     Should missing observations in the dataset be predicted. Default is TRUE.
#' @param control 
#'     control parameters: see 'loess.control'.
#' @param ...
#'     arguments passed to or from other methods.
#' @details
#'     This spaloess function is the first wrapper of the spatial loess fitting procedure. It checks
#'     all the validity of all input arguments, and formats arguments like drop,square, parametric.
#'     Also generate other important arguments, like iteration, and pass all arguments into the 
#'     second wrapper function: newsimpleLoess 
#' @author 
#'     Xiaosu Tong, based on 'loess' function of B. D. Ripley, and 'cloess' package of Cleveland,
#'     Grosse and Shyu.  
#' @export
#' @examples
#'     set.seed(66)
#'     x1 <- rnorm(100, mean=-100, sd=10)
#'     x2 <- rnorm(100, mean=38, sd=4)
#'     y <- 0.1*x1 + 1*x2 - 10 + rnorm(100, 0, 1.3)
#'     testdata <- data.frame(LON = x1, LAT = x2, tmax = y)
#'     cars.lo <- spaloess(tmax ~ LON + LAT, testdata, distance = "Latlong")


spaloess <- function (formula, data, weights, subset, na.action, model = FALSE, napred = TRUE, 
    span = 0.75, enp.target, degree = 2L, parametric = FALSE, distance = "Latlong",
    drop.square = FALSE, normalize = FALSE, family = c("gaussian", 
        "symmetric"), method = c("loess", "model.frame"), control = loess.control(...), 
    ...) 
{
  family <- match.arg(family)
  method <- match.arg(method)

  ## returns a call in which all of the specified arguments are specified by their full names.
  mf <- match.call(expand.dots = FALSE)
  mf$model <- mf$span <- mf$enp.target <- mf$degree <- mf$parametric <- mf$distance <- mf$napred <- NULL
  mf$drop.square <- mf$normalize <- mf$family <- mf$method <- mf$control <- mf$... <- NULL
  
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  if (match.arg(method) == "model.frame") return(mf)
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  w <- model.weights(mf)

  ## Set the weights
  if (is.null(w)) w <- rep(1, length(y))

  ## nmx is the independent variables name.
  nmx <- as.character(attr(mt, "variables"))[-(1L:2)]
  ## nmy is the response variable name
  nmy <- as.character(attr(mt, "variables"))[2]

  x <- mf[, nmx, drop = FALSE]

  ## get the locations of NAs
  if(napred) {
    index <- is.na(data[,nmy])
    xna <- data[index, nmx]
  }


  if (any(sapply(x, is.factor))) stop("predictors must all be numeric")
  
  if (distance == "Latlong" | distance == "L") {
    if("la" %in% tolower(substr(nmx, 1, 2)) & "lo" %in% tolower(substr(nmx, 1, 2))) {
      for(ii in c("la","lo")) {
        indx <- grep(ii, tolower(substr(nmx, 1, 2)))
        x[, indx] <- 2 * pi * x[, indx]/360
      }
    } else {
      stop("predictors must be longitude and latitude for great circle distance") 
    }
  }
      
  x <- as.matrix(x)
  D <- ncol(x) 
  nmx <- colnames(x)
  names(nmx) <- nmx

  drop.square <- match(nmx, nmx[drop.square], 0L) > 0L 
  parametric <- match(nmx, nmx[parametric], 0L) > 0L 
  
  if (!match(degree, 0L:2L, 0L)) stop("'degree' must be 0, 1 or 2")
  
  ## iteration times for robust fit. If family is gaussian, no robust fit
  iterations <- if (family == "gaussian") 1 else control$iterations

  ## calculate the span argument if enp.target is specified
  if (!missing(enp.target)) { 
    if (!missing(span)) 
        warning("both 'span' and 'enp.target' specified: 'span' will be used")
    else {
        tau <- switch(degree + 1L, 1, D + 1, (D + 1) * (D + 
          2)/2) - sum(drop.square)
        span <- 1.2 * tau/enp.target #calculate the span based on enp.target 
    }
  }
  ## Check the control arguments
  if (!is.list(control) || !is.character(control$surface) || 
      !is.character(control$statistics) || !is.character(control$trace.hat) || 
      !is.numeric(control$cell) || !is.numeric(iterations)) { 
    stop("invalid 'control' argument")
  } 
 
 ## After checking all arguments' validity, pass all valid arguments into newsimpleLoess function. 
  fit <- newsimpleLoess(
    y, x, w, span, degree, distance, parametric, drop.square,
    normalize, control$statistics, control$surface, control$cell, 
    iterations, control$trace.hat
  )

  fit$call <- match.call()
  fit$terms <- mt
  fit$xnames <- nmx
  fit$x <- x
  fit$y <- y
  fit$weights <- w

  if (model) fit$model <- mf 
  fit$na.action <- attr(mf, "na.action")
  
  if (napred) {
    naPrediction <- predloess(fit, newdata = xna)
    fit$pred <- cbind(xna, fitted = naPrediction)
  } else {
    fit$pred <- NULL
  }
  rownames(fit$pred) <- NULL

  fit

}

