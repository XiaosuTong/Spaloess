\name{my.loess1}
\alias{my.loess1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
my.loess1(formula, data, weights, subset, na.action, model = FALSE, span = 0.75, enp.target, degree = 2L, parametric = FALSE, drop.square = FALSE, normalize = TRUE, family = c("gaussian", "symmetric"), method = c("loess", "model.frame"), control = loess.control(...), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{weights}{
%%     ~~Describe \code{weights} here~~
}
  \item{subset}{
%%     ~~Describe \code{subset} here~~
}
  \item{na.action}{
%%     ~~Describe \code{na.action} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{span}{
%%     ~~Describe \code{span} here~~
}
  \item{enp.target}{
%%     ~~Describe \code{enp.target} here~~
}
  \item{degree}{
%%     ~~Describe \code{degree} here~~
}
  \item{parametric}{
%%     ~~Describe \code{parametric} here~~
}
  \item{drop.square}{
%%     ~~Describe \code{drop.square} here~~
}
  \item{normalize}{
%%     ~~Describe \code{normalize} here~~
}
  \item{family}{
%%     ~~Describe \code{family} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{control}{
%%     ~~Describe \code{control} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (formula, data, weights, subset, na.action, model = FALSE, 
    span = 0.75, enp.target, degree = 2L, parametric = FALSE, 
    drop.square = FALSE, normalize = TRUE, family = c("gaussian", 
        "symmetric"), method = c("loess", "model.frame"), control = loess.control(...), 
    ...) 
{
    family <- match.arg(family)
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
            span <- 1.2 * tau/enp.target
        }
    if (!is.list(control) || !is.character(control$surface) || 
        !is.character(control$statistics) || !is.character(control$trace.hat) || 
        !is.numeric(control$cell) || !is.numeric(iterations)) 
        stop("invalid 'control' argument")
    fit <- my.simple1(y, x, w, span, degree, parametric, drop.square, 
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
