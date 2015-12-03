\name{my.predict.loess}
\alias{my.predict.loess}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
my.predict.loess(object, newdata = NULL, se = FALSE, na.action = na.pass, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{newdata}{
%%     ~~Describe \code{newdata} here~~
}
  \item{se}{
%%     ~~Describe \code{se} here~~
}
  \item{na.action}{
%%     ~~Describe \code{na.action} here~~
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
function (object, newdata = NULL, se = FALSE, na.action = na.pass, 
    ...) 
{
    if (!inherits(object, "loess")) 
        stop("first argument must be a \"loess\" object")
    if (is.null(newdata) && !se) 
        return(fitted(object))
    newx <- if (is.null(newdata)) 
        object$x
    else if (is.data.frame(newdata)) 
        as.matrix(model.frame(delete.response(terms(object)), 
            newdata, na.action = na.action))
    else as.matrix(newdata)
    res <- with(object, my.predLoess(y, x, newx, s, weights, 
        pars$robust, pars$span, pars$degree, pars$normalize, 
        pars$parametric, pars$drop.square, pars$surface, pars$cell, 
        pars$family, kd, divisor, se = se))
    if (!is.null(out.attrs <- attr(newdata, "out.attrs"))) {
        if (se) {
            res$fit <- array(res$fit, out.attrs$dim, out.attrs$dimnames)
            res$se.fit <- array(res$se.fit, out.attrs$dim, out.attrs$dimnames)
        }
        else res <- array(res, out.attrs$dim, out.attrs$dimnames)
    }
    if (se) 
        res$df <- object$one.delta^2/object$two.delta
    res
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
