\name{my.simple1}
\alias{my.simple1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
my.simple1(y, x, weights, span = 0.75, degree = 2L, parametric = FALSE, drop.square = FALSE, normalize = TRUE, statistics = "approximate", surface = "interpolate", cell = 0.2, iterations = 1L, trace.hat = "exact")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{weights}{
%%     ~~Describe \code{weights} here~~
}
  \item{span}{
%%     ~~Describe \code{span} here~~
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
  \item{statistics}{
%%     ~~Describe \code{statistics} here~~
}
  \item{surface}{
%%     ~~Describe \code{surface} here~~
}
  \item{cell}{
%%     ~~Describe \code{cell} here~~
}
  \item{iterations}{
%%     ~~Describe \code{iterations} here~~
}
  \item{trace.hat}{
%%     ~~Describe \code{trace.hat} here~~
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
function (y, x, weights, span = 0.75, degree = 2L, parametric = FALSE, 
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
    s <- NULL
    pars <- list(robust = robust, span = span, degree = degree, 
        normalize = normalize, parametric = parametric, drop.square = drop.square, 
        surface = surface, cell = cell, family = if (iterations <= 
            1L) "gaussian" else "symmetric", iterations = iterations)
    fit <- list(n = N, s = s, divisor = divisor)
    fit$pars <- pars
    if (surface == "interpolate") 
        fit$kd <- fit.kd
    class(fit) <- "loess"
    fit
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
