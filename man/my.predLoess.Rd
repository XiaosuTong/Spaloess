\name{my.predLoess}
\alias{my.predLoess}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
my.predLoess(y, x, newx, s, weights, robust, span, degree, normalize, parametric, drop.square, surface, cell, family, kd, divisor, se = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{newx}{
%%     ~~Describe \code{newx} here~~
}
  \item{s}{
%%     ~~Describe \code{s} here~~
}
  \item{weights}{
%%     ~~Describe \code{weights} here~~
}
  \item{robust}{
%%     ~~Describe \code{robust} here~~
}
  \item{span}{
%%     ~~Describe \code{span} here~~
}
  \item{degree}{
%%     ~~Describe \code{degree} here~~
}
  \item{normalize}{
%%     ~~Describe \code{normalize} here~~
}
  \item{parametric}{
%%     ~~Describe \code{parametric} here~~
}
  \item{drop.square}{
%%     ~~Describe \code{drop.square} here~~
}
  \item{surface}{
%%     ~~Describe \code{surface} here~~
}
  \item{cell}{
%%     ~~Describe \code{cell} here~~
}
  \item{family}{
%%     ~~Describe \code{family} here~~
}
  \item{kd}{
%%     ~~Describe \code{kd} here~~
}
  \item{divisor}{
%%     ~~Describe \code{divisor} here~~
}
  \item{se}{
%%     ~~Describe \code{se} here~~
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
function (y, x, newx, s, weights, robust, span, degree, normalize, 
    parametric, drop.square, surface, cell, family, kd, divisor, 
    se = FALSE) 
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
                fit = double(M), L = double(N * M))[c("fit", 
                "L")]
            fit[!nas] <- z$fit
            ses <- (matrix(z$L^2, M, N)/rep(weights, rep(M, N))) \%*\% 
                rep(1, N)
            se.fit[!nas] <- drop(s * sqrt(ses))
        }
        else {
            fit[!nas] <- .C("loess_dfit", y, x, as.double(x.evaluate), 
                as.double(weights * robust), as.double(span), 
                as.integer(degree), as.integer(nonparametric), 
                as.integer(order.drop.sqr), as.integer(sum.drop.sqr), 
                as.integer(D), as.integer(N), as.integer(M), 
                fit = double(M))$fit
        }
    }
    else {
        inside <- matrix(FALSE, M, ncol = D)
        ranges <- apply(x, 2L, range)
        inside <- (x.evaluate <= rep(ranges[2L, ], rep(M, D))) & 
            (x.evaluate >= rep(ranges[1L, ], rep(M, D)))
        inside <- inside \%*\% rep(1, D) == D
        inside[is.na(inside)] <- FALSE
        M1 <- sum(inside)
        fit <- rep(NA_real_, M)
        if (any(inside)) 
            fit[inside] <- .C("loess_ifit", as.integer(kd$parameter), 
                as.integer(kd$a), as.double(kd$xi), as.double(kd$vert), 
                as.double(kd$vval), as.integer(M1), as.double(x.evaluate[inside, 
                  ]), fit = double(M1))$fit
        if (se) {
            se.fit <- rep(NA_real_, M)
            if (any(inside)) {
                L <- .C("loess_ise", y, x, as.double(x.evaluate[inside, 
                  ]), as.double(weights), as.double(span), as.integer(degree), 
                  as.integer(nonparametric), as.integer(order.drop.sqr), 
                  as.integer(sum.drop.sqr), as.double(span * 
                    cell), as.integer(D), as.integer(N), as.integer(M1), 
                  double(M1), L = double(N * M1))$L
                tmp <- (matrix(L^2, M1, N)/rep(weights, rep(M1, 
                  N))) \%*\% rep(1, N)
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
