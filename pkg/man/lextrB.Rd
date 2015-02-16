\name{lextrB}
\alias{lextrB}
\alias{lmaxB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find extreme eigenvalues of binary symmetric spatial weights}
\description{The functions find extreme eigenvalues of binary symmetric spatial weights, particularly when these form planar graphs. \code{lmaxB} finds the largest eigenvalue using Rayleigh quotient methods. \code{lextrB} first calls \code{lmaxB}, and uses its output to find the smallest eigenvalue in addition.}
\usage{
lextrB(lw, zero.policy = TRUE, control = list())
lmaxB(lw, zero.policy=TRUE, control=list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lw}{a binary symmetric \code{listw} object from, for example, \code{nb2listw} with style \dQuote{B}}
  \item{zero.policy}{default NULL, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}
  \item{control}{a list of control arguments}
}
\section{Control arguments}{
\describe{
  \item{trace}{report values in while loops, default NULL assuming FALSE; logical}
  \item{tol}{tolerance for breaking while loops, default \code{.Machine$double.eps^(1/2)}; numeric}
  \item{maxiter}{maximum number of iterations in while loops, default \code{6 * (length(lw$neighbours) - 2}; integer}
  \item{useC}{use C code, default TRUE, logical (only in \code{lextrB})}
}
}
\value{
The functions return approximations to the extreme eigenvalues with the eigenvectors returned as attributes of this object.
}
\references{Griffith, D. A. (2004). Extreme eigenfunctions of adjacency matrices for planar graphs employed in spatial analyses. \emph{Linear Algebra and its Applications}, 388:201–219.}
\author{Roger Bivand, Yongwan Chun, Daniel Griffith}

\note{It may be necessary to modify control arguments if warnings about lack of convergence are seen.}

\examples{
data(boston)
ab.listb <- nb2listw(boston.soi, style="B")
er <- range(eigenw(ab.listb))
er
res_1 <- lextrB(ab.listb)
c(res_1)
if (require(igraph)) {
  B <- as(ab.listb, "symmetricMatrix")
  n <- length(boston.soi)
  f2 <- function(x, extra=NULL) {as.vector(B \%*\% x)}
  ar1 <- arpack(f2, sym=TRUE, options=list(n=n, nev=2, ncv=8,
    which="BE", maxiter=200))
  print(ar1$values)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}
