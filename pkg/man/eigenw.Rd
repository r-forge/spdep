% Copyright 2002-2011 by Roger S. Bivand
\name{eigenw}
\alias{eigenw}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Spatial weights matrix eigenvalues}
\description{
  The function returns a numeric vector of eigenvalues of the weights matrix
generated from the spatial weights object \code{listw}. The eigenvalues are 
used to speed the computation of the Jacobian in spatial SAR model estimation:

\deqn{\log(\det[I - \rho W]) = \log\prod_{i=1}^{n}(1 - \rho \lambda_i)}

where \eqn{W}{W} is the n by n spatial weights matrix, and \eqn{\lambda_i}{lambda[i]} are the
eigenvalues of \eqn{W}{W}.
}
\usage{
eigenw(listw, quiet=NULL)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{quiet}{default NULL, use global !verbose option value; set to FALSE for short summary}
}

\value{
  a numeric vector of eigenvalues of the weights matrix generated from the
spatial weights object \code{listw}.
}
\references{Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 155;
Ord, J. K. 1975 Estimation methods for models of spatial interaction, Journal
of the American Statistical Association, 70, 120-126.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{eigen}}%, \code{\link{logSpwdet}}
}

\examples{
data(oldcol)
W.eig <- eigenw(nb2listw(COL.nb, style="W"))
1/range(W.eig)
S.eig <- eigenw(nb2listw(COL.nb, style="S"))
1/range(S.eig)
B.eig <- eigenw(nb2listw(COL.nb, style="B"))
1/range(B.eig)
crds <- cbind(COL.OLD$X, COL.OLD$Y)
k6 <- knn2nb(knearneigh(crds, k=6))
is.symmetric.nb(k6)
k6eig <- eigenw(nb2listw(k6, style="W"))
is.complex(k6eig)
rho <- 0.5
Jc <- sum(log(1 - rho * k6eig))
Jc
W <- as(as_dgRMatrix_listw(nb2listw(k6, style="W")), "CsparseMatrix")
I <- diag(length(k6))
Jl <- sum(log(abs(diag(slot(lu(I - rho * W), "U")))))
Jl
Jr <- sum(log(1 - rho * Re(k6eig)))
Jr
Rev <- Re(k6eig)[which(Im(k6eig) == 0)]
Cev <- k6eig[which(Im(k6eig) != 0)]
pCev <- Cev[Im(Cev) > 0]
}
\keyword{spatial}
