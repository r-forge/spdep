% Copyright 2001-3 by Roger S. Bivand
\name{lag.listw}
\alias{lag.listw}
\title{Spatial lag of a numeric vector}
\description{
 Using a \code{listw} sparse representation of a spatial weights matrix, compute the lag vector \eqn{V x}
}
\usage{
lag.listw(x, var, zero.policy=FALSE, NAOK=FALSE, ...)
}

\arguments{
  \item{x}{a \code{listw} object created for example by \code{nb2listw}}
  \item{var}{a numeric vector the same length as the neighbours list in listw}
  \item{zero.policy}{if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}
  \item{NAOK}{if 'TRUE' then any 'NA' or 'NaN' or 'Inf' values in var are passed on to the foreign function.  If 'FALSE', the presence of 'NA' or 'NaN' or 'Inf' values is regarded as an error.}
  \item{...}{additional arguments}
}
\value{
  a numeric vector the same length as var
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{nb2listw}}}

\examples{
data(oldcol)
Vx <- lag.listw(nb2listw(COL.nb, style="W"), COL.OLD$CRIME)
plot(Vx, COL.OLD$CRIME)
plot(ecdf(COL.OLD$CRIME))
plot(ecdf(Vx), add=TRUE, col.points="red", col.hor="red")
}
\keyword{spatial}