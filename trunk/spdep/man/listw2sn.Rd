% Copyright 2001-2 by Roger S. Bivand
\name{listw2sn}
\alias{listw2sn}
\alias{spwdet}
\alias{log.spwdet}
\title{Spatial neighbour sparse representation}
\description{
The function makes a \code{spatial neighbour} object representation (similar to the S-PLUS spatial statististics module representation of a \code{listw} spatial weights object. The object is used in calculating the Jacobian of the likelihood functions of spatial autoregressive models when the "sparse" method is chosen in \code{errorsarlm} and \code{lagsarlm}. 
}
\usage{
listw2sn(listw)
spwdet(sparseweights, rho)
log.spwdet(sparseweights, rho)

}
\arguments{
  \item{listw}{a \code{listw} object from for example \code{nb2listw}}
  \item{sparseweights}{a \code{spatial.neighbour} object}
  \item{rho}{spatial autoregressive parameter value}
}

\details{sparse, B&P}

\value{
A data frame with three columns, and with class \code{spatial.neighbour}:
  \item{from}{region number id for the start of the link (S-PLUS row.id)}
  \item{to}{region number id for the end of the link (S-PLUS col.id)}
  \item{weights}{weight for this link}
\code{log.spwdet} returns log det(I - rho * W).
}

\references{sparse, B&P}

\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{nb2listw}}, \code{\link{errorsarlm}},
 \code{\link{lagsarlm}}}

\examples{
data(columbus)
col.listw <- nb2listw(col.gal.nb)
col.listw$neighbours[[1]]
col.listw$weights[[1]]
col.sn <- listw2sn(col.listw)
col.sn[col.sn[,1] == 1,]
rho <- seq(-0.8, 0.9, 0.1)
for (i in rho) print(paste("rho:", i, "log(det(I - rho*W))",
  log.spwdet(col.sn, i)), quote=TRUE)
}
\keyword{spatial}
