% Copyright 2001 by Roger S. Bivand
\name{listw2sn}
\alias{listw2sn}
\title{Spatial neighbour sparse representation}
\description{
The function makes a \code{spatial neighbour} object representation (similar to the S-PLUS spatial statististics module representation of a \code{listw} spatial weights object.
}
\usage{
listw2sn(listw)
}
\arguments{
  \item{listw}{a \code{listw} object from for example \code{nb2listw}}
}
\value{
A data frame with three columns, and with class \code{spatial.neighbour}:
  \item{from}{region number id for the start of the link (S-PLUS row.id)}
  \item{to}{region number id for the end of the link (S-PLUS col.id)}
  \item{weights}{weight for this link}
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{nb2listw}}}

\examples{
data(columbus)
col.listw <- nb2listw(col.gal.nb, style="C")
col.listw$neighbours[[1]]
col.listw$weights[[1]]
col.sn <- listw2sn(col.listw)
col.sn[col.sn[,1] == 1,]
}
\keyword{spatial}
