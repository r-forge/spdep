\name{BH}
\alias{BH}
\alias{BH.polys}
\alias{BH.bbs}
\alias{BH.xy}
\alias{BH.nb}
\non_function{}
\title{Belo Horizont 1994 homicide data}
\usage{data(BH)}
\description{
The \code{BH} data frame has 81 rows and 4 columns for 81 planning units (UP) in the city. The dataset also includes \code{BH.polys} UP polygon boundaries, \code{BH.xy} UP polygon centroids, \code{BH.bbs} UP polygon bounding boxes, and \code{BH.nb} the neighbours list used in Assuncao and Reis (1999), 
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{UP}{UP identifiers}
    \item{ID}{numeric identifier}
    \item{POP}{UP population 1994}
    \item{HOMICIDES}{UP homicides counts 1994}
  }
}
\details{
 The dataset was originally used to motivate and test a new proposal for adjusting Moran's I when the variable being tested for spatial autocorrelation is a rate (see \code{EBImoran.mc()}}. The examples do not replicate the printed results because of an unfortunate re-ordering of the data - see the \code{moran.test()} with the re-ordered rate \code{p} to retrieve a printed result), but the substantive conclusions are sustained.
}
\source{
 Data obtained from paper authors directly.
}
\references{
Assuncao RM, Reis EA 1999 A new proposal to adjust Moran's I for population density. Statistics in Medicine 18, pp. 2147--2162
}
\examples{
data(BH)
plotpolys(BH.polys, BH.bbs)
plot(BH.nb, BH.xy, add=TRUE)
p <- BH$HOMICIDES / BH$POP
p.f <- ordered(cut(p*100000, breaks=c(0,1,6,15,40,625),
 include.lowest=TRUE))
plotpolys(BH.polys, BH.bbs,
 col=rev(c("black", heat.colors(4)))[codes(p.f)])
legend(c(618000, 624000), c(781000, 787000),
 legend=c("0", levels(p.f)[2:4], "625"),
 fill=rev(c("black", heat.colors(4))))
title("Homicide rates per 100,000, Belo Horizonte, 1994")
joincount.mc(p.f, nb2listw(BH.nb, style="B"), nsim=999)
moran.test(p, nb2listw(BH.nb, style="B"), randomisation=FALSE)
moran.test(p[order(BH$UP)], nb2listw(BH.nb, style="B"),
 randomisation=FALSE)
moran.test(p, nb2listw(BH.nb, style="B"))
moran.mc(p, nb2listw(BH.nb, style="B"), nsim=999)
EBImoran.mc(BH$HOMICIDES, BH$POP, nb2listw(BH.nb, style="B"), nsim=999)
outlier <- which(p*100000 == 625)
BH[outlier,]
plot(log10(BH$POP), BH$HOMICIDES,
 pch=ifelse(p*100000 == 625, 1, 18), xaxt="n",
 xlab="UP population", ylab="Homicides",
 main="Homicides, Belo Horizonte, 1994")
axis(1, at=c(2,3,4), labels=c(100,1000, 10000))
text(log10(BH$POP)[outlier], BH$HOMICIDES[outlier],
 pos=1, adj=0.5, paste("UP", BH$UP[outlier]))
moran.mc(p[-outlier],
 nb2listw(subset(BH.nb, !(attr(BH.nb, "region.id") == BH$UP[outlier])),
 style="B"), nsim=999)
EBImoran.mc(BH$HOMICIDES[-outlier], BH$POP[-outlier],
 nb2listw(subset(BH.nb, !(attr(BH.nb, "region.id") == BH$UP[outlier])),
 style="B"), nsim=999)
}
\keyword{datasets}
