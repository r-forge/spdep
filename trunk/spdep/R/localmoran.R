# Copyright 2001 by Roger Bivand 
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#

localmoran <- function(x, listw, zero.policy=FALSE)
{
	if (class(listw) != "listw")
		stop(paste(deparse(substitute(listw)), "is not a listw object"))
	if (!is.null(attr(listw$neighbours, "self.included")) &&
		attr(listw$neighbours, "self.included"))
		stop("Self included among neighbours")
	n <- length(listw$neighbours)
	if (!is.numeric(x))
		stop(paste(deparse(substitute(x)), "is not a numeric vector"))
	if (any(is.na(x))) stop(paste("NA in ", deparse(substitute(x))))
	if (n != length(x))stop("Different numbers of observations")
	res <- data.frame(matrix(nrow=n,ncol=4))
	colnames(res) <- c("Ii", "E.Ii", "Var.Ii", "Z.Ii")
	z <- scale(x, scale=F)
	lz <- lag.listw(listw, z, zero.policy=zero.policy)
	s2 <- sum(z^2)/n
	res[,1] <- (z/s2) * lz
	Wi <- sapply(listw$weights, sum)
	res[,2] <- -Wi / (n-1)
	b2 <- (sum(z^4)/n)/(s2^2)
	A <- (n-b2) / (n-1)
	B <- (2*b2 - n) / ((n-1)*(n-2))
	C <- Wi^2 / ((n-1)^2)
	Wi2 <- sapply(listw$weights, function(x) sum(x^2))
	Wikh2 <- sapply(listw$weights, function(x) {1 - crossprod(x,x)})
	res[,3] <- A*Wi2 + B*Wikh2 - C
	res[,4] <- (res[,1] - res[,2]) / sqrt(res[,3])
	attr(res, "call") <- match.call()
	class(res) <- c(class(res), "localmoran")
	invisible(res)
}


