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

localG <- function(x, listw, zero.policy=FALSE) {
	if (class(listw) != "listw")
		stop(paste(deparse(substitute(listw)), "is not a listw object"))
	if (!is.numeric(x))
		stop(paste(deparse(substitute(x)), "is not a numeric vector"))
	if (any(is.na(x))) stop(paste("NA in ", deparse(substitute(x))))
	n <- length(listw$neighbours)
	if (n != length(x))stop("Different numbers of observations")
	gstari <- FALSE
	if (!is.null(attr(listw$neighbours, "self.included")) &&
		attr(listw$neighbours, "self.included")) gstari <- TRUE
	lx <- lag.listw(listw, x, zero.policy=zero.policy)
	if (gstari) {
		xibar <- rep(mean(x), n)
		si2 <- rep(sum(scale(x, scale=FALSE)^2)/n, n)
	} else {
		xibar <- (rep(sum(x),n) - x) / (n - 1)
		si2 <- ((rep(sum(x^2),n) - x^2) / (n - 1)) - xibar^2
	}
	Wi <- sapply(listw$weights, sum)
	S1i <- sapply(listw$weights, function(x) sum(x^2))
	res <- (lx - Wi*xibar)
	if (gstari) {
		res <- res / sqrt(si2*((n*S1i - Wi^2)/(n-1)))
	} else {
		res <- res / sqrt(si2*(((n-1)*S1i - Wi^2)/(n-2)))
	}
	attr(res, "gstari") <- gstari
	attr(res, "call") <- match.call()
	class(res) <- "localG"
	invisible(res)
}
