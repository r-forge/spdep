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

moran.plot <- function(x, listw, zero.policy=FALSE, labels=NULL, xlab=NULL, ylab=NULL, ...)
{
	if (class(listw) != "listw") stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	xname <- deparse(substitute(x))
	if (!is.numeric(x)) stop(paste(xname, "is not a numeric vector"))
	if (any(is.na(x))) stop("NA in X")
	n <- length(listw$neighbours)
	if (n != length(x)) stop("objects of different length")
	labs <- TRUE
	if (is.logical(labels) && !labels) labs <- FALSE
	if (is.null(labels) || length(labels) != n)
		labels <- as.character(1:n)
	wx <- lag.listw(listw, x, zero.policy)
	if (is.null(xlab)) xlab <- xname
	if (is.null(ylab)) ylab <- paste("spatially lagged", xname)
	plot(x, wx, xlab=xlab, ylab=ylab, ...)
	if (zero.policy)
		points(x[wx == 0.0], wx[wx == 0.0], col="orange", pch=19)
	xwx.lm <- lm(wx ~ x)
	abline(xwx.lm)
	abline(h=mean(wx), lty=2)
	abline(v=mean(x), lty=2)
	infl.xwx <- influence.measures(xwx.lm)
	is.inf <- which(apply(infl.xwx$is.inf, 1, any))
	points(x[is.inf], wx[is.inf], pch=9, col="red")
	if (labs)
	    text(x[is.inf], wx[is.inf], labels=labels[is.inf], pos=2, cex=0.7)
	rownames(infl.xwx$infmat) <- labels
	summary(infl.xwx)
	invisible(infl.xwx)
}

