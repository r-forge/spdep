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


nb2mat <- function(neighbours, glist=NULL, style="W", zero.policy=FALSE)
{
	if(class(neighbours) != "nb") stop("Not a neighbours list")
	listw <- nb2listw(neighbours, glist=glist, style=style,
		zero.policy=zero.policy)
	res <- listw2mat(listw)
	attr(res, "call") <- match.call()
	invisible(res)
}

listw2mat <- function(listw) {
	n <- length(listw$neighbours)
	cardnb <- card(listw$neighbours)
	if (any(is.na(unlist(listw$weights))))
		stop ("NAs in general weights list")
	res <- matrix(0, nrow=n, ncol=n)
	for (i in 1:n)
	    if (cardnb[i] > 0)
		res[i, listw$neighbours[[i]]] <- listw$weights[[i]]
	invisible(res)
}

invIrM <- function(neighbours, rho, glist=NULL, style="W") {
	if(class(neighbours) != "nb") stop("Not a neighbours list")
	n <- length(neighbours)
	V <- nb2mat(neighbours, glist, style)
	feasible <- 1/(range(eigen(V, only.values=TRUE)$values))
	if (rho <= feasible[1] || rho >= feasible[2])
		stop(paste("Rho outside feasible range:", feasible))
	mat <- diag(n) - rho * V
	res <- solve(mat)
	attr(res, "call") <- match.call()
	invisible(res)
}
