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

listw2sn <- function(listw) {
	if(class(listw) != "listw") stop("not a listw object")
	z <- .Call("listw2sn", listw$neighbours, listw$weights,
		PACKAGE="spdep")
	res <- as.data.frame(list(from=z[[1]], to=z[[2]], weights=z[[3]]))
	class(res) <- c(class(res), "spatial.neighbour")
	attr(res, "region.id") <- attr(listw, "region.id")
	neighbours.attrs <- names(attributes(listw$neighbours))
	attr(res, "neighbours.attrs") <- neighbours.attrs
	weights.attrs <- names(attributes(listw$weights))
	attr(res, "weights.attrs") <- weights.attrs
	attr(res, "listw.call") <- attr(listw, "call")
	invisible(res)
}
