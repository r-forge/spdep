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


tri2nb <- function(coords, row.names = NULL) {
	require(tripack)
	n <- nrow(coords)
    	if (!is.null(row.names)) if(length(row.names) != n)
        	stop("row.names wrong length")
    	if (is.null(row.names)) row.names <- as.character(1:n)
	tri <- tri.mesh(x=coords[,1], y=coords[,2])
	nb <- neighbours(tri)
 	attr(nb, "region.id") <- row.names
	class(nb) <- "nb"
	attr(nb, "tri") <- TRUE
	attr(nb, "call") <- match.call()
	nb <- sym.attr.nb(nb)
	invisible(nb)
}

