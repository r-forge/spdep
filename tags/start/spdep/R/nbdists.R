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


nbdists <- function(nb, coords) {
	if (class(nb) != "nb") 
        	stop("Not a neighbours list")
	if (!is.matrix(coords)) 
            stop("Data not in matrix form")
        if (any(is.na(coords))) 
            stop("Data include NAs")
	n.nb <- length(nb)
	np <- nrow(coords)
        if (np != n.nb) 
            stop("Number of coords not equal to number of regions")
        dimension <- ncol(coords)
        dlist <- .Call("nbdists", nb, as.matrix(coords), as.integer(np), 
            as.integer(dimension))
	attr(dlist[[1]], "call") <- match.call()
	invisible(dlist[[1]])
}

