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

droplinks <- function(nb, drop, sym=TRUE) {
	if (class(nb) != "nb") stop("Not a neighbours list")
	n <- length(nb)
	if (is.logical(drop)) {
		if(length(drop) != n) stop("Argument lengths differ")
		idrop <- which(drop == TRUE)
	} else if(is.character(drop)) {
		row.names <- as.character(attr(nb, "region.id"))
		idrop <- match(drop, row.names)
		if(any(is.na(idrop))) stop("Region to drop not found")
	} else {
		idrop <- match(drop, 1:n)
		if(any(is.na(idrop))) stop("Region to drop not found")
	}
	if((attr(nb, "sym") == FALSE) && (sym == TRUE)) {
		warning("setting sym to FALSE")
		sym <- FALSE
	}
	for (i in idrop) {
		if (sym) {
			for (j in nb[[i]])
				nb[[j]] <- nb[[j]][nb[[j]] != i]
		}
		nb[[i]] <- 0
	}
	invisible(nb)
}

