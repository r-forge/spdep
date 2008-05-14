# Copyright 2001-6 by Roger Bivand 
#

droplinks <- function(nb, drop, sym=TRUE) {
	if (class(nb) != "nb") stop("Not a neighbours list")
	n <- length(nb)
	cnb <- card(nb)
	if (n < 1) stop("non-positive length of nb")
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
		if (sym && cnb[i] > 0) {
			for (j in nb[[i]])
				nb[[j]] <- nb[[j]][nb[[j]] != i]
		}
		nb[[i]] <- as.integer(0)
	}
	nb <- sym.attr.nb(nb)
	nb
}
