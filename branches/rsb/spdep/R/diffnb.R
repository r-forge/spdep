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


diffnb <- function(x, y, verbose=TRUE) {
	if (class(x) != "nb") stop("not a neighbours list")
	if (class(y) != "nb") stop("not a neighbours list")
	n <- length(x)
	if(n != length(y)) stop("lengths differ")
	res <- vector(mode="list", length=n)
	for (i in 1:n) {
		xi <- x[[i]]
		yi <- y[[i]]
		xt <- xi %in% yi
		yt <- yi %in% xi
		if (!(all(xt) && all(yt))) {
			res[[i]] <- sort(unique(c(xi[which(!xt)],
				yi[which(!yt)])))
			if(verbose)
				cat("Neighbour difference for region:",
				i, "in relation to:", res[[i]], "\n")
		}
	}
	class(res) <- "nb"
	attr(res, "region.id") <- attr(x, "region.id")
	attr(res, "call") <- match.call()
	invisible(res)
}	
	
