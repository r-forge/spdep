# Copyright 2001-2002 by Roger Bivand
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


nblag <- function(neighbours, maxlag)
{
	if (!is.null(attributes(neighbours)$self.included) &&
		(as.logical(attributes(neighbours)$self.included)))
		stop("No lags for neighbours lists including self")
	n <- length(neighbours)
	lags <- vector(mode="list", length=maxlag)
	lags[[1]] <- neighbours
	cds <- card(neighbours)
	for (thislag in 2:maxlag)
		lags[[thislag]] <- vector(mode="list", length=n)
	for (i in 1:n) {
		already <- i
		new <- neighbours[[i]]
		for (thislag in 2:maxlag) {
			if (cds[i] > 0) {
				already <- c(already, new)
				active <- new
				new <- NULL
				for (j in active)
					new <- c(new, neighbours[[j]])
				new <- sort(unique(new))
				res <- new[-which(new %in% already)]
				if(length(res) == 0) lags[[thislag]][[i]] <- 0
				else lags[[thislag]][[i]] <- res
			}
			else lags[[thislag]][[i]] <- 0
		}
	}
	for (i in 1:maxlag) {
		class(lags[[i]]) <- "nb"
		attr(lags[[i]], "region.id") <- attr(neighbours, "region.id")
		lags[[i]] <- sym.attr.nb(lags[[i]])
	}
	attr(lags, "call") <- match.call()
	invisible(lags)
}
