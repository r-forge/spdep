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


plot.nb <- function(x, coords, col="black", points=TRUE, add=FALSE, ...) {
	nb <- x
	x <- coords[,1]
	y <- coords[,2]
	n <- length(nb)
	xlim <- range(x)
	ylim <- range(y)
	if (!add) {
		plot.new()
        	plot.window(xlim = xlim, ylim = ylim, log="", asp=1)
	}
	for (i in 1:n) {
        	inb <- nb[[i]]
        	for (j in inb)
			lines(c(x[i], x[j]), c(y[i], y[j]),
				col=col, ...)
	}
	if (points) points(x, y, ...)
}
