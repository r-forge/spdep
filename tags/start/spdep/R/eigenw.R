# Copyright 2002 by Roger Bivand 
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

eigenw <- function(listw, quiet=TRUE)
{
	if(class(listw) != "listw") stop("not a listw object")
	w <- listw2mat(listw)
	e <- eigen(w, only.values=TRUE)$values
	if (is.complex(e)) e <- Re(e)
	if (!quiet) {
		cat("Largest eigenvalue:", max(e),
		"Sum of eigenvalues:", sum(e), "\n")
	}
	e
}

