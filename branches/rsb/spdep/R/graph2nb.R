# Copyright 2001 by Nicholas Lewin-Koh
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


graph2nb <- function(gob, row.names=NULL,sym=FALSE) {
	if (class(gob) != "Graph") stop("Not a Graph object")
	res <- vector(mode="list", length=gob$np)
    	if (!is.null(row.names)) if(length(row.names) != gob$np)
        	stop("row.names wrong length")
    	if (is.null(row.names)) row.names <- as.character(1:gob$np)
        if(sym){
          for (i in 1:gob$np) {
		res[[i]] <- sort(unique(c(gob$to[gob$from==i],
                                       gob$from[gob$to==i])))
	  	if(length(res[[i]]) == 0) res[[i]] <- 0
	  }
        }
        else{
	  for (i in 1:gob$np) {
		res[[i]] <- sort(gob$to[gob$from==i])
	  	if(length(res[[i]]) == 0) res[[i]] <- 0
	  }
        }
        attr(res, "region.id") <- row.names
 	attr(res, "call") <- attr(gob, "call")
 	attr(res, "type") <- attr(gob, "type")
	class(res) <- "nb"
	res <- sym.attr.nb(res)
	invisible(res)
}
