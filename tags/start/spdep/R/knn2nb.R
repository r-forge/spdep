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


knn2nb <- function(knn, row.names=NULL, sym=FALSE) {
	if (class(knn) != "knn") stop("Not a knn object")
	res <- vector(mode="list", length=knn$np)
    	if (!is.null(row.names)) if(length(row.names) != knn$np)
        	stop("row.names wrong length")
    	if (is.null(row.names)) row.names <- as.character(1:knn$np)
        if(sym){
          to<-as.vector(knn$nn)
          from<-rep(1:knn$np,knn$k)
          for (i in 1:knn$np)res[[i]] <- sort(unique(c(to[from==i],
                                                       from[to==i]))) 
        } else {
          for (i in 1:knn$np) res[[i]] <- sort(knn$nn[i,])
        }
 	attr(res, "region.id") <- row.names
 	attr(res, "call") <- attr(knn, "call")
        attr(res, "sym") <- sym
	attr(res, "type") <- "knn"
 	attr(res, "knn-k") <- knn$k
	class(res) <- "nb"
	invisible(res)
}
