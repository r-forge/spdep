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

relativeneigh <- function(coords) {
    x <- coords
    if (!is.matrix(x)) stop("Data not in matrix form")
    if (any(is.na(x))) stop("Data cannot include NAs")
    np <- nrow(x)
    if(ncol(x)!=2) stop("Only works in 2d")
    g1<-g2<-rep(0,np*3)
    nogab <- 0
    z <- .C("compute_relative", np=as.integer(np), from=as.integer(g1),
             to=as.integer(g2), nedges=as.integer(nogab),
             x=as.double(x[,1]), y=as.double(x[,2]))
    z$from<-z$from[1:z$nedges]
    z$to<-z$to[1:z$nedges]
    attr(z, "call") <- match.call()
    class(z)<-c("Graph","relative")
    invisible(z)
}

plot.relative<-function(x, show.points=FALSE, add=FALSE, linecol=par(col),...){
  if(!add) plot(rel$x,rel$y,type='n')
  segments(rel$x[rel$from],rel$y[rel$from],
           rel$x[rel$to],rel$y[rel$to],col=linecol)
  if(show.points) points(rel$x,rel$y,...)
}
