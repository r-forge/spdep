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


soi.graph <- function(tri.nb,coords){
  x <- coords
  if (!is.matrix(x)) stop("Data not in matrix form")
  if (any(is.na(x))) stop("Data cannot include NAs")
  np<-length(tri.nb)
  noedges<-0
  rad<-nearneigh<-rep(0,np)
  neigh<-unlist(tri.nb)  
  noneigh<-unlist(lapply(tri.nb,length))
  g1<-g2<-rep(0,sum(noneigh))
  answ<-.C("compute_soi", np=as.integer(np), from=as.integer(g1),
     to=as.integer(g2), nedges=as.integer(noedges),
     notri.nb=as.integer(noneigh), tri.nb=as.integer(neigh),
     nn=as.integer(nearneigh), 
     circles=as.double(rad), x=as.double(x[,1]), y=as.double(x[,2]),
     PACKAGE="spdep")
  answ$from<-answ$from[1:answ$nedges]
  answ$to<-answ$to[1:answ$nedges]
  answ<-list(np=answ$np,nedges=answ$nedges,
             from=answ$from,to=answ$to,circles=answ$circ)
  attr(answ, "call") <- match.call()
  class(answ)<-c("Graph","SOI")
  invisible(answ)
}
