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


n.comp.nb <- function(nb.obj){
  if(!inherits(nb.obj,"nb"))stop("not a neighbours list")
  nb.obj <- make.sym.nb(nb.obj)
  comp <- rep(0,length(nb.obj))
  comp <- .Call("g_components", nb.obj, as.integer(comp), PACKAGE="spdep")
  answ <- list(nc=length(unique(comp)), comp.id=comp)
  answ
}

