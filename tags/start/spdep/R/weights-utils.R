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

write.nb.gal <- function(nb, file) {
	if(class(nb) != "nb") stop("not an nb object")
	n <- length(nb)
	con <- file(file, open="w")
	open(con, open="w")
	writeLines(paste(n), con)
	for (i in 1:n) {
		writeLines(paste(i, length(nb[[i]]), collapse=" "), con)
		writeLines(paste(nb[[i]], collapse=" "), con)
	}
	close(con)
}

is.symmetric.nb <- function(nb, verbose=TRUE, force=FALSE)
{
	if(class(nb) != "nb") stop("Not neighbours list")
	nbsym <- attr(nb, "sym")
	if(!is.null(nbsym)) res <- nbsym
	if(force || is.null(nbsym)) {
		res <- .Call("symtest", nb=nb, card=as.integer(card(nb)),
			verbose=as.logical(verbose))
	}
	if(!res && verbose) cat("Non-symmetric neighbours list\n")
	invisible(res)
}

sym.attr.nb <- function(nb) {
	if(class(nb) != "nb") stop("Not neighbours list")
	nbsym <- attr(nb, "sym")
	if(is.null(nbsym))
		attr(nb, "sym") <- is.symmetric.nb(nb, verbose=FALSE,
			force=TRUE)
	invisible(nb)
}

include.self <- function(nb) {
	if (!is.null(attributes(nb)$self.included) &&
		(as.logical(attributes(nb)$self.included)))
		stop("Self already included")
	n <- length(nb)
	for (i in 1:n) nb[[i]] <- sort(c(i, nb[[i]]))
	attr(nb, "self.included") <- TRUE
	invisible(nb)
}

# Copyright 2001 by Nicholas Lewin-Koh 

make.sym.nb <- function(nb){
	if(class(nb) != "nb") stop("Not neighbours list")
	if (is.symmetric.nb(nb, FALSE, TRUE)) {
		res <- nb
	} else {
        	k <- unlist(lapply(nb,length))
        	to <- unlist(nb)
        	from <- NULL
        	res <- vector(mode="list", length=length(nb))
        	for(i in 1:length(nb)){
        		from <- c(from,rep(i,k[i]))
        	}
        	for(i in 1:length(nb)){
        		res[[i]] <- sort(unique(c(to[from==i],from[to==i])))
        		if(length(res[[i]]) == 0) res[[i]] <- 0
        	}
        	attr(res, "region.id") <- attr(nb,"region.id")
        	attr(res, "call") <- attr(nb, "call")
        	attr(res, "type") <- attr(nb, "type")
        	attr(res, "sym") <- TRUE
        	class(res) <- "nb"
	}
	invisible(res)
}

