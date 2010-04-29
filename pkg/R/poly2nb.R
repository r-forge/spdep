# Copyright 2001-2010 by Roger Bivand 
#
#
# Modified by Micah Altman 2010
	


poly2nb <- function(pl, row.names=NULL, snap=sqrt(.Machine$double.eps),
	queen=TRUE) {
	if (!inherits(pl, "polylist")) {
		if (extends(class(pl), "SpatialPolygons"))
			pl <- maptools:::.SpP2polylist(pl)
		else stop("Not a polygon list")
	}
	if (inherits(pl, "multiparts")) stop("Convert to newer polylist format")
	n <- length(pl)
	if (n < 1) stop("non-positive number of entities")
	if (is.null(row.names)) regid <- attr(pl, "region.id")
	else regid <- NULL
	if (is.null(regid)) {
		if(is.null(row.names)) regid <- as.character(1:n)
		else {
			if(length(row.names) != n)
				stop("row.names wrong length")
			else if (length(unique(row.names)) != length(row.names))
	    			stop("non-unique row.names given")
			else regid <- row.names
		}
	}
        if (snap < 0) snap <- abs(snap)
        if (snap < .Machine$double.eps) {
            bsnap <- .Machine$double.eps
        } else { 
            bsnap <- snap
        }

        genBBIndex<-function(pl,snap=bsnap) { 
            poly2bbs <- function(pl) {
                n <- length(pl)
                if (n < 1) stop("non-positive number of entities")
                res <- matrix(0, nrow = n, ncol = 4)
                for (i in 1:n) res[i, ] <- attr(pl[[i]], "bbox")
                res
            }
            bb<-poly2bbs(pl)
            if (storage.mode(bb) != "double") storage.mode(bb) <- "double"
            dsnap <- as.double(snap)
            bb[, 1] <- bb[, 1] - dsnap
            bb[, 2] <- bb[, 2] - dsnap
            bb[, 3] <- bb[, 3] + dsnap
            bb[, 4] <- bb[, 4] + dsnap

            bxv <- as.vector(bb[,c(1,3)])
            byv <- as.vector(bb[,c(2,4)])
            obxv <- order(bxv)
            rbxv <- c(1:(length(pl)*2))[obxv]
            mbxv <- match(1:(length(pl)*2),obxv)
            obyv <- order(byv)
            rbyv <- c(1:(length(pl)*2))[obyv]
            mbyv <- match(1:(length(pl)*2),obyv)

            return(list(bb=bb, bxv=bxv, byv=byv, obxv=obxv, obyv=obyv, 
                mbxv=mbxv, mbyv=mbyv, rbyv=rbyv, rbxv=rbxv))
        }
	dbsnap <- as.double(bsnap)
        dsnap <- as.double(snap)
	BBindex <- genBBIndex(pl)
	bb <- BBindex$bb

	nrs <- integer(n)
	for (i in 1:n) {
		pl[[i]] <- na.omit(pl[[i]][-1,])
		nrs[i] <- as.integer(nrow(pl[[i]]))
		pl[[i]] <- as.double(pl[[i]])
	}
	
#        findInBox <- function(i, sp, bigger=TRUE) {
#	    n <- dim(sp$bb)[1]
#	
#	# ! i1 > j3 --> i1 <= j3
#	    tmp1 <- sp$rbxv[sp$mbxv[i]:(n*2)] 
#	    tmp1 <- tmp1[which(tmp1>n)] - n
#	# ! i2 > j4 --> i2 <= bj4
#	    tmp2 <- sp$rbyv[sp$mbyv[i]:(n*2)] 
#	    tmp2 <- tmp2[which(tmp2>n)] - n
#	# ! i3 < j1 -> i3 >= j1
#	    tmp3 <- sp$rbxv[1:sp$mbxv[i+n]] 
#	    tmp3 <- tmp3[which(tmp3<=n)]
#	# ! i4 < j2 -> i4 >= j2
#	    tmp4 <- sp$rbyv[1:sp$mbyv[i+n]] 
#	    tmp4 <- tmp4[which(tmp4<=n)]
#	    result <- intersect(intersect(tmp1,tmp2), intersect(tmp3,tmp4))
#	    if (bigger) {
#		result <- result[which(result>i)]
#	    }
#	    return(sort(result))
#        }

# faster findInBox
       qintersect<-function(x,y) {
	    # streamlined intersect function for unique vectors
	    y[match(x, y, 0L)]
       }
        findInBox<-function(i,sp,bigger=TRUE) {
            n<-dim(sp$bb)[1]

# use index structure to identify which other BB's fall in i's BB
# by getting id's of polygons with BBmin_j < BBmax_i, BBmax_j > BBmin_i for x and y 
# then taking the intersection of these four lists of id's

	    tmp<-vector(mode="list", length=4)
        # ! i1 > j3 --> i1 <= j3
            tmp[[1]] <- sp$rbxv[sp$mbxv[i]:(n*2)]
            tmp[[1]]<- tmp[[1]][which(tmp[[1]]>n)] - n
        # ! i2 > j4 --> i2 <= bj4
            tmp[[2]] <- sp$rbyv[sp$mbyv[i]:(n*2)]
            tmp[[2]]<- tmp[[2]][which(tmp[[2]]>n)] - n
        # ! i3 < j1 -> i3 >= j1
            tmp[[3]] <- sp$rbxv[1:sp$mbxv[i+n]]
            tmp[[3]] <- tmp[[3]][which(tmp[[3]]<=n)]
        # ! i4 < j2 -> i4 >= j2
            tmp[[4]] <- sp$rbyv[1:sp$mbyv[i+n]]
            tmp[[4]]<- tmp[[4]][which(tmp[[4]]<=n)]

	# for performance, order the comparison of the lists

	    lentmp <- order(sapply(tmp,length))

	# use qintersect, since these are already vectors and unique 
	    result <- qintersect(tmp[[lentmp[2]]],tmp[[lentmp[1]]])
	    result <- qintersect(tmp[[lentmp[3]]],result)
	    result <- qintersect(tmp[[lentmp[4]]],result)

            if (bigger) {
                result<-result[which(result>i)]
            }
            return(sort(result))
        }

	polypoly2 <- function(poly1, nrs1, poly2, nrs2, snap) {
		if (any(nrs1 == 0 || nrs2 == 0)) return(as.integer(0))
		res <- .Call("polypoly", poly1, nrs1, poly2, 
			nrs2, snap, PACKAGE="spdep")
		res
	}

	ans <- vector(mode="list", length=n)
	for (i in 1:n) ans[[i]] <- integer(0)
	criterion <- ifelse(queen, 0, 1)
	for (i in 1:(n-1)) {
		#for (j in (i+1):n) {
		for (j in findInBox(i,BBindex)) {
			jhit <- .Call("spOverlap", bb[i,], 
				bb[j,], PACKAGE="spdep")
			if (jhit > 0) {
			    khit <- 0
			    khit <- polypoly2(pl[[i]], nrs[i], pl[[j]], 
				nrs[j], dsnap)

			    if (khit > criterion) {
				ans[[i]] <- c(ans[[i]], j)
				ans[[j]] <- c(ans[[j]], i)
			    }
			}
		}
	}
	for (i in 1:n) ans[[i]] <- sort(ans[[i]])
	class(ans) <- "nb"
	attr(ans, "region.id") <- regid
	attr(ans, "call") <- match.call()
	if (queen) attr(ans, "type") <- "queen"
	else attr(ans, "type") <- "rook"
	ans <- sym.attr.nb(ans)
	ans
}	


#poly2nb <- function(pl, row.names=NULL, snap=sqrt(.Machine$double.eps),
#	queen=TRUE) {
#	if (!inherits(pl, "polylist")) {
#		if (extends(class(pl), "SpatialPolygons"))
#			pl <- maptools:::.SpP2polylist(pl)
#		else stop("Not a polygon list")
#	}
#	if (inherits(pl, "multiparts")) stop("Convert to newer polylist format")
#	n <- length(pl)
#	if (n < 1) stop("non-positive number of entities")
#	if (is.null(row.names)) regid <- attr(pl, "region.id")
#	else regid <- NULL
#	if (is.null(regid)) {
#		if(is.null(row.names)) regid <- as.character(1:n)
#		else {
#			if(length(row.names) != n)
#				stop("row.names wrong length")
#			else if (length(unique(row.names)) != length(row.names))
#	    			stop("non-unique row.names given")
#			else regid <- row.names
#		}
#	}
#	poly2bbs <- function(pl) {
#		n <- length(pl)
#		if (n < 1) stop("non-positive number of entities")
#		res <- matrix(0, nrow=n, ncol=4)
#		for (i in 1:n) res[i,] <- attr(pl[[i]], "bbox")
#		res
#	}
#	bb <- poly2bbs(pl)
#	if (storage.mode(bb) != "double") storage.mode(bb) <- "double"
#	dsnap <- as.double(snap)
#	bb[,1] <- bb[,1] - dsnap
#	bb[,2] <- bb[,2] - dsnap
#	bb[,3] <- bb[,3] + dsnap
#	bb[,4] <- bb[,4] + dsnap
#	nrs <- integer(n)
#	for (i in 1:n) {
#		pl[[i]] <- na.omit(pl[[i]][-1,])
#		nrs[i] <- as.integer(nrow(pl[[i]]))
#		pl[[i]] <- as.double(pl[[i]])
#	}
#	
#
#
#	polypoly2 <- function(poly1, nrs1, poly2, nrs2, snap) {
#		if (any(nrs1 == 0 || nrs2 == 0)) return(as.integer(0))
#		res <- .Call("polypoly", poly1, nrs1, poly2, 
#			nrs2, snap, PACKAGE="spdep")
#		res
#	}
#
#	ans <- vector(mode="list", length=n)
#	for (i in 1:n) ans[[i]] <- integer(0)
#	criterion <- ifelse(queen, 0, 1)
#	for (i in 1:(n-1)) {
#		for (j in (i+1):n) {
#			jhit <- .Call("spInsiders", bb[i,], 
#				bb[j,], PACKAGE="spdep")
#			if (jhit > 0) {
#			    khit <- 0
#			    khit <- polypoly2(pl[[i]], nrs[i], pl[[j]], 
#				nrs[j],dsnap)
#
#			    if (khit > criterion) {
#				ans[[i]] <- c(ans[[i]], j)
#				ans[[j]] <- c(ans[[j]], i)
#			    }
#			}
#		}
#	}
#	for (i in 1:n) ans[[i]] <- sort(ans[[i]])
#	class(ans) <- "nb"
#	attr(ans, "region.id") <- regid
#	attr(ans, "call") <- match.call()
#	if (queen) attr(ans, "type") <- "queen"
#	else attr(ans, "type") <- "rook"
#	ans <- sym.attr.nb(ans)
#	ans
#}
#

