nbcosts <- function(nb, data, method=c("euclidean", "maximum", "manhattan",
                                "canberra", "binary", "minkowski",
                                "mahalanobis"), p=2, cov, inverted=FALSE) {
#  if ((!require(parallel)) | (length(nb)<300))
#    clist <- lapply(1:length(nb), function(i)
#                    nbcost(data, i, nb[[i]], method,
#                           p, cov, inverted))
#  else {
#    if (.Platform$OS.type == "windows") {
#      cl <- makeCluster(getOption("cl.cores", 2))
#      clusterEvalQ(cl, library(spdep))
    if (any(card(nb) == 0L)) stop("nbcosts: no-neighbour nodes")
    nc <- n.comp.nb(nb)$nc
    if (nc > 1) stop("nbcosts:", nc, "disjoint connected subgraphs")
    if (missing(cov)) cov <- NULL
    cores <- get.coresOption()
    if (is.null(cores)) {
        parallel <- "no"
    } else {
        parallel <- ifelse (get.mcOption(), "multicore", "snow")
    }
    ncpus <- ifelse(is.null(cores), 1L, cores)
    cl <- NULL
    if (parallel == "snow") {
        cl <- get.ClusterOption()
        if (is.null(cl)) {
            parallel <- "no"
            warning("no cluster in ClusterOption, parallel set to no")
        }
    }
    if (length(nb)<300) parallel <- "no"
    if (parallel == "snow") {
        parallel <- "no"
        warning("no parallel calculations available")
    }
    
    if (parallel == "snow") {
#        require(parallel)
#        method <- match.arg(method)
#        data <- as.matrix(data)
#        sI <- splitIndices(length(nb), length(cl))
#        if (method=="mahalanobis") {
#            out <- parLapply(cl, X = sI, fun=function(I) {lapply(I,
#                FUN=function(i) {mahalanobis(data[nb[[i]],,drop=FALSE],
#                    data[i,,drop=FALSE], cov, inverted)})}, data=data,
#                     nb=nb, p=p, cov=cov, inverted=inverted)
#        } else {
#            out <- parLapply(cl, X = sI, fun=function(I) {lapply(I,
#                FUN=function(i) {id.neigh <- nb[[i]];
#                dt <- rbind(data[i,,drop=FALSE], data[id.neigh,,drop=FALSE]);
#                dist(dt, method=method, p=p)[1:length(id.neigh)]})},
#                data=data, nb=nb, method=method, p=p)
#        }
#        clist <- do.call("c", out)
    } else if (parallel == "multicore") {
        require(parallel)
        sI <- splitIndices(length(nb), ncpus)
        out <- mclapply(sI, FUN=lapply, function(i) {nbcost(data, i, nb[[i]],
            method, p, cov, inverted)}, mc.cores=ncpus)
        clist <- do.call("c", out)
    } else {
      clist <- lapply(1:length(nb),
                   function(i) nbcost(data, i, nb[[i]], method,
                           p, cov, inverted))
    }
    attr(clist, "call") <- match.call()
    attr(clist, "class") <- "nbdist"
    return(clist)
}

nbcost <- function(data, id, id.neigh,
                   method=c("euclidean", "maximum", "manhattan",
                     "canberra", "binary", "minkowski",
                     "mahalanobis"), p=2, cov, inverted=FALSE) {
  if (is.function(method))
    return(method(data, id, id.neigh))
  else {
    method <- match.arg(method)
    data <- as.matrix(data)
    if (method=="mahalanobis")
      return(mahalanobis(data[id.neigh,,drop=FALSE], data[id,,drop=FALSE],
cov, inverted))
    else
      return(dist(rbind(data[id,,drop=FALSE], data[id.neigh,,drop=FALSE]),
method=method,
                p=p)[1:length(id.neigh)])
  }
}


