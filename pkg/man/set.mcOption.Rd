\name{set.mcOption}
\alias{set.ClusterOption}
\alias{get.ClusterOption}
\alias{set.mcOption}
\alias{get.coresOption}
\alias{set.coresOption}
\alias{get.mcOption}
\alias{set.rlecuyerSeedOption}
\alias{get.rlecuyerSeedOption}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Options for parallel support}
\description{
  Provides support . 
}
\usage{
set.mcOption(value)
get.mcOption()
set.coresOption(value)
get.coresOption()
set.ClusterOption(cl)
get.ClusterOption()
set.rlecuyerSeedOption(seed)
get.rlecuyerSeedOption()
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{value}{valid replacement value}
  \item{cl}{a cluster object created by \code{makeCluster} in \pkg{parallel}}
  \item{seed}{an integer vector of length 6 to seed the Lecuyer streams for cluster RNG}
}
\details{
 .
}
\value{
  .
}

\author{Roger Bivand \email{Roger.Bivand@nhh.no}}
\note{
.
}

 
\examples{
library(parallel)
nc <- detectCores(logical=FALSE)
nc
coresOpt <- get.coresOption()
coresOpt
if (!is.na(nc)) {
 invisible(set.coresOption(nc))
 print(exists("aple.mc"))
 if(.Platform$OS.type == "windows") {
# forking not permitted on Windows - start cluster
  print(get.mcOption())
  cl <- makeCluster(get.coresOption())
  print(clusterEvalQ(cl, exists("aple.mc")))
  set.ClusterOption(cl)
  print(clusterEvalQ(cl, exists("aple.mc")))
  set.ClusterOption(NULL)
  print(clusterEvalQ(cl, exists("aple.mc")))
  stopCluster(cl)
 } else {
  mcOpt <- get.mcOption()
  print(mcOpt)
  print(mclapply(1:get.coresOption(), function(i) exists("aple.mc"),
   mc.cores=get.coresOption()))
  invisible(set.mcOption(FALSE))
  cl <- makeCluster(nc)
  print(clusterEvalQ(cl, exists("aple.mc")))
  set.ClusterOption(cl)
  print(clusterEvalQ(cl, exists("aple.mc")))
  set.ClusterOption(NULL)
  print(clusterEvalQ(cl, exists("aple.mc")))
  stopCluster(cl)
  invisible(set.mcOption(mcOpt))
 }
 invisible(set.coresOption(coresOpt))
}
}
\keyword{ spatial }
