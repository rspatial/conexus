\name{helpers}

\alias{dst2svect}
\alias{dst2graph}

  
\title{Helper functions}

\description{
Helper functions to create a SpatVector of graph from the distances data.frame returned by \code{\link{XC}}
}

\usage{
dst2svect(x)
dst2graph(x)
}

\arguments{
\item{x}{data.frame "dst" returned by \code{\link{XC}}}
}

\value{SpatVector or graph}

\examples{
rnge <- terra::rast(system.file("ex/range.tif", package="conexus"))
land <- terra::readRDS(system.file("ex/land.rds", package="conexus"))
sp <- terra::readRDS(system.file("ex/sp.rds", package="conexus"))

x <- get_range(rnge, sp, land, include=50, exclude=250)
k <- get_samplesize(x, omega=1/40)
xy <- c(terra::init(x, "x"), terra::init(x, "y"))
set.seed(321)
zones <- make_zones(xy, k$range, k$n, spread=TRUE)

x <- XC(zones, sp[sp$sample_type=="seed"])
x$XC 
head(x$dist)

v <- dst2svect(x$dist)
g <- dst2graph(x$dist)
}

\keyword{spatial}

