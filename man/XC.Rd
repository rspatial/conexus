\name{XC}

\alias{XC}

\title{Compute the XC score}

\description{
The ex-situ conservation score expresses the degree to which a taxon is conserved ex-situ. XC ranges from zero (no conservation) to 1 (highly conserved).

The XC needs the taxon's range that is divided into a number of zones, and optionally a function to combine the geographic and environmental distance between zones into a single distance metric (otherwise only geographic distance is used). The distance metric serves as weights on the vertices of a network connecting neighboring zones. 

The ex-situ conservation score is the sum of the length of all the remaining shortest paths from or to the vertices in the network, after setting the weights between zones from which seed samples are conserved to zero, and halving the weights between neighboring zones of which only one has a seed sample, divided by the sum of all shortest path lengths. 

The score can be adjusted for small sample sizes and for samples with unknown geographic origin.
}

\usage{
XC(regions, sample, env=NULL, envfun=NULL, minssize=10, maxdist=1500)
}

\arguments{
\item{regions}{SpatVector or polygons}
\item{sample}{SpatVector or points}
\item{env}{SpatRaster with environmental data or NULL}
\item{envfun}{function that transforms environmental distance to geographic distance}
\item{minssize}{postive integer. The minimal sample size needed. Below that threshold, a penalty is applied}
\item{maxdist}{postive number. The maximum distance between zones. Larger distances are clamped to this number}
}

\value{
list with a numeric value \code{XC}, the ex-situ conservation score; and a data.frame with the a row for all adjacent zones with the following variables:

 \itemize{
    \code{from} "from" zone ID 
	\code{to}  "to" zone ID
	\code{adj} Are the zones adjacent?
	\code{xf} "from" x coordinate
	\code{yf} "from" y coordinate
	\code{xt} "to" x coordinate
	\code{yt} "to" y coordinate
	\code{w} relative vertex weight. 0 (both zones have a sample), 0.5 (one zone has a sample) and 1 (neither zone has a sample)
	\code{dst} total distance between vertices
	\code{envdst} environmental distance between vertices
	\code{geodst} geographic distance between vertices
	\code{weight} edge weight, computd as \code{w * dst} 
}
 
This data.frame can be transformed to a SpatVector with \code{\link{dst2svect}} and to a graph with \code{\link{dst2graph}}
}

\examples{
rnge <- terra::rast(system.file("ex/range.tif", package="conexus"))
land <- terra::readRDS(system.file("ex/land.rds", package="conexus"))
sp <- terra::readRDS(system.file("ex/sp.rds", package="conexus"))
mprc <- readRDS(system.file("ex/m_prc.rds", package="conexus"))
mtmp <- readRDS(system.file("ex/m_tmp.rds", package="conexus"))
env <- terra::rast(system.file("ex/env.tif", package="conexus"))

x <- get_range(rnge, sp, land, include=50, exclude=250)
k <- get_samplesize(x, omega=1/40)
xy <- c(terra::init(x, "x"), terra::init(x, "y"))

set.seed(321)
zones <- make_zones(xy, k$range, k$n, spread=TRUE)

envdist <- function(x) {
	x$tmp[x$tmp > 13] <- 13 
	x$prc[x$prc > 2000] <- 2000
	p_tmp <- predict(mtmp, x)
	p_pr <- predict(mprc, x)
	rowMeans(cbind(p_tmp, p_pr))
}

x <- XC(zones, sp[sp$sample_type=="seed"], env=NULL, envfun=NULL)
x$XC 
head(x$dist)

y <- XC(zones, sp[sp$sample_type=="seed"], env=env, envfun=envdist)
y$XC 
head(y$dist)
}

\keyword{spatial}

