\name{FCS}

\alias{GRex}
\alias{ERex}
\alias{SRex}
\alias{FCex}

\title{FCS}

\description{
Implementation of the FCS ex-situ method following 

Ramírez-Villegas, J., Khoury, C., Jarvis, A., Debouck, D. G., & Guarino, L. (2010). A gap analysis methodology for collecting crop genepools: a case study with Phaseolus beans. PloS one, 5(10), e13497.

Carver, D., Sosa, C. C., Khoury, C. K., Achicanoy, H. A., Diaz, M. V., Sotelo, S., Castañeda-Álvarez, N. P., & Ramirez-Villegas, J. (2021). GapAnalysis: an R package to calculate conservation indicators using spatial information. Ecography, 44(7), 1000–1009. doi.org/10.1111/ecog.05430
}

\usage{
GRex(srange, ca, inrange)
ERex(srange, ca, ecoregions, inrange)
SRex(s, h)
FCex(seed, nherbarium, srange, ecoregions, bsize=50, nseed_nogeo=0, inrange=TRUE)
}


\arguments{
\item{srange}{SpatRaster with the Species range (\code{1} for preseence, \code{NA} in other cells)}
\item{ca}{SpatVector with circular area (buffer) around each seed observation}
\item{ecoregions}{SpatRaster of ecoregions}
\item{inrange}{logical. If \code{TRUE} the area of occurrence is restricted to the area of srange}
\item{s}{number of seed samples}
\item{h}{number of non-seed (e.g. herbarium) samples}
\item{seed}{SpatVector of points with seed sample locations}
\item{nherbarium}{non negative integer. The number of herbarium (and other non-seed) records}
\item{bsize}{Positive number, buffer size in km}
\item{nseed_nogeo}{non-negative integer. The number of seed records without georeference}
}

\value{
numeric 
}

\examples{

srnge <- terra::rast(system.file("ex/range.tif", package="conexus"))
#this is a bad example for ecoregions
ereg <- terra::readRDS(system.file("ex/land.rds", package="conexus"))
sp <- terra::readRDS(system.file("ex/sp.rds", package="conexus"))

ecoreg <- terra::rasterize(ereg, srnge, 1:nrow(ereg))

seed <- sp[sp$sample_type == "seed", ]
herb <- nrow(sp[sp$sample_type != "seed", ])

FCex(seed, herb, srnge, ecoreg, bsize=50, inrange=FALSE)

# the components 
ca <- terra::buffer(seed, 50000)
ca <- terra::rasterize(ca, srnge)
GRex(srnge, ca, FALSE)
ERex(srnge, ca, ecoreg, FALSE)
SRex(nrow(seed), nrow(herb))
}


\keyword{spatial}

