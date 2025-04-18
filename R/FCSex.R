#In GapAnalysis::GRSex the CA50 areas outside of the SDM range are considered for the genebank area and in GapAnalysis::ERSex the number of ecoregions considered for the genebank samples includes areas that are outside the SDM. 
# In both cases you could end up with a score > 1. That can be corrected for, but this is inconsistent. Areas outside the range (SDM) should not be included in the analysis.
# They are not included in the analysis below (but can be with "inrange=FALSE").


GRex <- function(srange, ca, inrange) {
	if (inrange) {
		bufr <- terra::mask(srange, ca)
		e <- terra::expanse(c(srange, bufr), unit="km")
		e[2,2] / e[1,2]
	} else {
		bufr <- terra::expanse(ca, unit="km")
		e <- terra::expanse(srange, unit="km")
		out <- bufr[1,2] / e[1,2]
		min(1, out)
	}
}

ERex <- function(srange, ca, ecoregions, inrange) {
	eco <- terra::mask(ecoregions, srange) 
	ueco <- nrow(terra::unique(eco))
	if (inrange) {
		seedeco <- terra::mask(eco, ca)
	} else {
		seedeco <- terra::mask(ecoregions, ca)	
	}
	seco <- nrow(terra::unique(seedeco))
	if (is.null(seco)) return(0)
	out <- seco / ueco
	if (inrange) {
		out 
	} else {
		min(1, out)
	}
}

SRex <- function(s, h) {
	if (s == 0) return(0)
	if (h == 0) return(1)
	min(1, s / h)
}

FCex <- function(seed, nherbarium, srange, ecoregions, bsize=50, nseed_nogeo=0, inrange=TRUE) {
	srange <- terra::subst(srange, 0, NA)
	ca <- terra::buffer(seed, bsize*1000) |> terra::aggregate()  
	ca <- terra::rasterize(ca, srange)
	g <- GRex(srange, ca, inrange)
	e <- ERex(srange, ca, ecoregions, inrange)
	nseed <- NROW(seed) + nseed_nogeo[1]
	s <- SRex(nseed, nherbarium)
	c(nseed=NROW(seed), nherbarium=NROW(nherbarium), nseed_nogeo=nseed_nogeo, GRex=g, ERex=e, SRex=s, FCex=mean(c(g, e, s)))
}

