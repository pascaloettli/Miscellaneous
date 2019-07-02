# Loading required packages
libs <- c("raster", "xts", "lubridate", "foreach", "doMC", "Deducer", "rasterVis", "latticeExtra", "rgeos", "rgdal")
lapply(libs, library, character.only = TRUE)

# Period of work
startYear <- 1983
endYear <- 2014
Period <- paste(startYear,'-01-01/',endYear,'-12-31',sep='')

# Season of work (easier to keep this one and adjust...
Seas <- 1:3

# ...the time delay here)
mths <- 0


# Selection of the area of interest
e <- extent(33, 293, -33, 33)

# Reading SSTA from NetCDF file
# Directory and file name need to be adjusted
anom.sst <- brick('sst.mnanom_clm8110.v4.nc')

# Dates from the SSTA NetCDF file
Dates.anom.sst <- as.Date(getZ(anom.sst))
Dates.anom.sst <- Dates.anom.sst - day(Dates.anom.sst) + days(1)

# Period selection
dstart <- paste(startYear,'01','01',sep='-'); dstart <- grep(dstart, Dates.anom.sst)
dend <- paste(endYear,'12','01',sep='-'); dend <- grep(dend, Dates.anom.sst)
Dates.anom.sst <- Dates.anom.sst[dstart:dend]

# Data selection based on period
anom.sst <- subset(anom.sst, dstart:dend)

# Data selection based on area
anom.sst <- crop(anom.sst, e)

# Extract data from Raster* and convert to time-serie
anom.sst.ts <- getValues(anom.sst)
anom.sst.ts <- xts(t(anom.sst.ts), order.by = Dates.anom.sst)
anom.sst.ts <- anom.sst.ts[Period]

# YR0 is the selection of years for the composite
# YR1 is the remaining years
YR0 <- c(1991,1992,1998,2002,2010)
YR1 <- (1987:2014); YR1 <- YR1[!YR1 %in% YR0]

# Selection of correspondig seasons
# Better to adjust the selection with mths, i.e. -1 selects Dec, Jan and Feb, +1 selects Feb, Mar and Apr.
seas0 <- as.Date(paste(rep(YR0, each = length(Seas)), rep(Seas, length(YR0)), 1, sep = "-")) %m+% months(mths)
seas1 <- as.Date(paste(rep(YR1, each = length(Seas)), rep(Seas, length(YR1)), 1, sep = "-")) %m+% months(mths)

# Seasonal average
s.pos <- anom.sst.ts[seas0]
ep <- c(0, seq(length(Seas), NROW(s.pos), length(Seas)))
s.pos <- period.apply(s.pos, ep, mean)
s.pos <- t(as.matrix(s.pos))

s.neut <- anom.sst.ts[seas1]
ep <- c(0, seq(length(Seas), NROW(s.neut), length(Seas)))
s.neut <- period.apply(s.neut, ep, mean)
s.neut <- t(as.matrix(s.neut))

# Number of replication for the t-test
B <- 500

# Resgister parallel backend
dc <- detectCores()
registerDoMC(cores = dc)

# t-test
out <-
  foreach(pg = seq_len(NROW(s.pos)), .export = "perm.t.test") %dopar% {
    obj <- try(perm.t.test(s.pos[pg,], s.neut[pg,], statistic = "t", alternative = "two.sided", B = B), silent=TRUE)
    if (is(obj, "try-error")) return(NA) else return(obj)
  }

# Get p-value
idx <- seq_along(out)
out.pval <- rep(NA, length(idx))
idx <- idx[!is.na(out)]
for(ii in idx){
  out.pval[ii] <- out[[ii]]$p.value
}

# Significance
sig <- 0.05

# Mean and p-value to raster 
pos.st <- pos.st.pvalue <- init(anom.sst, v='cell')
pos.st <- setValues(pos.st, rowMeans(s.pos))
pos.st.pvalue <- setValues(pos.st.pvalue, ifelse(out.pval <= sig, 1, NA))

# Land 
# From https://www.naturalearthdata.com/downloads/110m-physical-vectors/110m-land/
# Directory need to be adjusted
msk <- readOGR("NaturalEarthData","ne_110m_land")
CP1 <- as(extent(xmin(e),180,ymin(e),ymax(e)), "SpatialPolygons")
proj4string(CP1) <- CRS(proj4string(msk))
CP2 <- as(extent(-180,xmax(e)-360+1,ymin(e),ymax(e)), "SpatialPolygons")
proj4string(CP2) <- CRS(proj4string(msk))
ne_10m_land1 <- gIntersection(msk, CP1, byid=TRUE)
ne_10m_land2 <- gIntersection(msk, CP2, byid=TRUE)
ne_10m_land2 <- shift(ne_10m_land2, x = 360)

# Map options
x.scale <- list(fontfamily='FreeSans',
  font=12,
  fontface=2,
  cex=1.2, 
  alternating=1)
y.scale <- list(fontfamily='FreeSans',
  font=12,
  fontface=2,
  cex=1.2, 
  alternating=1)

# Map
levelplot(pos.st.pvalue, at = c(0,1), col.regions = "grey", maxpixels = 1e8, 
  colorkey = FALSE, margin = FALSE, main = "", panel = "panel.levelplot",
  scales=list(x=x.scale, y=y.scale), xlab = "", ylab = "") +
  contourplot(pos.st, at = seq(-5,-0.2,0.2), labels = list(col = "blue"), label.style = "align", lty = 2, lwd = 1) +
  contourplot(pos.st, at = seq(0.2,5,0.2), labels = list(col = "blue"), label.style = "align", lty = 1, lwd = 1) +
  layer(sp.polygons(ne_10m_land1, col="black", fill="white")) +
  layer(sp.polygons(ne_10m_land2, col="black", fill="white")) 
