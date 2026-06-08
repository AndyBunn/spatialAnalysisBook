# buildBirdRichnessMexico.R
# Adds WorldClim bioclimatic variables to the Mexico bird richness dataset.
# Source: WorldClim v2.1, 10 arc-minute resolution (Fick & Hijmans 2017)
# Variables extracted:
#   mat          BIO1  Annual Mean Temperature (degrees C * 10)
#   tempSeason   BIO4  Temperature Seasonality (SD * 100)
#   tempRange    BIO7  Temperature Annual Range (degrees C * 10)
#   map          BIO12 Annual Precipitation (mm)
#   precipSeason BIO15 Precipitation Seasonality (CV)
#   precipDryQ   BIO17 Precipitation of Driest Quarter (mm)
#
# Output: data/birdRichnessMexico.rds (overwrites original)

library(sf)
library(terra)
library(geodata)

birds <- readRDS("data/birdRichnessMexico.rds")

# Download WorldClim bioclim at 10 arc-min resolution
wc <- worldclim_global(var = "bio", res = 10, path = tempdir())

wc_sub <- wc[[c(1, 4, 7, 12, 15, 17)]]
names(wc_sub) <- c("mat", "tempSeason", "tempRange", "map", "precipSeason", "precipDryQ")

# Extract to bird locations (transform to WGS84 for extraction)
birds_wgs <- st_transform(birds, 4326)
pts <- vect(birds_wgs)
clim <- extract(wc_sub, pts, ID = FALSE)

# Sanity check
stopifnot(nrow(clim) == nrow(birds))
stopifnot(all(colSums(is.na(clim)) == 0))

# Bind back to original sf object (keeping original projection)
birds <- cbind(birds, clim)

saveRDS(birds, "data/birdRichnessMexico.rds")
cat("Saved enriched birdRichnessMexico.rds with", ncol(birds) - 1, "variables +  geometry\n")
cat("Variables:", paste(names(birds)[names(birds) != "geometry"], collapse = ", "), "\n")
