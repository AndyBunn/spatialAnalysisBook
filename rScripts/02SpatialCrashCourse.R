## ----include=FALSE------------------------------------------------------------
library(tidyverse)
library(sf)
library(terra)


## -----------------------------------------------------------------------------
californiaOzonePoints <- read_csv("../data/californiaOzonePoints.csv")
head(californiaOzonePoints)


## -----------------------------------------------------------------------------
californiaOzonePointsSF <- st_as_sf(californiaOzonePoints, coords = c("longitude", "latitude"), crs = 4326)
californiaOzonePointsSF


## -----------------------------------------------------------------------------
plot(californiaOzonePointsSF)


## -----------------------------------------------------------------------------
californiaElevation <- rast("../data/californiaElevation.tif")
californiaElevation


## -----------------------------------------------------------------------------
plot(californiaElevation)


## -----------------------------------------------------------------------------
californiaElevationFt <- californiaElevation * 3.281
plot(californiaElevationFt)


## -----------------------------------------------------------------------------
elevDF <- as.data.frame(californiaElevation, xy = TRUE)
names(elevDF)[3] <- "elevation"

ggplot() +
  geom_raster(data = elevDF, aes(x = x, y = y, fill = elevation)) +
  geom_contour(data = elevDF, aes(x = x, y = y, z = elevation),
               color = "white", alpha = 0.4, bins = 15) +
  geom_sf(data = californiaOzonePointsSF, color = "red", size = 1.5) +
  scale_fill_viridis_c(name = "Elevation (m)") +
  coord_sf() +
  theme_minimal() +
  labs(x = NULL, y = NULL)


## -----------------------------------------------------------------------------
elevation_at_ozone_points <- terra::extract(californiaElevation, californiaOzonePointsSF)
head(elevation_at_ozone_points)


## -----------------------------------------------------------------------------
californiaOzonePointsAlbers <- st_transform(californiaOzonePointsSF, crs = 3310)
californiaOzonePointsAlbers


## -----------------------------------------------------------------------------
californiaElevationAlbers <- project(californiaElevation, "EPSG:3310")
californiaElevationAlbers

