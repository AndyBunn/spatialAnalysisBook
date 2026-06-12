## ----include=FALSE------------------------------------------------------------
library(tidyverse)
library(sf)
library(terra)
library(tidyterra)


## -----------------------------------------------------------------------------
californiaOzonePoints <- read_csv("../data/californiaOzonePoints.csv")
head(californiaOzonePoints)


## -----------------------------------------------------------------------------
californiaOzonePointsSF <- st_as_sf(californiaOzonePoints, coords = c("longitude", "latitude"), crs = 4326)
head(californiaOzonePointsSF)


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
ggplot() +
  geom_spatraster(data = californiaElevation) +
  scale_fill_terrain_c(name = "Elevation (m)") +
  geom_spatraster_contour(
    data = californiaElevation,
    color = "black", alpha = 0.4, bins = 10
  ) +
  geom_sf(
    data = californiaOzonePointsSF, pch = 21, color = "grey",
    fill = "white",
    size = 2, alpha = 0.5
  ) +
  theme_minimal()


## -----------------------------------------------------------------------------
elevationAtOzonePoints <- terra::extract(californiaElevation, californiaOzonePointsSF)
head(elevationAtOzonePoints)


## -----------------------------------------------------------------------------
californiaOzonePointsSF <- californiaOzonePointsSF %>%
  add_column(elevation = elevationAtOzonePoints$elev)
head(californiaOzonePointsSF)


## -----------------------------------------------------------------------------
californiaOzonePointsSF %>%
  ggplot(mapping = aes(x = elevation, y = ozone)) +
  geom_point() +
  geom_smooth() +
  labs(x = "Elevation (m)", y = "Ozone (ppb)") +
  theme_minimal()


## -----------------------------------------------------------------------------
californiaOzonePointsAlbers <- st_transform(californiaOzonePointsSF, crs = 3310)
californiaOzonePointsAlbers


## -----------------------------------------------------------------------------
californiaElevationAlbers <- project(californiaElevation, "EPSG:3310")
californiaElevationAlbers

