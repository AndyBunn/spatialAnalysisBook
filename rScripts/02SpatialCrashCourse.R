## -----------------------------------------------------------------------------
#| label: packages
#| message: false
library(tidyverse)
library(sf)
library(terra)
library(tidyterra)


## -----------------------------------------------------------------------------
#| label: ca-ozone-load
californiaOzonePoints <- read_csv("../data/californiaOzonePoints.csv")
head(californiaOzonePoints)


## -----------------------------------------------------------------------------
#| label: ca-ozone-sf
californiaOzonePointsSF <- st_as_sf(californiaOzonePoints, coords = c("longitude", "latitude"), crs = 4326)
head(californiaOzonePointsSF)


## -----------------------------------------------------------------------------
#| label: ca-ozone-plot
plot(californiaOzonePointsSF)


## -----------------------------------------------------------------------------
#| label: ca-elevation-load
californiaElevation <- rast("../data/californiaElevation.tif")
californiaElevation


## -----------------------------------------------------------------------------
#| label: ca-elevation-plot
plot(californiaElevation)


## -----------------------------------------------------------------------------
#| label: ca-elevation-feet
californiaElevationFt <- californiaElevation * 3.281
plot(californiaElevationFt)


## -----------------------------------------------------------------------------
#| label: ca-elevation-map
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
#| label: extract-elevation
elevationAtOzonePoints <- terra::extract(californiaElevation, californiaOzonePointsSF)
head(elevationAtOzonePoints)


## -----------------------------------------------------------------------------
#| label: attach-elevation
californiaOzonePointsSF <- californiaOzonePointsSF %>%
  add_column(elevation = elevationAtOzonePoints$elev)
head(californiaOzonePointsSF)


## -----------------------------------------------------------------------------
#| label: elevation-ozone-plot
californiaOzonePointsSF %>%
  ggplot(mapping = aes(x = elevation, y = ozone)) +
  geom_point() +
  geom_smooth() +
  labs(x = "Elevation (m)", y = "Ozone (ppb)") +
  theme_minimal()


## -----------------------------------------------------------------------------
#| label: ca-ozone-albers
californiaOzonePointsAlbers <- st_transform(californiaOzonePointsSF, crs = 3310)
californiaOzonePointsAlbers


## -----------------------------------------------------------------------------
#| label: ca-elevation-albers
californiaElevationAlbers <- project(californiaElevation, "EPSG:3310")
californiaElevationAlbers

