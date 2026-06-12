## ----echo=FALSE, message=FALSE, warning=FALSE, results='hide'-----------------
set.seed(733)


## ----message=FALSE, warning=FALSE---------------------------------------------
library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(gstat)
library(cowplot)


## -----------------------------------------------------------------------------
meuse2 <- readRDS("../data/meuse2.Rds")
meuse.grid2 <- readRDS("../data/meuse.grid2.Rds")


## -----------------------------------------------------------------------------
meusePointsSf <- st_as_sf(meuse2, coords = c("x", "y"), crs = 28992)
meuseGridSf <- st_as_sf(meuse.grid2, coords = c("x", "y"), crs = 28992)

covarsSf <- meuseGridSf %>% select(river_dist_m)
meusePointsSf <- st_join(meusePointsSf, covarsSf, join = st_nearest_feature)


## -----------------------------------------------------------------------------
meusePointsSf <- meusePointsSf %>%
  mutate(log_dist = log(river_dist_m)) %>%
  drop_na(om)
meuseGridSf <- meuseGridSf %>%
  mutate(log_dist = log(river_dist_m))


## -----------------------------------------------------------------------------
omLm <- lm(om ~ log_dist, data = meusePointsSf)
summary(omLm)


## -----------------------------------------------------------------------------
meuseGridSf$omHat <- predict(omLm, newdata = meuseGridSf)

# same sf2Rast function from the IDW notes -- write it once, use it everywhere
sf2Rast <- function(sfObject, variable2get = 1) {
  tmp <- sfObject[, variable2get] %>% st_drop_geometry()
  dfObject <- data.frame(st_coordinates(sfObject), z = tmp)
  rastObject <- rast(dfObject, crs = crs(sfObject))
  return(rastObject)
}

meuseGridRast <- sf2Rast(meuseGridSf, variable2get = "omHat")
ggplot() +
  geom_spatraster(
    data = meuseGridRast,
    mapping = aes(fill = omHat)
  ) +
  scale_fill_terrain_c() +
  labs(fill = "omHat")


## ----echo=FALSE---------------------------------------------------------------
w <- spdep::knn2nb(spdep::knearneigh(meusePointsSf, k = 4))
foo <- spdep::moran.test(residuals(omLm), spdep::nb2listw(w))


## -----------------------------------------------------------------------------
omVar <- variogram(om ~ 1, meusePointsSf)
plot(omVar,
  pch = 20, cex = 1.5, col = "black",
  ylab = expression("Semivariance (" * gamma * ")"),
  xlab = "Distance (m)", main = "% Soil Organic Matter"
)


## -----------------------------------------------------------------------------
omGstat <- gstat(
  id = "omModel", formula = om ~ log_dist,
  data = meusePointsSf
)
omGstatObsVariogram <- variogram(omGstat)
plot(omGstatObsVariogram,
  pch = 20, cex = 1.5, col = "black",
  ylab = expression("Semivariance (" * gamma * ")"),
  xlab = "Distance (m)", main = "Model Residuals"
)


## -----------------------------------------------------------------------------
omGauVariogramModel <- vgm(psill = 2, model = "Gau", range = 500, nugget = 4)
omGauFittedVariogram <- fit.variogram(
  object = omGstatObsVariogram,
  model = omGauVariogramModel
)
plot(omGstatObsVariogram, omGauFittedVariogram,
  pch = 20, cex = 1.5, col = "black",
  ylab = expression("Semivariance (" * gamma * ")"),
  xlab = "Distance (m)", main = "Model Residuals"
)


## -----------------------------------------------------------------------------
# Update the gstat object with the variogram:
omGstatWVariogram <- gstat(omGstat, id = "omModel", model = omGauFittedVariogram)
# And predict
omHatSf <- predict(omGstatWVariogram, newdata = meuseGridSf)
omHatSf


## -----------------------------------------------------------------------------
#| fig-width: 10
#| fig-height: 9
omHatRast <- sf2Rast(omHatSf, variable2get = "omModel.pred")

# kriged residual surface = RK prediction minus OLS trend
omHatSf$krigedResid <- omHatSf$omModel.pred - meuseGridSf$omHat
omResidRast <- sf2Rast(omHatSf, variable2get = "krigedResid")

pOls <- ggplot() +
  geom_spatraster(data = meuseGridRast, mapping = aes(fill = omHat)) +
  scale_fill_terrain_c(limits = c(0, 18)) +
  labs(title = "OLS trend", fill = "om (%)") +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

pResid <- ggplot() +
  geom_spatraster(data = omResidRast, mapping = aes(fill = krigedResid)) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0, na.value = "transparent"
  ) +
  labs(title = "Kriged residuals", fill = "resid") +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

pRk <- ggplot() +
  geom_spatraster(data = omHatRast, mapping = aes(fill = omModel.pred)) +
  scale_fill_terrain_c(limits = c(0, 18)) +
  labs(fill = "om (%)") +
  theme(axis.text = element_blank(), axis.ticks = element_blank())


eq <- ggdraw() + draw_label("OLS trend  +  Kriged residuals  =  RK prediction",
  fontface = "italic", size = 11
)

topRow <- plot_grid(pOls, pResid, ncol = 2)
plot_grid(topRow, eq, pRk, nrow = 3, rel_heights = c(1, 0.1, 1))


## -----------------------------------------------------------------------------
omHatVarRast <- sf2Rast(omHatSf, variable2get = "omModel.var")
omHatVarRast <- omHatVarRast %>% mutate(omModel.var.sqrt = sqrt(omModel.var))
ggplot() +
  geom_spatraster(
    data = omHatVarRast,
    mapping = aes(fill = omModel.var.sqrt)
  ) +
  scale_fill_terrain_c() +
  labs(fill = "Organic Matter SD (%)")

