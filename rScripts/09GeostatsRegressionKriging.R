## ----echo=FALSE, message=FALSE, warning=FALSE, results='hide'-----------------
#| label: setup
set.seed(733)


## ----message=FALSE, warning=FALSE---------------------------------------------
#| label: packages
library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(gstat)
library(cowplot)


## -----------------------------------------------------------------------------
#| label: meuse-data
meuse2 <- readRDS("../data/meuse2.Rds")
meuse.grid2 <- readRDS("../data/meuse.grid2.Rds")


## -----------------------------------------------------------------------------
#| label: extract-covariates
meusePointsSf <- st_as_sf(meuse2, coords = c("x", "y"), crs = 28992)
meuseGridSf <- st_as_sf(meuse.grid2, coords = c("x", "y"), crs = 28992)

covarsSf <- meuseGridSf %>% select(river_dist_m)
meusePointsSf <- st_join(meusePointsSf, covarsSf, join = st_nearest_feature)


## -----------------------------------------------------------------------------
#| label: transforms
meusePointsSf <- meusePointsSf %>%
  mutate(logLead = log(lead), logDist = log(river_dist_m)) %>%
  drop_na(logLead)
meuseGridSf <- meuseGridSf %>%
  mutate(logDist = log(river_dist_m))


## -----------------------------------------------------------------------------
#| label: ols-fit
leadLm <- lm(logLead ~ logDist, data = meusePointsSf)
summary(leadLm)


## -----------------------------------------------------------------------------
#| label: ols-surface
meuseGridSf$leadHat <- predict(leadLm, newdata = meuseGridSf)

# same sf2Rast function from the IDW notes -- write it once, use it everywhere
sf2Rast <- function(sfObject, variable2get = 1) {
  tmp <- sfObject[, variable2get] %>% st_drop_geometry()
  dfObject <- data.frame(st_coordinates(sfObject), z = tmp)
  rastObject <- rast(dfObject, crs = crs(sfObject))
  return(rastObject)
}

meuseGridRast <- sf2Rast(meuseGridSf, variable2get = "leadHat")
ggplot() +
  geom_spatraster(
    data = meuseGridRast,
    mapping = aes(fill = leadHat)
  ) +
  scale_fill_terrain_c() +
  labs(fill = "leadHat")


## ----echo=FALSE---------------------------------------------------------------
#| label: ols-moran
w <- spdep::knn2nb(spdep::knearneigh(meusePointsSf, k = 4))
foo <- spdep::moran.test(residuals(leadLm), spdep::nb2listw(w))


## -----------------------------------------------------------------------------
#| label: lead-variogram
leadVar <- variogram(logLead ~ 1, meusePointsSf)
plot(leadVar,
  pch = 20, cex = 1.5, col = "black",
  ylab = expression("Semivariance (" * gamma * ")"),
  xlab = "Distance (m)", main = "Log Lead"
)


## -----------------------------------------------------------------------------
#| label: resid-variogram
leadGstat <- gstat(
  id = "leadModel", formula = logLead ~ logDist,
  data = meusePointsSf
)
leadGstatObsVariogram <- variogram(leadGstat)
plot(leadGstatObsVariogram,
  pch = 20, cex = 1.5, col = "black",
  ylab = expression("Semivariance (" * gamma * ")"),
  xlab = "Distance (m)", main = "Model Residuals"
)


## -----------------------------------------------------------------------------
#| label: resid-variogram-fit
leadGauVariogramModel <- vgm(psill = 0.15, model = "Gau", range = 500, nugget = 0.15)
leadGauFittedVariogram <- fit.variogram(
  object = leadGstatObsVariogram,
  model = leadGauVariogramModel
)
plot(leadGstatObsVariogram, leadGauFittedVariogram,
  pch = 20, cex = 1.5, col = "black",
  ylab = expression("Semivariance (" * gamma * ")"),
  xlab = "Distance (m)", main = "Model Residuals"
)


## -----------------------------------------------------------------------------
#| label: rk-predict
# Add the fitted variogram to the existing gstat object:
leadGstatWVariogram <- gstat(leadGstat, id = "leadModel", model = leadGauFittedVariogram)
# And predict
leadHatSf <- predict(leadGstatWVariogram, newdata = meuseGridSf)
leadHatSf


## -----------------------------------------------------------------------------
#| label: rk-surfaces
#| fig-width: 10
#| fig-height: 9
leadHatRast <- sf2Rast(leadHatSf, variable2get = "leadModel.pred")

# kriged residual surface = RK prediction minus OLS trend
leadHatSf$krigedResid <- leadHatSf$leadModel.pred - meuseGridSf$leadHat
leadResidRast <- sf2Rast(leadHatSf, variable2get = "krigedResid")

pOls <- ggplot() +
  geom_spatraster(data = meuseGridRast, mapping = aes(fill = leadHat)) +
  scale_fill_terrain_c(limits = c(3.5, 7.5)) +
  labs(title = "OLS trend", fill = "log lead") +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

pResid <- ggplot() +
  geom_spatraster(data = leadResidRast, mapping = aes(fill = krigedResid)) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0, na.value = "transparent"
  ) +
  labs(title = "Kriged residuals", fill = "resid") +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

pRk <- ggplot() +
  geom_spatraster(data = leadHatRast, mapping = aes(fill = leadModel.pred)) +
  scale_fill_terrain_c(limits = c(3.5, 7.5)) +
  labs(fill = "log lead") +
  theme(axis.text = element_blank(), axis.ticks = element_blank())


eq <- ggdraw() + draw_label("OLS trend  +  Kriged residuals  =  RK prediction",
  fontface = "italic", size = 11
)

topRow <- plot_grid(pOls, pResid, ncol = 2)
plot_grid(topRow, eq, pRk, nrow = 3, rel_heights = c(1, 0.1, 1))


## -----------------------------------------------------------------------------
#| label: rk-variance
leadHatVarRast <- sf2Rast(leadHatSf, variable2get = "leadModel.var")
leadHatVarRast <- leadHatVarRast %>% mutate(leadModel.var.sqrt = sqrt(leadModel.var))
ggplot() +
  geom_spatraster(
    data = leadHatVarRast,
    mapping = aes(fill = leadModel.var.sqrt)
  ) +
  scale_fill_terrain_c() +
  labs(fill = "Log Lead SD")

