## ----echo=FALSE, message=FALSE, warning=FALSE, results='hide'-----------------
set.seed(184)


## ----message=FALSE------------------------------------------------------------
library(sf)
library(gstat)
library(tidyverse)
library(terra)
library(tidyterra)
library(automap)


## -----------------------------------------------------------------------------
# load
meuse2 <- readRDS("../data/meuse2.Rds")
meuse.grid2 <- readRDS("../data/meuse.grid2.Rds")

# make a variable to work with
meuse2$logLead <- log(meuse2$lead)
# make into sf
meuseSf <- st_as_sf(meuse2, coords = c("x", "y")) %>%
  st_set_crs(value = 28992)

meuseGridSf <- st_as_sf(meuse.grid2,
  coords = c("x", "y"),
  crs = st_crs(meuseSf)
)
meuseGridSf

p1 <- ggplot(data = meuseSf) +
  geom_sf(aes(fill = logLead),
    size = 4,
    shape = 21, color = "white", alpha = 0.8
  ) +
  scale_fill_continuous(type = "viridis", name = "log(ppm)") +
  labs(title = "Lead concentrations")
p1


## -----------------------------------------------------------------------------
leadVar <- variogram(logLead ~ 1, meuseSf)
plot(leadVar,
  pch = 20, cex = 1.5, col = "black",
  ylab = expression("Semivariance (" * gamma * ")"),
  xlab = "Distance (m)", main = "Lead concentrations (log(ppm))"
)


## -----------------------------------------------------------------------------
# note our initial estimates for the partial sill, range, and nugget.
sph.model <- vgm(psill = 0.5, model = "Sph", range = 750, nugget = 0.05)
sph.fit <- fit.variogram(object = leadVar, model = sph.model)
sph.fit # look at the fitted values
plot(leadVar,
  model = sph.fit, pch = 20, cex = 1.5, col = "black",
  ylab = expression("Semivariance (" * gamma * ")"),
  xlab = "Distance (m)", main = "Lead concentrations (log(ppm))",
  sub = "Points: Empirical, Line: Spherical Model"
)


## -----------------------------------------------------------------------------
# note our initial estimates for the sill, range, and nugget
exp.model <- vgm(psill = 0.5, model = "Exp", range = 750, nugget = 0.05)
exp.fit <- fit.variogram(object = leadVar, model = exp.model)
plot(leadVar,
  model = exp.fit, pch = 20, cex = 1.5, col = "black",
  ylab = expression("Semivariance (" * gamma * ")"),
  xlab = "Distance (m)", main = "Lead concentrations (log(ppm))",
  sub = "Points: Empirical, Line: Exponential Model"
)


## -----------------------------------------------------------------------------
foo <- data.frame(
  x = c(1, 3, 1, 4, 5),
  y = c(5, 4, 3, 5, 1),
  z = c(100, 105, 105, 100, 115)
)
foo
p2 <- ggplot() +
  geom_point(data = foo, aes(x = x, y = y, size = z)) +
  lims(x = c(0, 6), y = c(0, 6))
p2


## -----------------------------------------------------------------------------
p2 <- p2 +
  geom_point(aes(x = 2, y = 4), color = "red", size = 10, shape = 0) +
  geom_point(aes(x = 2, y = 4), color = "red", size = 6, shape = 63)
p2


## ----echo=FALSE---------------------------------------------------------------
d2s0 <- as.matrix(dist(cbind(c(foo$x, 2), c(foo$y, 4))))[1:5, 6]


## ----echo=FALSE---------------------------------------------------------------
sphGamma <- function(d, psill = 10, range = 6) {
  psill * (1.5 * (d / range) - 0.5 * (d / range)^3)
}
g <- sphGamma(d2s0)


## ----echo=FALSE---------------------------------------------------------------
dMat <- as.matrix(dist(cbind(foo$x, foo$y)))
G <- sphGamma(dMat)
diag(G) <- 0


## ----echo=FALSE---------------------------------------------------------------
lambda <- solve(G) %*% g
zhat <- sum(lambda[, 1] * foo$z)


## ----echo=FALSE---------------------------------------------------------------
0.2626 * (100) + 0.4985 * (105) + 0.2652 * (105) - 0.0165 * (100) - 0.0353 * (115)


## -----------------------------------------------------------------------------
leadVar <- variogram(logLead ~ 1, meuseSf)
# with initial gueses at the parameters psill, range, and nugget
leadModel <- vgm(psill = 0.6, model = "Sph", range = 750, nugget = 0.05)
# or try to let the estimation happen without initial guesses -- same result
leadModel <- vgm(model = "Sph", nugget = TRUE)

leadFit <- fit.variogram(object = leadVar, model = leadModel)
leadGstat <- gstat(
  formula = logLead ~ 1, locations = meuseSf,
  model = leadFit
)
leadKrigeSf <- predict(leadGstat, newdata = meuseGridSf)
leadKrigeSf


## -----------------------------------------------------------------------------
sf2Rast <- function(sfObject, variableIndex = 1) {
  # coerce sf to a data.frame
  dfObject <- data.frame(st_coordinates(sfObject),
    z = as.data.frame(sfObject)[, variableIndex]
  )
  # coerce data.frame to SpatRaster
  rastObject <- rast(dfObject, crs = crs(sfObject))

  names(rastObject) <- names(sfObject)[variableIndex]

  return(rastObject)
}

leadKrigeRast <- sf2Rast(leadKrigeSf)

leadKrigeRast

# and plot
ggplot() +
  geom_spatraster(data = leadKrigeRast, mapping = aes(fill = var1.pred), alpha = 0.8) +
  scale_fill_continuous(type = "viridis", name = "log(ppm)", na.value = "transparent") +
  labs(title = "Lead concentrations") +
  theme_minimal()


## -----------------------------------------------------------------------------
leadKrigeSf$var1.var.sqrt <- sqrt(leadKrigeSf$var1.var)

leadKrigeRast <- sf2Rast(leadKrigeSf, variableIndex = 4)
leadKrigeRast

# and plot
ggplot() +
  geom_spatraster(data = leadKrigeRast, mapping = aes(fill = var1.var.sqrt), alpha = 0.8) +
  scale_fill_continuous(type = "viridis", name = "log(ppm)", na.value = "transparent") +
  labs(title = "Variance of lead concentrations") +
  theme_minimal()


## -----------------------------------------------------------------------------
vgm()


## -----------------------------------------------------------------------------
#| fig-width: 9
#| fig-height: 4
show.vgms(models = c("Exp", "Mat", "Gau", "Sph"))


## -----------------------------------------------------------------------------
leadKrigeLOOCVSf <- krige.cv(
  formula = logLead ~ 1,
  locations = meuseSf,
  model = leadFit, verbose = FALSE
)
leadKrigeLOOCVSf
# CV R2 anbd RMSE
rsq <- cor(leadKrigeLOOCVSf$observed, leadKrigeLOOCVSf$var1.pred)^2
rmse <- sqrt(mean((leadKrigeLOOCVSf$observed - leadKrigeLOOCVSf$var1.pred)^2))
c(rsq = rsq, rmse = rmse)


## -----------------------------------------------------------------------------
library(automap)
leadVar <- autofitVariogram(formula = logLead ~ 1, input_data = meuseSf)
summary(leadVar)
plot(leadVar)


## -----------------------------------------------------------------------------
leadAutoKrigeLOOCV <- autoKrige.cv(
  formula = logLead ~ 1, input_data = meuseSf,
  verbose = c(FALSE, FALSE)
)
summary(leadAutoKrigeLOOCV)
# R2
cor(leadAutoKrigeLOOCV$krige.cv_output$observed, leadAutoKrigeLOOCV$krige.cv_output$var1.pred)^2


## -----------------------------------------------------------------------------
# precip point data
prcpCA <- readRDS("../data/prcpCA.rds")
# empty grid to interpolate into
gridCA <- readRDS("../data/gridCA.rds")

prcpCAsf <- prcpCA %>%
  st_as_sf(coords = c("X", "Y")) %>%
  st_set_crs(value = 3310)

prcpCAsf %>% ggplot() +
  geom_sf(aes(fill = ANNUAL, size = ANNUAL),
    color = "white",
    shape = 21, alpha = 0.8
  ) +
  scale_fill_continuous(type = "viridis", name = "mm") +
  labs(title = "Total Annual Precipitation") +
  scale_size(guide = "none")

