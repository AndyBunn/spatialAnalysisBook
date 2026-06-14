## ----echo=FALSE, message=FALSE, warning=FALSE, results='hide'-----------------
#| label: setup
set.seed(184)


## ----message=FALSE------------------------------------------------------------
#| label: packages
library(sf)
library(gstat)
library(tidyverse)
library(terra)
library(tidyterra)
library(automap)


## -----------------------------------------------------------------------------
#| label: meuse-data
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

ggplot(data = meuseSf) +
  geom_sf(aes(fill = logLead),
    size = 4,
    shape = 21, color = "white", alpha = 0.8
  ) +
  scale_fill_continuous(type = "viridis", name = "log(ppm)") +
  labs(title = "Lead concentrations") +
  theme_minimal()


## -----------------------------------------------------------------------------
#| label: lead-variogram
leadVar <- variogram(logLead ~ 1, meuseSf)
plot(leadVar,
  pch = 20, cex = 1.5, col = "black",
  ylab = expression("Semivariance (" * gamma * ")"),
  xlab = "Distance (m)", main = "Lead concentrations (log(ppm))"
)


## -----------------------------------------------------------------------------
#| label: variogram-anatomy
#| echo: false
#| fig-width: 7
#| fig-height: 5

# a schematic spherical variogram to label the three terms
nugget <- 0.2
psill <- 0.8 # the partial sill, sill minus nugget
sill <- nugget + psill
rng <- 10

sphSchematic <- function(d) {
  ifelse(d < rng,
    nugget + psill * (1.5 * (d / rng) - 0.5 * (d / rng)^3),
    sill
  )
}

vgCurve <- data.frame(d = seq(0, 18, length.out = 400))
vgCurve$gamma <- sphSchematic(vgCurve$d)

ggplot(vgCurve, aes(d, gamma)) +
  # reference lines for sill, nugget, and range
  geom_hline(yintercept = sill, linetype = "dashed", color = "grey55") +
  geom_hline(yintercept = nugget, linetype = "dashed", color = "grey55") +
  geom_vline(xintercept = rng, linetype = "dashed", color = "grey55") +
  # the model curve and the nugget intercept
  geom_line(linewidth = 1.3, color = "#440154") +
  annotate("point", x = 0, y = nugget, color = "#440154", size = 3) +
  # brackets for nugget and partial sill, parked in the flat region
  annotate("segment",
    x = 15.5, xend = 15.5, y = 0, yend = nugget,
    arrow = arrow(ends = "both", length = unit(0.07, "in")), color = "grey30"
  ) +
  annotate("segment",
    x = 15.5, xend = 15.5, y = nugget, yend = sill,
    arrow = arrow(ends = "both", length = unit(0.07, "in")), color = "grey30"
  ) +
  # labels
  annotate("text", x = 16, y = nugget / 2, label = "nugget", hjust = 0, size = 4) +
  annotate("text", x = 16, y = (nugget + sill) / 2, label = "partial\nsill", hjust = 0, size = 4) +
  annotate("text", x = 0.3, y = sill + 0.05, label = "sill", hjust = 0, size = 4) +
  annotate("text", x = rng + 0.3, y = 0.05, label = "range", hjust = 0, size = 4) +
  scale_x_continuous(breaks = NULL, expand = expansion(mult = c(0.01, 0.12))) +
  scale_y_continuous(breaks = NULL, limits = c(0, sill * 1.15)) +
  labs(
    x = "Distance", y = expression("Semivariance (" * gamma * ")"),
    title = "Anatomy of a variogram"
  ) +
  theme_minimal()


## -----------------------------------------------------------------------------
#| label: spherical-fit
# note our initial estimates for the partial sill, range, and nugget.
sphericalModel <- vgm(psill = 0.5, model = "Sph", range = 750, nugget = 0.05)
sphericalFit <- fit.variogram(object = leadVar, model = sphericalModel)
sphericalFit # look at the fitted values


## -----------------------------------------------------------------------------
plot(leadVar,
  model = sphericalFit, pch = 20, cex = 1.5, col = "black",
  ylab = expression("Semivariance (" * gamma * ")"),
  xlab = "Distance (m)", main = "Lead concentrations (log(ppm))",
  sub = "Points: Empirical, Line: Spherical Model"
)


## -----------------------------------------------------------------------------
#| label: exponential-fit
# note our initial estimates for the sill, range, and nugget
exponentialModel <- vgm(psill = 0.5, model = "Exp", range = 750, nugget = 0.05)
exponentialFit <- fit.variogram(object = leadVar, model = exponentialModel)
plot(leadVar,
  model = exponentialFit, pch = 20, cex = 1.5, col = "black",
  ylab = expression("Semivariance (" * gamma * ")"),
  xlab = "Distance (m)", main = "Lead concentrations (log(ppm))",
  sub = "Points: Empirical, Line: Exponential Model"
)


## -----------------------------------------------------------------------------
#| label: toy-data
#| fig-width: 5
#| fig-height: 5
foo <- data.frame(
  x = c(1, 3, 1, 4, 5),
  y = c(5, 4, 3, 5, 1),
  z = c(100, 105, 105, 100, 115)
)
foo

toyPlot <- ggplot() +
  geom_point(data = foo, aes(x = x, y = y, size = z)) +
  lims(x = c(0, 6), y = c(0, 6)) +
  coord_equal() +
  theme_minimal()
toyPlot


## -----------------------------------------------------------------------------
#| label: toy-unknown-point
#| fig-width: 5
#| fig-height: 5
toyPlot +
  # crosshair reticle marking the location we want to predict
  annotate("segment", x = 1.4, xend = 2.6, y = 4, yend = 4, color = "#21908C", linewidth = 0.6) +
  annotate("segment", x = 2, xend = 2, y = 3.4, yend = 4.6, color = "#21908C", linewidth = 0.6) +
  annotate("point", x = 2, y = 4, shape = 1, size = 13, color = "#21908C", stroke = 1) +
  annotate("point", x = 2, y = 4, shape = 21, fill = "#21908C", color = "white", size = 2.5, stroke = 0.6) +
  annotate("text", x = 2.55, y = 4.6, label = "s[0]", parse = TRUE, color = "#21908C", size = 5, hjust = 0, fontface = "italic")


## ----echo=FALSE---------------------------------------------------------------
#| label: toy-distances
d2s0 <- as.matrix(dist(cbind(c(foo$x, 2), c(foo$y, 4))))[1:5, 6]


## ----echo=FALSE---------------------------------------------------------------
#| label: toy-g
sphGamma <- function(d, psill = 10, range = 6) {
  psill * (1.5 * (d / range) - 0.5 * (d / range)^3)
}
g <- sphGamma(d2s0)


## ----echo=FALSE---------------------------------------------------------------
#| label: toy-gamma
dMat <- as.matrix(dist(cbind(foo$x, foo$y)))
G <- sphGamma(dMat)
diag(G) <- 0


## ----echo=FALSE---------------------------------------------------------------
#| label: toy-lambda
lambda <- solve(G) %*% g
zhat <- sum(lambda[, 1] * foo$z)


## ----echo=FALSE---------------------------------------------------------------
#| label: toy-prediction
0.2626 * (100) + 0.4985 * (105) + 0.2652 * (105) - 0.0165 * (100) - 0.0353 * (115)


## -----------------------------------------------------------------------------
#| label: lead-krige
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
#| label: lead-krige-surface
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
#| label: lead-krige-variance
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
#| label: vgm-table
vgm()


## -----------------------------------------------------------------------------
#| label: show-vgms
#| fig-width: 9
#| fig-height: 4
show.vgms(models = c("Exp", "Mat", "Gau", "Sph"))


## -----------------------------------------------------------------------------
#| label: lead-loocv
leadKrigeLOOCVSf <- krige.cv(
  formula = logLead ~ 1,
  locations = meuseSf,
  model = leadFit, verbose = FALSE
)
leadKrigeLOOCVSf
# CV R2 and RMSE
rsq <- cor(leadKrigeLOOCVSf$observed, leadKrigeLOOCVSf$var1.pred)^2
rmse <- sqrt(mean((leadKrigeLOOCVSf$observed - leadKrigeLOOCVSf$var1.pred)^2))
c(rsq = rsq, rmse = rmse)


## -----------------------------------------------------------------------------
#| label: automap-variogram
leadVar <- autofitVariogram(formula = logLead ~ 1, input_data = meuseSf)
summary(leadVar)
plot(leadVar)


## -----------------------------------------------------------------------------
#| label: automap-loocv
leadAutoKrigeLOOCV <- autoKrige.cv(
  formula = logLead ~ 1, input_data = meuseSf,
  verbose = c(FALSE, FALSE)
)
summary(leadAutoKrigeLOOCV)
# R2
cor(leadAutoKrigeLOOCV$krige.cv_output$observed, leadAutoKrigeLOOCV$krige.cv_output$var1.pred)^2


## -----------------------------------------------------------------------------
#| label: exercise-data
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

