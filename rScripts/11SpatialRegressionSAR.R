## ----echo=FALSE, include=FALSE------------------------------------------------
#| label: setup
set.seed(130)


## ----message=FALSE------------------------------------------------------------
#| label: packages
library(spatialreg)
library(spdep)
library(sf)
library(tidyverse)


## ----message=FALSE, warning=FALSE---------------------------------------------
#| label: toy-data
n <- 75
easting <- runif(n, 0, 100)
northing <- runif(n, 0, 100)
points <- cbind(easting, northing)

# Build spatially autocorrelated errors via a SAR process
dnb <- dnearneigh(x = points, d1 = 0, d2 = 150)
dsts <- nbdists(nb = dnb, coords = points)
p <- 2.25
idw <- lapply(dsts, function(x) x^-p)
wList <- nb2listw(neighbours = dnb, glist = idw, style = "W")

inv <- spatialreg::invIrW(x = wList, rho = 0.75)
epsilon <- inv %*% rnorm(n)
epsilon <- scale(epsilon[, 1])[, 1]

# Response variable: y is a function of x plus autocorrelated noise
x <- rnorm(n)
B0 <- 3
B1 <- 0.5
y <- B0 + B1 * x + epsilon

dat <- data.frame(easting, northing, y, x)


## -----------------------------------------------------------------------------
#| label: toy-weights
nbK8 <- knn2nb(knearneigh(points, k = 8))
W <- nb2listw(nbK8, style = "W")


## -----------------------------------------------------------------------------
#| label: toy-ols
ols <- lm(y ~ x, data = dat)
dat$olsResids <- residuals(ols)
moran.test(dat$olsResids, W)


## -----------------------------------------------------------------------------
#| label: toy-lmtests
lm.RStests(ols, listw = W, test = "all")


## -----------------------------------------------------------------------------
#| label: toy-sem
semFit <- errorsarlm(y ~ x, data = dat, listw = W)
summary(semFit)


## -----------------------------------------------------------------------------
#| label: toy-sem-resids
dat$semResids <- residuals(semFit)
moran.test(dat$semResids, W)


## -----------------------------------------------------------------------------
#| label: toy-slm
slmFit <- lagsarlm(y ~ x, data = dat, listw = W)
summary(slmFit)


## -----------------------------------------------------------------------------
#| label: toy-aic
AIC(semFit, slmFit)


## ----message=FALSE------------------------------------------------------------
#| label: bird-data
library(tmap)
birdsSf <- readRDS("../data/birdRichnessMexico.rds")
coordsBirds <- st_coordinates(birdsSf)

nbBirds <- knn2nb(knearneigh(coordsBirds, k = 8))
WBirds <- nb2listw(nbBirds, style = "W")


## -----------------------------------------------------------------------------
#| label: bird-richness-map
tmap_mode("plot")
tm_shape(birdsSf) +
  tm_symbols(
    col = "nSpecies", palette = "viridis",
    title.col = "Species", size = 0.4
  )


## -----------------------------------------------------------------------------
#| label: bird-ols
olsBirds <- lm(nSpecies ~ map + tempRange, data = birdsSf)
summary(olsBirds)


## -----------------------------------------------------------------------------
#| label: bird-ols-moran
birdsSf$ols_resids <- residuals(olsBirds)
moran.test(birdsSf$ols_resids, WBirds)


## -----------------------------------------------------------------------------
#| label: bird-resid-map
tm_shape(birdsSf) +
  tm_symbols(
    col = "ols_resids",
    palette = "-RdBu",
    midpoint = 0,
    title.col = "Residual",
    size = 0.4
  )


## -----------------------------------------------------------------------------
#| label: bird-lmtests
lm.RStests(olsBirds, listw = WBirds, test = "all")


## -----------------------------------------------------------------------------
#| label: bird-sem
semBirds <- errorsarlm(nSpecies ~ map + tempRange,
  data = birdsSf, listw = WBirds
)
summary(semBirds)


## -----------------------------------------------------------------------------
#| label: bird-sem-moran
birdsSf$sem_resids <- residuals(semBirds)
moran.test(birdsSf$sem_resids, WBirds)


## -----------------------------------------------------------------------------
#| label: bird-aic
AIC(olsBirds, semBirds)


## -----------------------------------------------------------------------------
#| label: slm-impacts
impacts(slmFit, listw = W, R = 500)


## ----eval=FALSE---------------------------------------------------------------
#| label: exercise-data
# meuse2 <- readRDS("../data/meuse2.Rds")
# meuse.grid2 <- readRDS("../data/meuse.grid2.Rds")
# meuseSf <- st_as_sf(meuse2, coords = c("x", "y"), crs = 28992)
# meuseSf$log_lead <- log(meuseSf$lead)
# # get river distance and flooding frequency from the grid
# covarsSf <- st_as_sf(meuse.grid2, coords = c("x", "y"), crs = 28992) %>%
#   select(ffreq, river_dist_m)
# meuseSf <- st_join(meuseSf, covarsSf, join = st_nearest_feature) %>%
#   mutate(ffreq = factor(ffreq))

