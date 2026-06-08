## ----echo=FALSE, include=FALSE------------------------------------------------
set.seed(130)


## ----message=FALSE------------------------------------------------------------
library(spatialreg)
library(spdep)
library(sf)
library(tidyverse)


## ----message=FALSE, warning=FALSE---------------------------------------------
n        <- 75
easting  <- runif(n, 0, 100)
northing <- runif(n, 0, 100)
points   <- cbind(easting, northing)

# Build spatially autocorrelated errors via a SAR process
dnb   <- dnearneigh(x = points, d1 = 0, d2 = 150)
dsts  <- nbdists(nb = dnb, coords = points)
p     <- 2.25
idw   <- lapply(dsts, function(x) x^-p)
wList <- nb2listw(neighbours = dnb, glist = idw, style = "W")

inv     <- spatialreg::invIrW(x = wList, rho = 0.75)
epsilon <- inv %*% rnorm(n)
epsilon <- scale(epsilon[,1])[,1]

# Response variable: y is a function of x plus autocorrelated noise
x  <- rnorm(n)
B0 <- 3
B1 <- 0.5
y  <- B0 + B1*x + epsilon

dat <- data.frame(easting, northing, y, x)


## -----------------------------------------------------------------------------
nb_k8 <- knn2nb(knearneigh(points, k = 8))
W <- nb2listw(nb_k8, style = "W")


## -----------------------------------------------------------------------------
ols <- lm(y ~ x, data = dat)
dat$olsResids <- residuals(ols)
moran.test(dat$olsResids, W)


## -----------------------------------------------------------------------------
lm.RStests(ols, listw = W, test = "all")


## -----------------------------------------------------------------------------
semFit <- errorsarlm(y ~ x, data = dat, listw = W)
summary(semFit)


## -----------------------------------------------------------------------------
dat$semResids <- residuals(semFit)
moran.test(dat$semResids, W)


## -----------------------------------------------------------------------------
slmFit <- lagsarlm(y ~ x, data = dat, listw = W)
summary(slmFit)


## -----------------------------------------------------------------------------
AIC(semFit, slmFit)


## ----message=FALSE------------------------------------------------------------
library(tmap)
birds_sf <- readRDS("../data/birdRichnessMexico.rds")
coords_birds <- st_coordinates(birds_sf)

nb_birds <- knn2nb(knearneigh(coords_birds, k = 8))
W_birds  <- nb2listw(nb_birds, style = "W")


## -----------------------------------------------------------------------------
tmap_mode("view")
tm_shape(birds_sf) +
  tm_symbols(col = "nSpecies", palette = "viridis",
             title.col = "Species", size = 0.4)


## -----------------------------------------------------------------------------
ols_birds <- lm(nSpecies ~ map + tempRange, data = birds_sf)
summary(ols_birds)


## -----------------------------------------------------------------------------
birds_sf$ols_resids <- residuals(ols_birds)
moran.test(birds_sf$ols_resids, W_birds)


## -----------------------------------------------------------------------------
tm_shape(birds_sf) +
  tm_symbols(col = "ols_resids",
             palette = "-RdBu",
             midpoint = 0,
             title.col = "Residual",
             size = 0.4)


## -----------------------------------------------------------------------------
lm.RStests(ols_birds, listw = W_birds, test = "all")


## -----------------------------------------------------------------------------
sem_birds <- errorsarlm(nSpecies ~ map + tempRange,
                        data = birds_sf, listw = W_birds)
summary(sem_birds)


## -----------------------------------------------------------------------------
birds_sf$sem_resids <- residuals(sem_birds)
moran.test(birds_sf$sem_resids, W_birds)


## -----------------------------------------------------------------------------
AIC(ols_birds, sem_birds)


## -----------------------------------------------------------------------------
impacts(slmFit, listw = W, R = 500)


## ----eval=FALSE---------------------------------------------------------------
# meuse2      <- readRDS("../data/meuse2.Rds")
# meuse.grid2 <- readRDS("../data/meuse.grid2.Rds")
# meuse_sf <- st_as_sf(meuse2, coords = c("x", "y"), crs = 28992)
# meuse_sf$log_lead <- log(meuse_sf$lead)
# # get river distance and flooding frequency from the grid
# covars_sf <- st_as_sf(meuse.grid2, coords = c("x","y"), crs = 28992) %>%
#   select(ffreq, river_dist_m)
# meuse_sf <- st_join(meuse_sf, covars_sf, join = st_nearest_feature) %>%
#   mutate(ffreq = factor(ffreq))

