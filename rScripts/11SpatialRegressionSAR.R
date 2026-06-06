## ----echo=FALSE, include=FALSE------------------------------------------------
set.seed(2112)


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
lm.LMtests(ols, listw = W, test = "all")


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


## ----eval=FALSE---------------------------------------------------------------
# dat_birds <- read.csv("../data/birdDiv.csv")
# datSF <- dat_birds %>%
#   st_as_sf(coords = c("UTME", "UTMN"), crs = 26912)
# 
# coords_birds <- st_coordinates(datSF)
# nb_birds     <- knn2nb(knearneigh(coords_birds, k = 8))
# W_birds      <- nb2listw(nb_birds, style = "W")


## ----eval=FALSE---------------------------------------------------------------
# ols_birds <- lm(birdDiv ~ plantDiv, data = dat_birds)
# moran.test(residuals(ols_birds), W_birds)


## ----eval=FALSE---------------------------------------------------------------
# lm.LMtests(ols_birds, listw = W_birds, test = "all")


## ----eval=FALSE---------------------------------------------------------------
# impacts(slmFit, listw = W, R = 500)


## ----eval=FALSE---------------------------------------------------------------
# data(meuse.all)
# meuse_sf <- st_as_sf(meuse.all, coords = c("x", "y")) %>%
#   st_set_crs(28992)
# meuse_sf$log_lead <- log(meuse_sf$lead)


## ----echo=FALSE, eval=FALSE---------------------------------------------------
# # --- instructor key ---
# data(meuse.all)
# meuse_sf <- st_as_sf(meuse.all, coords = c("x", "y")) %>%
#   st_set_crs(28992)
# meuse_sf$log_lead <- log(meuse_sf$lead)
# coords <- st_coordinates(meuse_sf)
# 
# nb_k8 <- knn2nb(knearneigh(coords, k = 8))
# W <- nb2listw(nb_k8, style = "W")
# 
# # OLS
# ols_lead <- lm(log_lead ~ dist + ffreq, data = as.data.frame(meuse_sf))
# moran.test(residuals(ols_lead), W)  # significant
# 
# # LM tests -- expect LMerr to dominate
# lm.LMtests(ols_lead, listw = W, test = "all")
# 
# # SEM
# sem_lead <- errorsarlm(log_lead ~ dist + ffreq,
#                        data = as.data.frame(meuse_sf), listw = W)
# summary(sem_lead)
# moran.test(residuals(sem_lead), W)  # should be clean
# 
# # SLM for comparison
# slm_lead <- lagsarlm(log_lead ~ dist + ffreq,
#                      data = as.data.frame(meuse_sf), listw = W)
# summary(slm_lead)
# 
# AIC(sem_lead, slm_lead)
# 
# # The lag model would imply that lead contamination at one location
# # spatially drives contamination at neighboring locations --
# # a diffusion or contagion process. For inert fluvial sediment,
# # that mechanism is not ecologically compelling. The error model
# # (shared unmeasured drivers: flood frequency variation, sediment
# # texture, proximity to historical channel features) is more defensible.

