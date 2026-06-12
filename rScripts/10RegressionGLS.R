## ----message=FALSE------------------------------------------------------------
library(automap)
library(gstat)
library(nlme)
library(spdep)
library(ncf)
library(tidyverse)


## ----echo=FALSE---------------------------------------------------------------
set.seed(3515) # 234 #104


## ----message=FALSE, warning=FALSE---------------------------------------------
# We will do n points
n <- 75

# CSR. Our mapped points.
easting <- runif(n, 0, 100)
northing <- runif(n, 0, 100)
points <- cbind(easting, northing)

# Generate the spatially autocorrelated error term, epsilon.
# Build neighbor objects with an inclusive distance threshold
dnb <- dnearneigh(x = points, d1 = 0, d2 = 150)
dsts <- nbdists(nb = dnb, coords = points)
# And follow an IDW process
p <- 2.25 # IDW power
idw <- lapply(dsts, function(x) x^-p)
# Weights matrix as a list
wList <- nb2listw(neighbours = dnb, glist = idw, style = "W")

# Compute SAR generating operator via a weights matrix
inv <- spatialreg::invIrW(x = wList, rho = 0.75)

# Calculate epsilon and rescale as z scores
epsilon <- inv %*% rnorm(n)
epsilon <- scale(epsilon[, 1])[, 1]


## -----------------------------------------------------------------------------
epsilonCorrelog <- spline.correlog(
  x = easting, y = northing,
  z = epsilon, resamp = 50, quiet = TRUE
)
plot(epsilonCorrelog, xlim = c(0, max(dist(points)) / 3))


## -----------------------------------------------------------------------------
# Generate x as noise.
x <- rnorm(n)
# Generate y as a function of x and epsilon
# (remember epsilon is spatially autocorrelated but we don't know that)
B0 <- 3 # intercept
B1 <- 0.5 # slope
y <- B0 + B1 * x + epsilon

# Make a data.frame
dat <- data.frame(easting, northing, y, x)


## -----------------------------------------------------------------------------
ggplot(data = dat, mapping = aes(
  x = easting, y = northing,
  size = x, fill = x
)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c() +
  guides(fill = guide_legend(), size = guide_legend()) +
  coord_fixed() +
  ggtitle("x")


## -----------------------------------------------------------------------------
ggplot(data = dat, mapping = aes(
  x = easting, y = northing,
  size = y, fill = y
)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c() +
  guides(fill = guide_legend(), size = guide_legend()) +
  coord_fixed() +
  ggtitle("y")


## -----------------------------------------------------------------------------
# Fit a model using gls
glsNaive <- gls(y ~ x, dat)
summary(glsNaive)


## -----------------------------------------------------------------------------
# Add residuals to dat
dat$glsNaiveResids <- residuals(glsNaive, type = "normalized")
ggplot(data = dat, mapping = aes(
  x = easting, y = northing,
  size = glsNaiveResids, fill = glsNaiveResids
)) +
  geom_point(shape = 21) +
  scale_fill_gradient2() +
  guides(fill = guide_legend(), size = guide_legend()) +
  coord_fixed()


## -----------------------------------------------------------------------------
# Test them for spatial autocorrelation using a correlogram.
residsI <- spline.correlog(
  x = dat$easting, y = dat$northing,
  z = dat$glsNaiveResids, resamp = 50, quiet = TRUE
)
plot(residsI, xlim = c(0, max(dist(points)) / 3))


## -----------------------------------------------------------------------------
datSF <- dat %>% st_as_sf(coords = c("easting", "northing"))
# A variogram of the residuals
plot(autofitVariogram(glsNaiveResids ~ 1, input_data = datSF, model = c("Gau", "Sph", "Exp")))


## -----------------------------------------------------------------------------
csSpatial <- corSpatial(form = ~ easting + northing, nugget = TRUE, type = "spherical")
glsUpdated <- update(glsNaive, correlation = csSpatial)
summary(glsUpdated)


## -----------------------------------------------------------------------------
# add the residuals to the spatial dat object
dat$glsResids <- residuals(glsUpdated, type = "normalized")
# Map the residuals from gls
ggplot(data = dat, mapping = aes(
  x = easting, y = northing,
  size = glsResids, fill = glsResids
)) +
  geom_point(shape = 21) +
  scale_fill_gradient2() +
  guides(fill = guide_legend(), size = guide_legend()) +
  coord_fixed()


## -----------------------------------------------------------------------------
residsI <- spline.correlog(
  x = dat$easting, y = dat$northing,
  z = dat$glsResids, resamp = 50, quiet = TRUE
)
plot(residsI, xlim = c(0, max(dist(points)) / 3))

