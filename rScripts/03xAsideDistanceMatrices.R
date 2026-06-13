## -----------------------------------------------------------------------------
#| label: setup
#| echo: false
set.seed(1984)


## -----------------------------------------------------------------------------
#| label: packages
#| message: false
library(tidyverse)
library(sf)


## -----------------------------------------------------------------------------
#| label: sample-points
dat <- data.frame(
  x = c(1, 3, 5, 2, 6, 4),
  y = c(1, 4, 2, 6, 5, 3)
)


## -----------------------------------------------------------------------------
#| label: sample-points-plot
#| echo: false
ggplot(dat, aes(x = x, y = y, label = 1:nrow(dat))) +
  geom_point(size = 4) +
  geom_text(nudge_y = 0.3, size = 3.5) +
  coord_fixed() +
  theme_minimal() +
  labs(title = "Six sample locations")


## -----------------------------------------------------------------------------
#| label: dist-object
D <- dist(dat)
D


## -----------------------------------------------------------------------------
#| label: class-dist
class(D)


## -----------------------------------------------------------------------------
#| label: dist-matrix
Dmat <- as.matrix(D)
Dmat


## -----------------------------------------------------------------------------
#| label: class-compare
class(Dmat)
class(D)


## -----------------------------------------------------------------------------
#| label: verify-symmetry
# Verify symmetry
all(Dmat == t(Dmat))

# Verify diagonal
all(diag(Dmat) == 0)


## -----------------------------------------------------------------------------
#| label: nearest-neighbor
Dnn <- Dmat
diag(Dnn) <- NA

# Index of the nearest neighbor for each point
nn1 <- apply(Dnn, 1, which.min)
nn1


## -----------------------------------------------------------------------------
#| label: knn
# 2 nearest neighbors for each point
nn2 <- apply(Dnn, 1, function(row) order(row)[1:2])
nn2


## -----------------------------------------------------------------------------
#| label: stream-cost
# cost[i,j] = effective distance traveling FROM site i TO site j
# Sites: 1=downstream, 2=mid, 3=upstream
# Upstream travel (against current) penalized 4x

streamCost <- matrix(
  c(
    0, 300 * 4, 700 * 4, # from site 1: going up to 2 or 3 is hard
    300, 0, 400 * 4, # from site 2: going down to 1 is easy, up to 3 is hard
    700, 400, 0 # from site 3: going down to 1 or 2 is easy
  ),
  nrow = 3, byrow = TRUE,
  dimnames = list(c("down", "mid", "up"), c("down", "mid", "up"))
)

streamCost


## -----------------------------------------------------------------------------
#| label: check-symmetry
isSymmetric(streamCost) # FALSE -- as expected
isSymmetric(Dmat) # TRUE  -- Euclidean is always symmetric


## -----------------------------------------------------------------------------
#| label: pair-counts
nVals <- c(10, 100, 500, 1000, 5000, 10000)
nPairs <- nVals * (nVals - 1) / 2

data.frame(n = nVals, pairs = format(nPairs, big.mark = ","))


## -----------------------------------------------------------------------------
#| label: dist-memory
sizes <- sapply(c(100, 500, 1000, 2000), function(n) {
  fake <- data.frame(x = runif(n), y = runif(n))
  object.size(dist(fake))
})

data.frame(
  n = c(100, 500, 1000, 2000),
  size = paste(round(sizes / 1e6, 2), "MB")
)


## -----------------------------------------------------------------------------
#| label: st-distance-compare
# First approach: dist() on raw coordinates
DDist <- dist(dat)

# Second approach: st_distance() on an sf object
datSf <- st_as_sf(dat, coords = c("x", "y"), crs = 32610) # UTM Zone 10N, meters
DSf <- st_distance(datSf)

# Compare: upper-left 3x3 corner
round(as.matrix(DDist)[1:3, 1:3], 3)
round(DSf[1:3, 1:3], 3)


## -----------------------------------------------------------------------------
#| label: mtbaker-sites
sites <- data.frame(
  name = c("Bellingham", "Glacier", "Concrete", "Anacortes"),
  lon  = c(-122.48, -121.93, -121.75, -122.61),
  lat  = c(48.75, 48.89, 48.54, 48.51)
)


## -----------------------------------------------------------------------------
#| label: latlon-distances
DLatlon <- as.matrix(dist(sites[, c("lon", "lat")]))
round(DLatlon, 3)


## -----------------------------------------------------------------------------
#| label: projected-distances
sitesSf <- st_as_sf(sites, coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(32610)

DProj <- st_distance(sitesSf)
rownames(DProj) <- sites$name
colnames(DProj) <- sites$name
round(DProj / 1000, 1) # convert to km

