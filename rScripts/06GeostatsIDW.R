## ----echo=FALSE, message=FALSE, warning=FALSE, results='hide'-----------------
#| label: setup
set.seed(184)


## ----message=FALSE------------------------------------------------------------
#| label: packages
library(tidyverse)
library(sf)
library(gstat)
library(terra)
library(tidyterra)


## -----------------------------------------------------------------------------
#| label: weight-curve
wCurve <- data.frame(d = 1:100, w = (1:100)^-1)
wCurve %>% ggplot(mapping = aes(x = d, y = w)) +
  geom_line() +
  labs(x = "Distance", y = "Weight")


## -----------------------------------------------------------------------------
#| label: weight-curve-log
wCurve %>% ggplot(mapping = aes(x = d, y = w)) +
  geom_line() +
  labs(x = "Distance", y = "Weight") +
  scale_y_log10()


## -----------------------------------------------------------------------------
#| label: power-range
wCurve <- rbind(
  data.frame(d = 1:100, w = (1:100)^0, p = "0"),
  data.frame(d = 1:100, w = (1:100)^-0.5, p = "0.5"),
  data.frame(d = 1:100, w = (1:100)^-1, p = "1"),
  data.frame(d = 1:100, w = (1:100)^-1.5, p = "1.5"),
  data.frame(d = 1:100, w = (1:100)^-2, p = "2"),
  data.frame(d = 1:100, w = (1:100)^-2.5, p = "2.5")
)
wCurve %>% ggplot(mapping = aes(x = d, y = w, color = p)) +
  geom_line() +
  labs(x = "Distance", y = "Weight")


## -----------------------------------------------------------------------------
#| label: power-range-log
wCurve %>% ggplot(mapping = aes(x = d, y = w, color = p)) +
  geom_line() +
  labs(x = "Distance", y = "Weight") +
  scale_y_log10()


## -----------------------------------------------------------------------------
#| label: toy-data
#| fig-width: 5
#| fig-height: 5
toyDat <- data.frame(
  x = c(1, 3, 1, 4, 5),
  y = c(5, 4, 3, 5, 1),
  z = c(100, 105, 105, 100, 115)
)
toyDat
toyPlot <- toyDat %>% ggplot() +
  geom_point(aes(x = x, y = y, size = z)) +
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
d2s0 <- as.matrix(dist(cbind(c(toyDat$x, 2), c(toyDat$y, 4))))[1:5, 6]


## ----echo=FALSE---------------------------------------------------------------
#| label: toy-weights
p <- 2
w <- d2s0^-p


## ----echo=FALSE---------------------------------------------------------------
#| label: toy-estimate
zhatS0 <- sum(w * toyDat$z) / sum(w)


## ----warning=FALSE------------------------------------------------------------
#| label: toy-idw
# stack the five known points with the unknown point s0 at (2, 4)
coords <- cbind(c(toyDat$x, 2), c(toyDat$y, 4))
# full pairwise distance matrix, then pull each known point's distance to s0 (column 6)
dmat <- as.matrix(dist(coords))
d2s0 <- dmat[1:5, 6]
# power
p <- 2
# weights
w <- d2s0^-p
# and the estimation itself
zhatS0 <- sum(w * toyDat$z) / sum(w)
zhatS0
# add it to the plot
toyPlot + geom_point(aes(x = 2, y = 4, size = zhatS0))


## -----------------------------------------------------------------------------
#| label: meuse-load
meuse2 <- readRDS("../data/meuse2.Rds")
glimpse(meuse2)
class(meuse2)
meuse2$logLead <- log(meuse2$lead)
# or for the tidyverse fans this is the same output
meuse2 <- meuse2 %>% mutate(logLead = log(lead))
# make into sf
meuseSf <- st_as_sf(meuse2, coords = c("x", "y")) %>%
  st_set_crs(value = 28992)

class(meuseSf) # note change in class from data.frame to sf and data.frame

p2 <- ggplot(data = meuseSf) +
  geom_sf(aes(fill = logLead),
    size = 4,
    shape = 21, color = "white", alpha = 0.8
  ) +
  scale_fill_continuous(type = "viridis", name = "log(ppm)") +
  labs(title = "Lead concentrations")
p2


## -----------------------------------------------------------------------------
#| label: meuse-grid-load
meuse.grid2 <- readRDS("../data/meuse.grid2.Rds")
head(meuse.grid2)


## -----------------------------------------------------------------------------
#| label: meuse-grid-sf
meuseGridSf <- st_as_sf(meuse.grid2,
  coords = c("x", "y"),
  crs = st_crs(meuseSf)
)
meuseGridSf


## -----------------------------------------------------------------------------
#| label: idw-p2-model
idwP2Model <- gstat(
  formula = logLead ~ 1,
  locations = meuseSf,
  set = list(idp = 2)
)
logLeadIDWP2Sf <- predict(idwP2Model, meuseGridSf)
logLeadIDWP2Sf


## -----------------------------------------------------------------------------
#| label: idw-p2-map
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

logLeadIDWP2Rast <- sf2Rast(logLeadIDWP2Sf)
logLeadIDWP2Rast
# and plot
ggplot() +
  geom_spatraster(data = logLeadIDWP2Rast, mapping = aes(fill = var1.pred), alpha = 0.8) +
  scale_fill_continuous(type = "viridis", name = "log(ppm)", na.value = "transparent") +
  labs(title = "Lead concentrations", subtitle = "IDW with p=2") +
  theme_minimal()


## -----------------------------------------------------------------------------
#| label: idw-p2-contour
ggplot() +
  geom_spatraster_contour_filled(
    data = logLeadIDWP2Rast,
    breaks = seq(from = 3.5, to = 6.5, by = 0.25),
    alpha = 0.9
  ) +
  scale_fill_discrete(name = "log(ppm)", na.value = "transparent") +
  labs(title = "Lead concentrations", subtitle = "IDW with p=2") +
  theme_minimal()


## -----------------------------------------------------------------------------
#| label: idw-insample-skill
obs <- meuseSf$logLead
preds <- extract(logLeadIDWP2Rast, meuseSf) %>% pull(var1.pred)
rsq <- cor(obs, preds)^2
rmse <- sqrt(mean((preds - obs)^2))
rsq
rmse
ggplot() +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(aes(x = obs, y = preds)) +
  coord_cartesian() +
  labs(
    x = "Observed Values",
    y = "Predicted Values",
    title = "Lead log(ppm)"
  )


## -----------------------------------------------------------------------------
#| label: idw-traintest
n <- nrow(meuseSf)
rows4test <- sample(x = 1:n, size = n * 0.25)
meuseTest <- meuseSf[rows4test, ]
meuseTrain <- meuseSf[-rows4test, ]

# note that we build the model with meuseTrain
idwP2Model <- gstat(
  formula = logLead ~ 1,
  locations = meuseTrain,
  set = list(idp = 2)
)
logLeadIDWP2Sf <- predict(idwP2Model, meuseGridSf)

logLeadIDWP2Rast <- sf2Rast(logLeadIDWP2Sf)

# and plot
ggplot() +
  geom_spatraster(data = logLeadIDWP2Rast, mapping = aes(fill = var1.pred), alpha = 0.8) +
  scale_fill_continuous(type = "viridis", name = "log(ppm)", na.value = "transparent") +
  labs(title = "Lead concentrations", subtitle = "IDW with p=2") +
  theme_minimal()


## -----------------------------------------------------------------------------
#| label: idw-oos-skill
# note use of meuseTest here
obs <- meuseTest$logLead
preds <- extract(logLeadIDWP2Rast, meuseTest) %>% pull(var1.pred)
rsq <- cor(obs, preds)^2
rmse <- sqrt(mean((preds - obs)^2))
rsq
rmse

ggplot() +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(aes(x = obs, y = preds)) +
  coord_fixed(ratio = 1, xlim = range(preds, obs), ylim = range(preds, obs)) +
  labs(
    x = "Observed Values",
    y = "Predicted Values",
    title = "Lead log(ppm)"
  )


## -----------------------------------------------------------------------------
#| label: null-model
rmseNULL <- sqrt(mean((mean(meuseTrain$logLead) - obs)^2))
rmseNULL
1 - (rmse / rmseNULL)


## -----------------------------------------------------------------------------
#| label: idw-p3-map
idwP3Model <- gstat(
  formula = logLead ~ 1,
  locations = meuseTrain,
  set = list(idp = 3)
)

logLeadIDWP3Sf <- predict(idwP3Model, meuseGridSf)

logLeadIDWP3Rast <- sf2Rast(logLeadIDWP3Sf)

# and plot
ggplot() +
  geom_spatraster(data = logLeadIDWP3Rast, mapping = aes(fill = var1.pred), alpha = 0.8) +
  scale_fill_continuous(type = "viridis", name = "log(ppm)", na.value = "transparent") +
  labs(title = "Lead concentrations") +
  theme_minimal()


## -----------------------------------------------------------------------------
#| label: idw-p3-skill
# note use of meuseTest here
obs <- meuseTest$logLead
preds <- extract(logLeadIDWP3Rast, meuseTest) %>% pull(var1.pred)
rsq <- cor(obs, preds)^2
rmse <- sqrt(mean((preds - obs)^2))
rsq
rmse
ggplot() +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(aes(x = obs, y = preds)) +
  coord_fixed(ratio = 1, xlim = range(preds, obs), ylim = range(preds, obs)) +
  labs(
    x = "Observed Values",
    y = "Predicted Values",
    title = "Lead log(ppm)"
  )


## -----------------------------------------------------------------------------
#| label: idw-p3-skillscore
1 - (rmse / rmseNULL)


## -----------------------------------------------------------------------------
#| label: prcp-ca
# precip point data
prcpCA <- readRDS("../data/prcpCA.rds")
# empty grid to interpolate into
gridCA <- readRDS("../data/gridCA.rds")

# make as sf
prcpCA <- prcpCA %>%
  st_as_sf(coords = c("X", "Y")) %>%
  st_set_crs(value = 3310)

# simple map
prcpCA %>% ggplot() +
  geom_sf(aes(fill = ANNUAL, size = ANNUAL),
    color = "white",
    shape = 21, alpha = 0.8
  ) +
  scale_fill_continuous(type = "viridis", name = "mm") +
  labs(title = "Total Annual Precipitation") +
  scale_size(guide = "none")

