## -----------------------------------------------------------------------------
#| label: setup
#| echo: false
#| include: false
set.seed(1984)


## -----------------------------------------------------------------------------
#| label: packages
#| message: false
library(sf)
library(tidyverse)
library(spdep)
library(tmap)


## -----------------------------------------------------------------------------
#| label: meuse-load
meuse2 <- readRDS("../data/meuse2.Rds")
meuseSf <- st_as_sf(meuse2, coords = c("x", "y")) %>%
  st_set_crs(value = 28992)
meuseSf$log_lead <- log(meuseSf$lead)


## -----------------------------------------------------------------------------
#| label: knn-weights
knn <- knearneigh(meuseSf, k = 8)
nb <- knn2nb(knn)
wList <- nb2listw(nb, style = "W")


## -----------------------------------------------------------------------------
#| label: moran-scatterplot
# standardize log lead
meuseSf$z <- as.vector(scale(meuseSf$log_lead))

# compute the spatial lag: weighted average of neighbors' z values
meuseSf$lag_z <- lag.listw(wList, meuseSf$z)

ggplot(meuseSf, aes(x = z, y = lag_z)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.7) +
  geom_smooth(
    method = "lm", formula = y ~ x - 1,
    se = FALSE, color = "firebrick"
  ) +
  labs(
    x = "Standardized log lead (z)",
    y = "Spatial lag of z",
    title = "Moran scatterplot: log lead, Meuse River"
  )


## -----------------------------------------------------------------------------
#| label: moran-slope-check
# slope of the Moran scatterplot
moranSlope <- coef(lm(lag_z ~ z - 1, data = as.data.frame(meuseSf)))
moranSlope

# global Moran's I
moran.test(meuseSf$log_lead, wList)$estimate[1]


## -----------------------------------------------------------------------------
#| label: sim-grid
nSide <- 8
gridPts <- expand.grid(
  x = 1:nSide,
  y = 1:nSide
)
gridPts$z <- rnorm(nrow(gridPts))

# plant the hot spot: upper-right 3x3 corner
hotIdx <- gridPts$x >= 6 & gridPts$y >= 6
gridPts$z[hotIdx] <- rnorm(sum(hotIdx), mean = 4, sd = 0.5)

gridSf <- st_as_sf(gridPts, coords = c("x", "y"))


## -----------------------------------------------------------------------------
#| label: sim-map
ggplot(gridSf) +
  geom_sf(aes(color = z), size = 6) +
  scale_color_viridis_c(name = "z") +
  labs(title = "Simulated data with planted hot spot")


## -----------------------------------------------------------------------------
#| label: sim-moran-global
knnSim <- knearneigh(gridSf, k = 4)
nbSim <- knn2nb(knnSim)
wListSim <- nb2listw(nbSim, style = "W")
moran.test(gridSf$z, wListSim)


## -----------------------------------------------------------------------------
#| label: sim-localmoran
lisaSim <- localmoran_perm(gridSf$z, listw = wListSim, nsim = 999)
head(lisaSim)


## -----------------------------------------------------------------------------
#| label: sim-attach
gridSf$Ii <- lisaSim[, "Ii"]
gridSf$pval <- lisaSim[, "Pr(z != E(Ii)) Sim"]


## -----------------------------------------------------------------------------
#| label: sim-cluster-map
zSim <- as.vector(scale(gridSf$z))
lagSim <- lag.listw(wListSim, zSim)

sig <- 0.05
gridSf$cluster <- case_when(
  zSim > 0 & lagSim > 0 & gridSf$pval < sig ~ "High-High",
  zSim < 0 & lagSim < 0 & gridSf$pval < sig ~ "Low-Low",
  zSim > 0 & lagSim < 0 & gridSf$pval < sig ~ "High-Low",
  zSim < 0 & lagSim > 0 & gridSf$pval < sig ~ "Low-High",
  TRUE ~ "Not significant"
)

ggplot(gridSf) +
  geom_sf(aes(color = cluster), size = 6) +
  scale_color_manual(
    values = c(
      "High-High"       = "#d7191c",
      "Low-Low"         = "#2c7bb6",
      "High-Low"        = "#fdae61",
      "Low-High"        = "#abd9e9",
      "Not significant" = "grey80"
    )
  ) +
  labs(
    title = "LISA cluster map: simulated data",
    color = "Cluster type"
  )


## -----------------------------------------------------------------------------
#| label: meuse-localmoran
lisa <- localmoran_perm(meuseSf$log_lead, listw = wList, nsim = 999)
head(lisa)


## -----------------------------------------------------------------------------
#| label: meuse-cluster
meuseSf$Ii <- lisa[, "Ii"]
meuseSf$pval <- lisa[, "Pr(z != E(Ii)) Sim"]

sig <- 0.01
meuseSf$cluster <- case_when(
  meuseSf$z > 0 & meuseSf$lag_z > 0 & meuseSf$pval < sig ~ "High-High",
  meuseSf$z < 0 & meuseSf$lag_z < 0 & meuseSf$pval < sig ~ "Low-Low",
  meuseSf$z > 0 & meuseSf$lag_z < 0 & meuseSf$pval < sig ~ "High-Low",
  meuseSf$z < 0 & meuseSf$lag_z > 0 & meuseSf$pval < sig ~ "Low-High",
  TRUE ~ "Not significant"
)

table(meuseSf$cluster)


## -----------------------------------------------------------------------------
#| label: meuse-lisa-map
#| message: false
tmap_mode("plot")
tm_basemap("Esri.WorldGrayCanvas") +
  tm_shape(meuseSf) +
  tm_symbols(
    fill = "cluster",
    fill.scale = tm_scale_categorical(
      values = c(
        "High-High"       = "#d7191c",
        "Low-Low"         = "#2c7bb6",
        "High-Low"        = "#fdae61",
        "Low-High"        = "#abd9e9",
        "Not significant" = "#d9d9d9"
      )
    ),
    size = 0.5,
    fill_alpha = 0.8,
    fill.legend = tm_legend(title = "LISA cluster")
  )


## -----------------------------------------------------------------------------
#| label: meuse-local-i-map
ggplot(meuseSf) +
  geom_sf(aes(color = Ii), size = 2.5) +
  scale_color_gradient2(
    low = "#2c7bb6", mid = "white", high = "#d7191c",
    midpoint = 0,
    name = expression(I[i])
  ) +
  labs(title = "Local Moran's I: log lead, Meuse River") +
  theme_minimal()


## -----------------------------------------------------------------------------
#| label: birds-load
#| eval: false
# birdsSf <- readRDS("../data/birdRichnessMexico.rds")

