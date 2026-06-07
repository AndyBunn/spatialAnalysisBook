## ----echo=FALSE, include=FALSE------------------------------------------------
set.seed(1984)


## ----message=FALSE------------------------------------------------------------
library(sf)
library(tidyverse)
library(spdep)
library(tmap)
library(gstat)


## -----------------------------------------------------------------------------
meuse2 <- readRDS("../data/meuse2.Rds")
meuse_sf <- st_as_sf(meuse2, coords = c("x", "y")) %>%
  st_set_crs(value = 28992)
meuse_sf$log_lead <- log(meuse_sf$lead)


## -----------------------------------------------------------------------------
k8 <- knearneigh(meuse_sf, k = 8)
nb_k8 <- knn2nb(k8)
W <- nb2listw(nb_k8, style = "W")


## -----------------------------------------------------------------------------
# standardize log lead
meuse_sf$z <- as.vector(scale(meuse_sf$log_lead))

# compute the spatial lag: weighted average of neighbors' z values
meuse_sf$lag_z <- lag.listw(W, meuse_sf$z)

ggplot(meuse_sf, aes(x = z, y = lag_z)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ x - 1, 
              se = FALSE, color = "firebrick") +
  labs(x = "Standardized log lead (z)", 
       y = "Spatial lag of z",
       title = "Moran scatterplot: log lead, Meuse River")


## -----------------------------------------------------------------------------
# slope of the Moran scatterplot
moranSlope <- coef(lm(lag_z ~ z - 1, data = as.data.frame(meuse_sf)))
moranSlope

# global Moran's I
moran.test(meuse_sf$log_lead, W)$estimate[1]


## -----------------------------------------------------------------------------
n_side <- 8
grid_pts <- expand.grid(
  x = 1:n_side,
  y = 1:n_side
)
grid_pts$z <- rnorm(nrow(grid_pts))

# plant the hot spot: upper-right 3x3 corner
hot_idx <- grid_pts$x >= 6 & grid_pts$y >= 6
grid_pts$z[hot_idx] <- rnorm(sum(hot_idx), mean = 4, sd = 0.5)

grid_sf <- st_as_sf(grid_pts, coords = c("x", "y"))


## -----------------------------------------------------------------------------
ggplot(grid_sf) +
  geom_sf(aes(color = z), size = 6) +
  scale_color_viridis_c(name = "z") +
  labs(title = "Simulated data with planted hot spot")


## -----------------------------------------------------------------------------
nb_sim <- knn2nb(knearneigh(grid_sf, k = 4))
W_sim <- nb2listw(nb_sim, style = "W")
moran.test(grid_sf$z, W_sim)


## -----------------------------------------------------------------------------
lisa_sim <- localmoran(grid_sf$z, listw = W_sim)
head(lisa_sim)


## -----------------------------------------------------------------------------
grid_sf$Ii <- lisa_sim[, "Ii"]
grid_sf$pval <- lisa_sim[, "Pr(z != E(Ii))"]


## -----------------------------------------------------------------------------
z_sim <- as.vector(scale(grid_sf$z))
lag_sim <- lag.listw(W_sim, z_sim)

sig <- 0.05
grid_sf$cluster <- case_when(
  z_sim > 0 & lag_sim > 0 & grid_sf$pval < sig ~ "High-High",
  z_sim < 0 & lag_sim < 0 & grid_sf$pval < sig ~ "Low-Low",
  z_sim > 0 & lag_sim < 0 & grid_sf$pval < sig ~ "High-Low",
  z_sim < 0 & lag_sim > 0 & grid_sf$pval < sig ~ "Low-High",
  TRUE ~ "Not significant"
)

ggplot(grid_sf) +
  geom_sf(aes(color = cluster), size = 6) +
  scale_color_manual(
    values = c(
      "High-High"      = "#d7191c",
      "Low-Low"        = "#2c7bb6",
      "High-Low"       = "#fdae61",
      "Low-High"       = "#abd9e9",
      "Not significant" = "grey80"
    )
  ) +
  labs(title = "LISA cluster map: simulated data",
       color = "Cluster type")


## -----------------------------------------------------------------------------
lisa <- localmoran(meuse_sf$log_lead, listw = W)
head(lisa)


## -----------------------------------------------------------------------------
meuse_sf$Ii <- lisa[, "Ii"]
meuse_sf$pval <- lisa[, "Pr(z != E(Ii))"]

sig <- 0.01
meuse_sf$cluster <- case_when(
  meuse_sf$z > 0 & meuse_sf$lag_z > 0 & meuse_sf$pval < sig ~ "High-High",
  meuse_sf$z < 0 & meuse_sf$lag_z < 0 & meuse_sf$pval < sig ~ "Low-Low",
  meuse_sf$z > 0 & meuse_sf$lag_z < 0 & meuse_sf$pval < sig ~ "High-Low",
  meuse_sf$z < 0 & meuse_sf$lag_z > 0 & meuse_sf$pval < sig ~ "Low-High",
  TRUE ~ "Not significant"
)

table(meuse_sf$cluster)


## -----------------------------------------------------------------------------
tmap_mode("view")
tm_shape(meuse_sf) +
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
ggplot(meuse_sf) +
  geom_sf(aes(color = Ii), size = 2.5) +
  scale_color_gradient2(
    low = "#2c7bb6", mid = "white", high = "#d7191c",
    midpoint = 0,
    name = expression(I[i])
  ) +
  labs(title = "Local Moran's I: log lead, Meuse River")


## ----eval=FALSE---------------------------------------------------------------
# birds_sf <- readRDS("../data/birdRichnessMexico.rds")


## ----echo=FALSE,eval=FALSE----------------------------------------------------
# # --- instructor key ---
# birds_sf <- readRDS("../data/birdRichnessMexico.rds")
# 
# nb_birds <- knn2nb(knearneigh(birds_sf, k = 8))
# W_birds <- nb2listw(nb_birds, style = "W")
# 
# # Moran scatterplot
# birds_sf$z <- as.vector(scale(birds_sf$nSpecies))
# birds_sf$lag_z <- lag.listw(W_birds, birds_sf$z)
# ggplot(birds_sf, aes(x = z, y = lag_z)) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   geom_point(alpha = 0.5) +
#   geom_smooth(method = "lm", formula = y ~ x - 1, se = FALSE, color = "firebrick") +
#   labs(x = "z (nSpecies)", y = "Spatial lag of z")
# 
# # LISA
# lisa_birds <- localmoran(birds_sf$nSpecies, listw = W_birds)
# birds_sf$Ii <- lisa_birds[, "Ii"]
# birds_sf$pval <- lisa_birds[, "Pr(z != E(Ii))"]
# 
# sig <- 0.01
# birds_sf$cluster <- case_when(
#   birds_sf$z > 0 & birds_sf$lag_z > 0 & birds_sf$pval < sig ~ "High-High",
#   birds_sf$z < 0 & birds_sf$lag_z < 0 & birds_sf$pval < sig ~ "Low-Low",
#   birds_sf$z > 0 & birds_sf$lag_z < 0 & birds_sf$pval < sig ~ "High-Low",
#   birds_sf$z < 0 & birds_sf$lag_z > 0 & birds_sf$pval < sig ~ "Low-High",
#   TRUE ~ "Not significant"
# )
# 
# tmap_mode("view")
# tm_shape(birds_sf) +
#   tm_symbols(
#     fill = "cluster",
#     fill.scale = tm_scale_categorical(
#       values = c(
#         "High-High"       = "#d7191c",
#         "Low-Low"         = "#2c7bb6",
#         "High-Low"        = "#fdae61",
#         "Low-High"        = "#abd9e9",
#         "Not significant" = "#d9d9d9"
#       )
#     ),
#     size = 0.3,
#     fill_alpha = 0.8
#   )
# # HH in southern Mexico (Veracruz, Oaxaca, Chiapas) -- expected from tropical biodiversity
# # LL in Baja, Sonoran desert -- also expected
# # Some HL outliers worth examining

