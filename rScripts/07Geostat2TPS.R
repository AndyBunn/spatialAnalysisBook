## ----echo=FALSE, message=FALSE, warning=FALSE, results='hide'-----------------
set.seed(184)


## -----------------------------------------------------------------------------
#| warning: false
#| message: false
library(fields)
library(tidyverse)
library(gstat)
library(sf)
library(terra)
library(tidyterra)


## -----------------------------------------------------------------------------
#| warning: false
#| message: false
meuse2 <- readRDS("../data/meuse2.Rds")
meuse2$logLead <- log(meuse2$lead)
meuseSf <- st_as_sf(meuse2, coords = c("x", "y")) %>%
  st_set_crs(value = 28992)

meuse.grid2 <- readRDS("../data/meuse.grid2.Rds")


## -----------------------------------------------------------------------------
ggplot(data = meuseSf) +
  geom_sf(aes(fill = logLead),
    size = 4,
    shape = 21, color = "white", alpha = 0.8
  ) +
  scale_fill_continuous(type = "viridis", name = "log(ppm)") +
  labs(title = "Lead concentrations")


## -----------------------------------------------------------------------------
logLeadTPSmodel <- Tps(x = meuse2[, 1:2], Y = meuse2$logLead)
logLeadTPSmodel


## -----------------------------------------------------------------------------
# Predict the model over all the coordinates in meuse.grid2
logLeadPreds <- c(predict(object = logLeadTPSmodel, x = meuse.grid2[, 1:2]))
# Store in a data.frame with the x,y coordinates
logLeadTPS <- data.frame(
  x = meuse.grid2[, 1],
  y = meuse.grid2[, 2],
  logLead = logLeadPreds
)
# And into SpatRaster
logLeadTPSRast <- rast(logLeadTPS, crs = crs(meuseSf))

# Plot
ggplot() +
  geom_spatraster(data = logLeadTPSRast, mapping = aes(fill = logLead), alpha = 0.8) +
  scale_fill_continuous(type = "viridis", name = "log(ppm)", na.value = "transparent") +
  labs(title = "Lead concentrations", subtitle = "TPS") +
  theme_minimal()


## -----------------------------------------------------------------------------
obs <- meuse2$logLead
preds <- extract(logLeadTPSRast, meuseSf) %>% pull(logLead)
rsq <- cor(obs, preds)^2
rmse <- sqrt(mean((preds - obs)^2))
c(rsq = rsq, rmse = rmse)


## -----------------------------------------------------------------------------
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
n <- nrow(meuse2)
rows4test <- sample(x = 1:n, size = n * 0.2)
meuseTest <- meuse2[rows4test, ]
meuseTrain <- meuse2[-rows4test, ]

# Note meuseTrain here
logLeadTPSmodel <- Tps(x = meuseTrain[, 1:2], Y = meuseTrain$logLead)
logLeadTPSmodel

# Predict the model over all the coordinates in meuse.grid2
logLeadPreds <- c(predict(object = logLeadTPSmodel, x = meuse.grid2[, 1:2]))
# Store in a data.frame with the x,y coordinates
logLeadTPS <- data.frame(
  x = meuse.grid2[, 1],
  y = meuse.grid2[, 2],
  logLead = logLeadPreds
)
# And into SpatRaster
logLeadTPSRast <- rast(logLeadTPS, crs = crs(meuseSf))

# Look at skill on withheld data. Note meuseTest:
obs <- meuseTest$logLead
preds <- extract(logLeadTPSRast, meuseTest[, 1:2]) %>% pull(logLead)
rsq <- cor(obs, preds)^2
rmse <- sqrt(mean((preds - obs)^2))
c(rsq = rsq, rmse = rmse)


## -----------------------------------------------------------------------------
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
rmseNULL <- sqrt(mean((mean(meuseTrain$logLead) - obs)^2))
rmseNULL


## -----------------------------------------------------------------------------
1 - (rmse / rmseNULL)

