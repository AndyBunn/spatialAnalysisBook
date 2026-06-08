## ----echo=FALSE, message=FALSE, warning=FALSE, results='hide'-----------------
set.seed(733)


## ----message=FALSE, warning=FALSE---------------------------------------------
library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(gstat)
library(cowplot)


## -----------------------------------------------------------------------------
meuse2 <- readRDS("../data/meuse2.Rds")
meuse.grid2 <- readRDS("../data/meuse.grid2.Rds")


## -----------------------------------------------------------------------------
meusePoints_sf <- st_as_sf(meuse2, coords = c("x","y"), crs=28992)
meuseGrid_sf   <- st_as_sf(meuse.grid2, coords = c("x","y"), crs=28992)

covars_sf <- meuseGrid_sf %>% select(river_dist_m)
meusePoints_sf <- st_join(meusePoints_sf, covars_sf, join = st_nearest_feature)


## -----------------------------------------------------------------------------
meusePoints_sf <- meusePoints_sf %>%
  mutate(log_dist = log(river_dist_m)) %>%
  drop_na(om)
meuseGrid_sf <- meuseGrid_sf %>%
  mutate(log_dist = log(river_dist_m))


## -----------------------------------------------------------------------------
om_lm <- lm(om ~ log_dist, data = meusePoints_sf)
summary(om_lm)


## -----------------------------------------------------------------------------
meuseGrid_sf$omHat <- predict(om_lm, newdata=meuseGrid_sf)

# same sf_2_rast function from the IDW notes -- write it once, use it everywhere
sf_2_rast <- function(sfObject, variable2get = 1){
  tmp <- sfObject[,variable2get] %>% st_drop_geometry()
  dfObject <- data.frame(st_coordinates(sfObject), z=tmp)
  rastObject <- rast(dfObject, crs=crs(sfObject))
  return(rastObject)
}

meuseGrid_rast <- sf_2_rast(meuseGrid_sf, variable2get = "omHat")
ggplot() + geom_spatraster(data = meuseGrid_rast,
                           mapping = aes(fill=omHat)) +
  scale_fill_terrain_c() +
  labs(fill = "omHat")



## ----echo=FALSE---------------------------------------------------------------
w <- spdep::knn2nb(spdep::knearneigh(meusePoints_sf, k=4))
foo <- spdep::moran.test(residuals(om_lm), spdep::nb2listw(w))


## -----------------------------------------------------------------------------
omVar <- variogram(om ~ 1, meusePoints_sf)
plot(omVar, pch=20, cex=1.5, col="black",
     ylab=expression("Semivariance ("*gamma*")"),
     xlab="Distance (m)", main = "% Soil Organic Matter")


## -----------------------------------------------------------------------------
omGstat <- gstat(id = "omModel", formula = om ~ log_dist,
                 data = meusePoints_sf)
omGstat_obsVariogram <- variogram(omGstat)
plot(omGstat_obsVariogram, pch=20, cex=1.5, col="black",
     ylab=expression("Semivariance ("*gamma*")"),
     xlab="Distance (m)", main = "Model Residuals")


## -----------------------------------------------------------------------------
omGau_variogramModel <- vgm(psill = 2, model = "Gau", range = 500, nugget = 4)
omGau_fittedVariogram <- fit.variogram(object = omGstat_obsVariogram, 
                                       model = omGau_variogramModel)
plot(omGstat_obsVariogram, omGau_fittedVariogram, pch=20, cex=1.5, col="black",
     ylab=expression("Semivariance ("*gamma*")"),
     xlab="Distance (m)", main = "Model Residuals")


## -----------------------------------------------------------------------------
# Update the gstat object with the variogram:
omGstat_w_variogram <- gstat(omGstat, id="omModel", model=omGau_fittedVariogram)
# And predict
omHat_sf <- predict(omGstat_w_variogram, newdata = meuseGrid_sf)
omHat_sf


## -----------------------------------------------------------------------------
#| fig-width: 10
#| fig-height: 9
omHat_rast <- sf_2_rast(omHat_sf, variable2get = "omModel.pred")

# kriged residual surface = RK prediction minus OLS trend
omHat_sf$krigedResid <- omHat_sf$omModel.pred - meuseGrid_sf$omHat
omResid_rast <- sf_2_rast(omHat_sf, variable2get = "krigedResid")

p_ols <- ggplot() +
  geom_spatraster(data = meuseGrid_rast, mapping = aes(fill=omHat)) +
  scale_fill_terrain_c(limits=c(0,18)) +
  labs(title="OLS trend", fill="om (%)") +
  theme(axis.text=element_blank(), axis.ticks=element_blank())

p_resid <- ggplot() +
  geom_spatraster(data = omResid_rast, mapping = aes(fill=krigedResid)) +
  scale_fill_gradient2(low="blue", mid="white", high="red", 
                     midpoint=0, na.value="transparent") +
  labs(title="Kriged residuals", fill="resid") +
  theme(axis.text=element_blank(), axis.ticks=element_blank())

p_rk <- ggplot() +
  geom_spatraster(data = omHat_rast, mapping = aes(fill=omModel.pred)) +
  scale_fill_terrain_c(limits=c(0,18)) +
  labs(fill="om (%)") +
  theme(axis.text=element_blank(), axis.ticks=element_blank())


eq <- ggdraw() + draw_label("OLS trend  +  Kriged residuals  =  RK prediction", 
                             fontface="italic", size=11)

top_row <- plot_grid(p_ols, p_resid, ncol=2)
plot_grid(top_row, eq, p_rk, nrow=3, rel_heights=c(1, 0.1, 1))



## -----------------------------------------------------------------------------
omHatVar_rast <- sf_2_rast(omHat_sf, variable2get = "omModel.var")
omHatVar_rast <- omHatVar_rast %>% mutate(omModel.var.sqrt = sqrt(omModel.var))
ggplot() + geom_spatraster(data = omHatVar_rast,
                           mapping = aes(fill=omModel.var.sqrt)) +
  scale_fill_terrain_c() +
  labs(fill = "Organic Matter SD (%)")

