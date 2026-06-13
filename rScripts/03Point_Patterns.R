## ----echo=FALSE, include=FALSE------------------------------------------------
set.seed(1984)


## ----message=FALSE------------------------------------------------------------
library(tidyverse)
library(spatstat)


## -----------------------------------------------------------------------------
#| fig-width: 9
#| fig-height: 6
n <- 100
patternA <- data.frame(x = rnorm(n), y = rnorm(n), id = "A")
patternB <- data.frame(x = runif(n), y = runif(n), id = "B")
patternC <- data.frame(
  expand.grid(
    x = seq(0, 1, length.out = sqrt(n)),
    y = seq(0, 1, length.out = sqrt(n))
  ),
  id = "C"
)

patternD <- data.frame(
  expand.grid(
    x = rep(seq(0, 1, length.out = n / 20), 2),
    y = rep(seq(0, 1, length.out = n / 20), 2)
  ),
  id = "D"
)
patternD$x <- jitter(patternD$x)
patternD$y <- jitter(patternD$y)

simDat <- bind_rows(patternA, patternB, patternC, patternD)

simDat <- simDat %>%
  group_by(id) %>%
  mutate(
    x = scales::rescale(x),
    y = scales::rescale(y)
  )

ggplot(simDat, aes(x = x, y = y)) +
  geom_point() +
  coord_fixed() +
  facet_wrap(~id) +
  theme_minimal()


## -----------------------------------------------------------------------------
data(japanesepines)
data(redwood)

summary(japanesepines)
summary(redwood)

oldPar <- par(no.readonly = TRUE)
par(mfrow = c(1, 2))
plot(japanesepines)
plot(redwood)
par(oldPar)


## -----------------------------------------------------------------------------
tibble(
  x = japanesepines$x * japanesepines$window$units$multiplier,
  y = japanesepines$y * japanesepines$window$units$multiplier
) %>%
  ggplot(mapping = aes(x = x, y = y)) +
  geom_point() +
  labs(
    x = paste("x (",
      japanesepines$window$units$plural,
      ")",
      sep = ""
    ),
    y = paste("y (",
      japanesepines$window$units$plural,
      ")",
      sep = ""
    )
  ) +
  coord_fixed() +
  theme_minimal()


## -----------------------------------------------------------------------------
x <- rnorm(n=1e3)
ggplot() +
  geom_histogram(mapping = aes(x = x, after_stat(density)), fill = "grey") +
  geom_density(mapping = aes(x = x)) +
  theme_minimal()


## -----------------------------------------------------------------------------
japanesepinesDensity <- density(japanesepines)
summary(japanesepinesDensity)


## -----------------------------------------------------------------------------
#| fig-width: 7
#| fig-height: 6
persp(japanesepinesDensity, theta = 30, phi = 30)


## -----------------------------------------------------------------------------
plot(japanesepinesDensity, main = NULL) # omit title
contour(japanesepinesDensity, add = TRUE, col = "white")
points(japanesepines, pch = 20, col = "white")


## -----------------------------------------------------------------------------
japanesepinesK <- envelope(japanesepines, fun = Kest, nsim = 1e3, verbose = FALSE)


## -----------------------------------------------------------------------------
plot(japanesepinesK)


## -----------------------------------------------------------------------------
redwoodK <- envelope(redwood, fun = Kest, nsim = 1e3, verbose = FALSE)
plot(redwoodK)


## -----------------------------------------------------------------------------
plot(envelope(redwood,
  fun = Kest, nsim = 100,
  verbose = FALSE, rmax = 0.5
), main = "")


## -----------------------------------------------------------------------------
redwoodL <- envelope(redwood, fun = Lest, nsim = 1e3, verbose = FALSE)
plot(redwoodL)


## -----------------------------------------------------------------------------
ggplot(redwoodK, mapping = aes(x = r, ymin = lo - pi * r^2, ymax = hi - pi * r^2)) +
  geom_ribbon(fill = "#56B4E9", alpha = 0.3) +
  geom_line(mapping = aes(y = theo - pi * r^2), col = "grey40", linetype = "dashed") +
  geom_line(mapping = aes(y = obs - pi * r^2), col = "#D55E00", linewidth = 0.8) +
  labs(y = expression(K(r) - pi ~ r^2), x = "r") +
  theme_minimal()


## -----------------------------------------------------------------------------
data(longleaf)
summary(longleaf)
plot(longleaf)


## -----------------------------------------------------------------------------
hist(longleaf$marks)


## -----------------------------------------------------------------------------
longleafL <- envelope(longleaf, fun = Lest, nsim = 1e3, verbose = FALSE)
plot(longleafL, main = "All Trees")


## -----------------------------------------------------------------------------
bigTrees <- subset(longleaf, marks > 50)
bigTreesL <- envelope(bigTrees, fun = Lest, nsim = 1e3, verbose = FALSE)
plot(bigTreesL, main = "Big Trees")


## -----------------------------------------------------------------------------
longleafAge <- cut(longleaf,
  breaks = c(0, 30, Inf),
  labels = c("Sapling", "Adult")
)
summary(longleafAge)
plot(longleafAge)


## -----------------------------------------------------------------------------
#| fig-height: 8
adults <- subset(longleafAge, marks == "Adult", drop = TRUE)
adultsL <- envelope(adults, fun = Lest, nsim = 1e3, verbose = FALSE)

saplings <- subset(longleafAge, marks == "Sapling", drop = TRUE)
saplingsL <- envelope(saplings, fun = Lest, nsim = 1e3, verbose = FALSE)


par(mfrow = c(2, 1))
plot(adultsL, main = "Adult")
plot(saplingsL, main = "Sapling")
par(oldPar)


## -----------------------------------------------------------------------------
adultsDensity <- density(adults)
plot(adultsDensity)
points(saplings, pch = 20)


## -----------------------------------------------------------------------------
saplingAdultL <- envelope(longleafAge, "Lcross",
  i = "Sapling", j = "Adult",
  verbose = FALSE
)
plot(saplingAdultL)


## ----eval=FALSE---------------------------------------------------------------
# demo(data)


## ----eval=FALSE---------------------------------------------------------------
# myPattern <- clickppp(n = 50) # click to make 50 points
# # KDE
# plot(density(myPattern))
# points(myPattern, pch = 20)
# # L
# myPatternL <- envelope(myPattern, Lest, nsim = 100)
# plot(myPatternL)


## ----eval=FALSE---------------------------------------------------------------
# patternD <- simDat %>% filter(id == "D")
# patternD <- ppp(
#   x = patternD$x, y = patternD$y,
#   xrange = c(0, 1), yrange = c(0, 1),
#   unitname = "km"
# )
# summary(patternD)
# plot(patternD)


## ----eval=FALSE---------------------------------------------------------------
# data(bei)
# summary(bei)
# plot(bei)
# 
# # First order: the intensity surface. Where are the trees dense?
# plot(density(bei))
# 
# # The bandwidth controls how smooth that surface is. The default
# # picks sigma for you (about 62 m here). Bracket it with a spiky
# # value and an over-smoothed one and watch the surface change:
# plot(density(bei, sigma = 10))
# plot(density(bei, sigma = 200))
# 
# # Second order: where does bei depart from CSR?
# plot(envelope(bei, fun = Lest, nsim = 99, verbose = FALSE))


## -----------------------------------------------------------------------------
data(sporophores)
summary(sporophores)
plot(sporophores, chars = c(16, 1, 2), cex = 0.6, leg.args = list(cex = 1.1))
points(0, 0, pch = 16, cex = 2)
text(15, 8, "Tree", cex = 0.75)


## -----------------------------------------------------------------------------
# using subset
hebPub <- subset(sporophores,
  marks %in% c("L pubescens", "Hebloma spp"),
  drop = TRUE
)
crossHebPub <- envelope(hebPub, "Lcross", verbose = FALSE)
# or directly in envelope (note i and j arguments)
crossHebPub <- envelope(sporophores, "Lcross", i = "L pubescens", j = "Hebloma spp", verbose = FALSE)

