## ----echo=FALSE---------------------------------------------------------------
set.seed(1984)


## ----message=FALSE------------------------------------------------------------
library(tidyverse)
library(spatstat)


## -----------------------------------------------------------------------------
n <- 12
area <- 100  # 10 x 10 = 100 square units
x <- runif(n = n, min = 0, max = 10)
y <- runif(n = n, min = 0, max = 10)

ggplot() +
  geom_point(aes(x = x, y = y), size = 3) +
  coord_fixed()


## -----------------------------------------------------------------------------
Dmat <- as.matrix(dist(cbind(x, y)))
Dmat


## -----------------------------------------------------------------------------
Dvec <- c(Dmat[upper.tri(Dmat)], Dmat[lower.tri(Dmat)])
Dvec


## -----------------------------------------------------------------------------
r <- seq(0, 2.5, by = 0.1)

Kr <- numeric(length(r))
for (i in seq_along(r)) {
  Kr[i] <- (area / (n * (n - 1))) * sum(Dvec <= r[i])
}


## -----------------------------------------------------------------------------
ggplot() +
  geom_line(aes(x = r, y = pi * r^2), color = "red", linetype = "dashed") +
  geom_line(aes(x = r, y = Kr), color = "blue", linewidth = 1) +
  labs(y = "K(r)", x = "r",
       caption = "Blue = observed K. Red dashed = theoretical K under CSR.")


## -----------------------------------------------------------------------------
xy <- as.ppp(cbind(x, y), W = c(0, 10, 0, 10))
xyK <- Kest(xy, r = r, correction = "none")


## -----------------------------------------------------------------------------
ggplot() +
  geom_line(aes(x = r, y = pi * r^2), color = "red", linetype = "dashed") +
  geom_line(aes(x = r, y = Kr), color = "blue", linewidth = 1) +
  geom_line(data = xyK, aes(x = r, y = un),
            linetype = "dashed", color = "white") +
  labs(y = "K(r)", x = "r")

