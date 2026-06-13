## -----------------------------------------------------------------------------
#| label: setup
#| echo: false
set.seed(1984)


## -----------------------------------------------------------------------------
#| label: packages
#| message: false
library(tidyverse)
library(spatstat)


## -----------------------------------------------------------------------------
#| label: point-pattern
n <- 12
area <- 100 # 10 x 10 = 100 square units
x <- runif(n = n, min = 0, max = 10)
y <- runif(n = n, min = 0, max = 10)

ggplot() +
  geom_point(aes(x = x, y = y), size = 3) +
  coord_fixed() +
  theme_minimal()


## -----------------------------------------------------------------------------
#| label: dist-matrix
Dmat <- as.matrix(dist(cbind(x, y)))
Dmat


## -----------------------------------------------------------------------------
#| label: dist-vector
Dvec <- c(Dmat[upper.tri(Dmat)], Dmat[lower.tri(Dmat)])
Dvec


## -----------------------------------------------------------------------------
#| label: k-by-hand
r <- seq(0, 2.5, by = 0.1)

Kr <- numeric(length(r))
for (i in seq_along(r)) {
  Kr[i] <- (area / (n * (n - 1))) * sum(Dvec <= r[i])
}


## -----------------------------------------------------------------------------
#| label: k-curve-plot
kCurves <- tibble(
  r = r,
  Observed = Kr,
  `Theoretical (CSR)` = pi * r^2
) %>%
  pivot_longer(-r, names_to = "series", values_to = "K")

ggplot(kCurves, aes(x = r, y = K, color = series, linetype = series)) +
  geom_line(linewidth = 1) +
  scale_color_manual(
    values = c("Observed" = "#D55E00", "Theoretical (CSR)" = "grey40")
  ) +
  scale_linetype_manual(
    values = c("Observed" = "solid", "Theoretical (CSR)" = "dashed")
  ) +
  labs(y = "K(r)", x = "r", color = NULL, linetype = NULL) +
  theme_minimal() +
  theme(legend.position = "top")


## -----------------------------------------------------------------------------
#| label: kest-check
xy <- as.ppp(cbind(x, y), W = c(0, 10, 0, 10))
xyK <- Kest(xy, r = r, correction = "none")


## -----------------------------------------------------------------------------
#| label: k-overlay-plot
ggplot() +
  geom_line(
    data = kCurves, aes(x = r, y = K, color = series, linetype = series),
    linewidth = 1
  ) +
  geom_point(
    data = xyK, aes(x = r, y = un, shape = "Kest()"),
    size = 1.8
  ) +
  scale_color_manual(
    values = c("Observed" = "#D55E00", "Theoretical (CSR)" = "grey40")
  ) +
  scale_linetype_manual(
    values = c("Observed" = "solid", "Theoretical (CSR)" = "dashed")
  ) +
  scale_shape_manual(values = c("Kest()" = 1)) +
  labs(y = "K(r)", x = "r", color = NULL, linetype = NULL, shape = NULL) +
  theme_minimal() +
  theme(legend.position = "top")

