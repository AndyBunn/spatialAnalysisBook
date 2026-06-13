## -----------------------------------------------------------------------------
#| label: setup
#| echo: false
#| message: false
#| warning: false
library(tidyverse)
set.seed(42)


## -----------------------------------------------------------------------------
#| label: make-data
n <- 80
x <- rnorm(n, mean = 50, sd = 10)
y <- 2 + 0.6 * x + rnorm(n, sd = 6)

ggplot(data.frame(x, y), aes(x = x, y = y)) +
  geom_point(alpha = 0.5) +
  theme_minimal()


## -----------------------------------------------------------------------------
#| label: correlation
cor(x, y)
cor(y, x)


## -----------------------------------------------------------------------------
#| label: two-regressions
# Regression of y on x
fitYx <- lm(y ~ x)
# Regression of x on y
fitXy <- lm(x ~ y)


## -----------------------------------------------------------------------------
#| label: rearrange-xy
a <- coef(fitXy)[1]
b <- coef(fitXy)[2]

# Rearranged: y = -a/b + (1/b) * x
interceptRearranged <- -a / b
slopeRearranged <- 1 / b


## -----------------------------------------------------------------------------
#| label: two-lines-plot
#| echo: false
ggplot(data.frame(x, y), aes(x = x, y = y)) +
  geom_point(alpha = 0.5) +
  geom_abline(
    intercept = coef(fitYx)[1], slope = coef(fitYx)[2],
    color = "steelblue", linewidth = 1, aes(linetype = "lm(y ~ x)")
  ) +
  geom_abline(
    intercept = interceptRearranged, slope = slopeRearranged,
    color = "firebrick", linewidth = 1, aes(linetype = "lm(x ~ y) rearranged")
  ) +
  scale_linetype_manual(
    name = "Model",
    values = c(
      "lm(y ~ x)" = "solid",
      "lm(x ~ y) rearranged" = "dashed"
    ),
    guide = guide_legend()
  ) +
  labs(
    title = "Two regressions, one scatter plot",
    subtitle = "Same data. Different lines."
  ) +
  theme_minimal()


## -----------------------------------------------------------------------------
#| label: slope-from-correlation
r <- cor(x, y)
slopeFormula <- r * (sd(y) / sd(x))
slopeLm <- coef(fitYx)[2]

round(c(from_formula = slopeFormula, from_lm = slopeLm), 6)


## -----------------------------------------------------------------------------
#| label: correlation-test
# Small correlation, large n
rSmall <- 0.11
nLarge <- 1000
tStat <- rSmall * sqrt(nLarge - 2) / sqrt(1 - rSmall^2)
pVal <- 2 * pt(-abs(tStat), df = nLarge - 2)
pVal

