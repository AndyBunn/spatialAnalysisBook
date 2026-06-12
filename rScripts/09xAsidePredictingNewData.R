## ----message=FALSE------------------------------------------------------------
library(tidyverse)


## -----------------------------------------------------------------------------
set.seed(42)
n <- 50
b0 <- 4
b1 <- 0.6
b2 <- 3
originalDat <- data.frame(x1 = rnorm(n), x2 = rnorm(n), epsilon = rnorm(n, sd = 2))
originalDat$y <- b0 + b1 * originalDat$x1 + b2 * originalDat$x2 + originalDat$epsilon


## -----------------------------------------------------------------------------
lm1 <- lm(y ~ x1 + x2, data = originalDat)
summary(lm1)


## -----------------------------------------------------------------------------
newDat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
head(newDat)


## -----------------------------------------------------------------------------
newDat$yhat <- predict(object = lm1, newdata = newDat)
head(newDat)

