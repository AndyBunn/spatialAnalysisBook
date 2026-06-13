## ----echo=FALSE---------------------------------------------------------------
#| label: setup
set.seed(2613)


## ----message=FALSE------------------------------------------------------------
#| label: packages
library(tidyverse)
library(sf)
library(gstat)
library(PNWColors)


## ----echo=FALSE---------------------------------------------------------------
#| label: palette
anem7 <- pnw_palette(name = "Anemone", n = 7, type = "discrete")


## -----------------------------------------------------------------------------
#| label: sim-data
n <- 50
x <- runif(n, min = -4, max = 4)
eps <- rnorm(n, sd = 15)
y <- 30 + 2 * x^3 - 5 * x + eps
dat <- data.frame(x = x, y = y)


## -----------------------------------------------------------------------------
#| label: sim-plot
ggplot(dat, aes(x = x, y = y)) +
  geom_point(size = 2.5, alpha = 0.8) +
  theme_minimal()


## -----------------------------------------------------------------------------
#| label: overfit-curve
testOvIdx <- sample(1:nrow(dat), size = round(nrow(dat) * 0.33))
trainOv <- dat[-testOvIdx, ]
testOv <- dat[testOvIdx, ]

degrees <- 1:10
rsqIn <- numeric(length(degrees))
rsqOut <- numeric(length(degrees))

for (d in degrees) {
  fitD <- lm(y ~ poly(x, d), data = trainOv)
  rsqIn[d] <- cor(trainOv$y, predict(fitD))^2
  rsqOut[d] <- cor(testOv$y, predict(fitD, newdata = testOv))^2
}

rsqByDegree <- data.frame(
  degree = rep(degrees, 2),
  rsq    = c(rsqIn, rsqOut),
  set    = rep(c("Training", "Test"), each = length(degrees))
)

ggplot(rsqByDegree, aes(x = degree, y = rsq, color = set)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = degrees) +
  scale_color_manual(values = c(anem7[7], anem7[1])) +
  labs(x = "Polynomial degree", y = bquote(R^2), color = NULL) +
  theme_minimal()


## -----------------------------------------------------------------------------
#| label: holdout-split
testIdx <- sample(1:nrow(dat), size = round(nrow(dat) * 0.25))
dat$split <- "train"
dat$split[testIdx] <- "test"

ggplot(dat, aes(x = x, y = y, color = split)) +
  geom_point(size = 3) +
  scale_color_manual(values = c(anem7[1], anem7[7])) +
  theme_minimal()


## -----------------------------------------------------------------------------
#| label: holdout-metrics
train <- dat[dat$split == "train", ]
test <- dat[dat$split == "test", ]

fit <- lm(y ~ poly(x, 3), data = train)
test$yhat <- predict(fit, newdata = test)

rsq <- cor(test$y, test$yhat)^2
rmse <- sqrt(mean((test$yhat - test$y)^2))
mae <- mean(abs(test$yhat - test$y))

c("rsq" = rsq, "rmse" = rmse, "mae" = mae)


## -----------------------------------------------------------------------------
#| label: holdout-obs-pred
ggplot(test, aes(x = y, y = yhat)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(size = 3, color = anem7[7]) +
  coord_fixed(
    xlim = range(test$y, test$yhat),
    ylim = range(test$y, test$yhat)
  ) +
  labs(x = "Observed", y = "Predicted") +
  theme_minimal()


## -----------------------------------------------------------------------------
#| label: null-model
rmseNull <- sqrt(mean((mean(train$y) - test$y)^2))
maeNull <- mean(abs(mean(train$y) - test$y))

rmseNull
1 - (rmse / rmseNull)


## -----------------------------------------------------------------------------
#| label: kfold
k <- 10
dat2 <- dat[sample(nrow(dat)), c("x", "y")] # shuffle, drop the split column
dat2$fold <- cut(seq(1, nrow(dat2)), breaks = k, labels = FALSE)

results <- data.frame(fold = 1:k, rsq = NA, rmse = NA, mae = NA)

for (i in 1:k) {
  trainI <- dat2[dat2$fold != i, ]
  testI <- dat2[dat2$fold == i, ]

  fitI <- lm(y ~ poly(x, 3), data = trainI)
  yhatI <- predict(fitI, newdata = testI)

  results$rsq[i] <- cor(testI$y, yhatI)^2
  results$rmse[i] <- sqrt(mean((yhatI - testI$y)^2))
  results$mae[i] <- mean(abs(yhatI - testI$y))
}

results


## -----------------------------------------------------------------------------
#| label: kfold-means
colMeans(results[, c("rsq", "rmse", "mae")])


## -----------------------------------------------------------------------------
#| label: kfold-hist
results %>%
  pivot_longer(
    cols = c(rsq, rmse, mae),
    names_to = "metric", values_to = "value"
  ) %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 6, fill = anem7[5], color = anem7[7]) +
  labs(x = NULL, y = NULL) +
  facet_wrap(~metric, scales = "free_x") +
  theme_minimal()


## ----message=FALSE------------------------------------------------------------
#| label: meuse-holdout
meuse2 <- readRDS("../data/meuse2.Rds")
meuseSf <- st_as_sf(meuse2, coords = c("x", "y"), crs = 28992)

testMIdx <- sample(1:nrow(meuseSf), size = round(nrow(meuseSf) * 0.25))
meuseSf$split <- "train"
meuseSf$split[testMIdx] <- "test"

ggplot(meuseSf) +
  geom_sf(aes(color = split, shape = split), size = 3) +
  scale_color_manual(values = c(anem7[1], anem7[7])) +
  theme_minimal() +
  labs(
    title = "Random 75/25 holdout",
    subtitle = "Test points scattered among training points"
  )

