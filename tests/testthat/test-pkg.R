test_that("spatial loess fit with Great-circle distance", {

  set.seed(66)
  x1 <- rnorm(100, mean=-100, sd=10)
  x2 <- rnorm(100, mean=38, sd=4)
  y <- 0.1*x1 + 1*x2 - 10 + rnorm(100, 0, 1.3)
  testdata <- data.frame(LON = x1, LAT = x2, tmax = y)

  cars.lo0 <- spaloess(tmax ~ LON + LAT, testdata, method = "model.frame", distance = "Latlong", napred = FALSE, alltree = TRUE)

  cars.lo1 <- spaloess(tmax ~ LON + LAT, testdata, distance = "Latlong", napred = FALSE, alltree = TRUE)

  cars.lo2 <- spaloess(tmax ~ LON + LAT, testdata, distance = "Latlong", napred = FALSE, control = loess.control(surf="direct"), alltree = FALSE)

  cars.lo3 <- spaloess(tmax ~ LON + LAT, testdata, family = "gaussian", distance = "Latlong", napred = FALSE, alltree = FALSE)

  cars.lo4 <- spaloess(tmax ~ LON + LAT, testdata, family = "gaussian", distance = "Latlong", napred = FALSE, control = loess.control(surf="direct"), alltree = FALSE)

  # prove it here!
  expect_true(TRUE)

})

test_that("spatial loess fit with Euclidean distance", {

  set.seed(66)
  x1 <- rnorm(100, mean=-100, sd=10)
  x2 <- rnorm(100, mean=38, sd=4)
  y <- 0.1*x1 + 1*x2 - 10 + rnorm(100, 0, 1.3)
  testdata <- data.frame(LON = x1, LAT = x2, tmax = y)

  cars.lo0 <- spaloess(tmax ~ LON + LAT, testdata, method = "model.frame", distance = "Euclid", napred = FALSE, alltree = TRUE)

  cars.lo1 <- spaloess(tmax ~ LON + LAT, testdata, distance = "Euclid", napred = FALSE, alltree = TRUE)

  cars.lo2 <- spaloess(tmax ~ LON + LAT, testdata, distance = "Euclid", napred = FALSE, control = loess.control(surf="direct"), alltree = FALSE)

  cars.lo3 <- spaloess(tmax ~ LON + LAT, testdata, family = "gaussian", distance = "Euclid", napred = FALSE, alltree = FALSE)

  cars.lo4 <- spaloess(tmax ~ LON + LAT, testdata, family = "gaussian", distance = "Euclid", napred = FALSE, control = loess.control(surf="direct"), alltree = FALSE)

  # prove it here!
  expect_true(TRUE)

})
