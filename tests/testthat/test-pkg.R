test_that("spatial loess fit with Great-circle distance", {

  set.seed(66)
  x1 <- rnorm(100, mean=-100, sd=10)
  x2 <- rnorm(100, mean=38, sd=4)
  y <- 0.1*x1 + 1*x2 - 10 + rnorm(100, 0, 1.3)
  training <- data.frame(LON = x1, LAT = x2, tmax = y)

  temp.lo0 <- spaloess(tmax ~ LON + LAT, training
    , method = "model.frame"
    , distance = "Latlong"
    , napred = FALSE
    , alltree = TRUE
  )

  temp.lo1 <- spaloess(tmax ~ LON + LAT, training
    , family = "symmetric"
    , distance = "Latlong"
    , napred = FALSE
    , alltree = TRUE
  )

  temp.lo2 <- spaloess(tmax ~ LON + LAT, training
    , distance = "Latlong"
    , family = "symmetric"
    , napred = FALSE
    , control = loess.control(surf="direct")
    , alltree = FALSE
  )

  temp.lo3 <- spaloess(tmax ~ LON + LAT, training
    , family = "gaussian"
    , distance = "Latlong"
    , napred = FALSE
    , alltree = FALSE
  )

  temp.lo4 <- spaloess(tmax ~ LON + LAT, training
    , family = "gaussian"
    , distance = "Latlong"
    , napred = FALSE
    , control = loess.control(surf="direct")
    , alltree = FALSE
  )

  y[1:5] <- NA
  training <- data.frame(LON = x1, LAT = x2, tmax = y)
  temp.lo5 <- spaloess(tmax ~ LON + LAT, training
    , family = "symmetric"
    , distance = "Latlong"
    , napred = TRUE
    , control = loess.control(surf="interpolate")
    , alltree = TRUE
  )

  set.seed(1)
  x1 <- rnorm(100, mean=-100, sd=10)
  x2 <- rnorm(100, mean=38, sd=4)
  testing <- data.frame(LON = x1, LAT = x2)

  rst4 <- predloess(
    object = temp.lo4, 
    newdata = data.frame(LON = testing$LON, LAT = testing$LAT)
  )

  rst3 <- predloess(
    object = temp.lo3, 
    newdata = data.frame(LON = testing$LON, LAT = testing$LAT)
  )

  rst2 <- predloess(
    object = temp.lo2, 
    newdata = data.frame(LON = testing$LON, LAT = testing$LAT)
  )

  rst1 <- predloess(
    object = temp.lo1, 
    newdata = data.frame(LON = testing$LON, LAT = testing$LAT)
  )

  # prove it here!
  expect_true(TRUE)

})

test_that("spatial loess fit with Euclidean distance", {

  set.seed(66)
  x1 <- rnorm(100, mean=-100, sd=10)
  x2 <- rnorm(100, mean=38, sd=4)
  y <- 0.1*x1 + 1*x2 - 10 + rnorm(100, 0, 1.3)
  training <- data.frame(LON = x1, LAT = x2, tmax = y)

  temp.lo0 <- spaloess(tmax ~ LON + LAT
    , training, method = "model.frame"
    , distance = "Euclid"
    , napred = FALSE
    , alltree = TRUE
  )

  temp.lo1 <- spaloess(tmax ~ LON + LAT, training
    , distance = "Euclid"
    , napred = FALSE
    , alltree = FALSE
    , normalize = TRUE
    , control = loess.control(statistics = "exact")
  )

  temp.lo2 <- spaloess(tmax ~ LON + LAT, training
    , distance = "Euclid"
    , napred = FALSE
    , control = loess.control(surf="direct", statistics = "exact")
    , alltree = FALSE
  )

  temp.lo3 <- spaloess(tmax ~ LON + LAT, training
    , family = "gaussian"
    , distance = "Euclid"
    , normalize = TRUE
    , control = loess.control(surf="interpolate", statistics = "approximate", trace.hat = "approximate")
    , napred = FALSE
    , alltree = FALSE
  )

  temp.lo4 <- spaloess(tmax ~ LON + LAT, training
    , family = "gaussian"
    , distance = "Euclid"
    , napred = FALSE
    , normalize = TRUE
    , control = loess.control(surf="direct", statistics = "exact", trace.hat = "approximate")
    , alltree = FALSE
  )

  set.seed(1)
  x1 <- rnorm(100, mean=-100, sd=10)
  x2 <- rnorm(100, mean=38, sd=4)
  testing <- data.frame(LON = x1, LAT = x2)

  rst4 <- predloess(
    object = temp.lo4, 
    newdata = data.frame(LON = testing$LON, LAT = testing$LAT),
    se = TRUE
  )

  rst3 <- predloess(
    object = temp.lo3, 
    newdata = data.frame(LON = testing$LON, LAT = testing$LAT)
  )

  rst2 <- predloess(
    object = temp.lo2, 
    newdata = data.frame(LON = testing$LON, LAT = testing$LAT),
    se = TRUE
  )

  rst1 <- predloess(
    object = temp.lo1, 
    newdata = data.frame(LON = testing$LON, LAT = testing$LAT)
  )

  # prove it here!
  expect_true(TRUE)

})
