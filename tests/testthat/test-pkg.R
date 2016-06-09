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

  temp.lo6 <- spaloess(tmax ~ LON + LAT, training
    , family = "gaussian"
    , distance = "Latlong"
    , napred = FALSE
    , enp.target = 40
    , degree = 2
    , control = loess.control(surf="direct")
    , alltree = FALSE
  )

  temp.lo7 <- spaloess(tmax ~ LON + LAT, training
    , family = "gaussian"
    , span = 2
    , distance = "Latlong"
    , napred = FALSE
    , enp.target = 40
    , degree = 2
    , control = loess.control(surf="direct")
    , alltree = FALSE
  )

  expect_warning(
    temp.lo8 <- spaloess(tmax ~ LON + LAT, training
      , family = "gaussian"
      , distance = "Latlong"
      , napred = FALSE
      , degree = 2
      , span = 0.01
      , control = loess.control(surf="direct", cel=0.001)
      , alltree = FALSE
    )
  )

  expect_error(
    temp.lo8 <- spaloess(tmax ~ LON + LAT, training
      , family = "gaussian"
      , distance = "Latlong"
      , napred = FALSE
      , degree = 2
      , span = 0.001
      , control = loess.control(surf="direct", cel=0.001)
      , alltree = FALSE
    )
  )

  expect_warning(
    temp.lo8 <- spaloess(tmax ~ LON + LAT, training
      , family = "symmetric"
      , distance = "Latlong"
      , napred = FALSE
      , degree = 2
      , span = 0.1
      , control = loess.control(surf="direct", cel=0.001)
      , alltree = FALSE
    )
  )

  y[1:5] <- NA
  training <- data.frame(LON = x1, LAT = x2, tmax = y)
  temp.lo8 <- spaloess(tmax ~ LON + LAT, training
    , family = "symmetric"
    , distance = "Latlong"
    , napred = TRUE
    , control = loess.control(surf="interpolate")
    , alltree = TRUE
    , model = TRUE
  )

  rst1 <- predloess(object = temp.lo4, newdata = NULL)

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

  expect_error(
    rst <- predloess(
      object = temp.lo4,
      newdata = data.frame(x = testing$LON, y = testing$LAT)  
    )
  )

  training <- data.frame(LON = x1, LAT = as.factor(x2), tmax = y)
  expect_error(
    temp.lo1 <- spaloess(tmax ~ LON + LAT, training
      , family = "symmetric"
      , distance = "Latlong"
      , napred = FALSE
      , alltree = TRUE
    ),
    "predictors must all be numeric"
  )

  training <- data.frame(x = x1, y = x2, tmax = y)
  expect_error(
    temp.lo1 <- spaloess(tmax ~ x + y, training
      , family = "symmetric"
      , distance = "Latlong"
      , napred = FALSE
      , alltree = TRUE
    ),
    "predictors must be longitude and latitude for great circle distance"
  )

  training <- data.frame(lat = x1, lon = x2, tmax = y)
  expect_error(
    temp.lo1 <- spaloess(tmax ~ lon + lat, training
      , family = "symmetric"
      , degree = 3
      , distance = "Latlong"
      , napred = FALSE
      , alltree = TRUE
    ),
    "'degree' must be 0, 1 or 2"
  )

  expect_error(
    temp.lo1 <- spaloess(tmax ~ lon + lat, training
      , family = "symmetric"
      , degree = 2
      , control = NULL
      , distance = "Latlong"
      , napred = FALSE
      , alltree = TRUE
    ),
    "invalid 'control' argument"
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

  temp.lo11 <- spaloess(tmax ~ LON + LAT, training
    , distance = "Euclid"
    , napred = FALSE
    , alltree = FALSE
    , normalize = TRUE
  )

  temp.lo12 <- spaloess(tmax ~ LON + LAT, training
    , distance = "Euclid"
    , family = "symmetric"
    , control = loess.control(surf="direct")
    , napred = FALSE
    , alltree = FALSE
    , normalize = TRUE
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

  temp.lo22 <- spaloess(tmax ~ LON + LAT, training
    , distance = "Euclid"
    , family = "gaussian"
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

  rst22 <- predloess(
    object = temp.lo22, 
    newdata = data.frame(LON = testing$LON, LAT = testing$LAT),
    se = TRUE
  )

  rst1 <- predloess(
    object = temp.lo1, 
    newdata = data.frame(LON = testing$LON, LAT = testing$LAT)
  )

  rst11 <- predloess(
    object = temp.lo11, 
    newdata = data.frame(LON = testing$LON, LAT = testing$LAT),
    se = TRUE
  )
  rst12 <- predloess(
    object = temp.lo12, 
    newdata = data.frame(LON = testing$LON, LAT = testing$LAT),
    se = TRUE
  )
  
  newdata <- as.matrix(data.frame(LON = testing$LON, LAT = testing$LAT))
  attr(newdata, "out.attrs") <- attributes(newdata)
  rst1 <- predloess(
    object = temp.lo2, 
    newdata = newdata,
    se = TRUE
  )
  rst1 <- predloess(
    object = temp.lo2, 
    newdata = newdata,
    se = FALSE
  )

  expect_error(
    rst <- predloess(
      object = loess(tmax ~ LON + LAT, training, method = "model.frame"),
      newdata = newdata  
    ),
    "first argument must be a \"spaloess\" object"
  )

  expect_warning(
    temp.lo2 <- spaloess(tmax ~ LON + LAT, training
      , distance = "Euclid"
      , span = 0.2
      , enp.target = 40
      , napred = FALSE
      , control = loess.control(surf="direct", statistics = "exact")
      , alltree = FALSE
    )
  )






  # prove it here!
  expect_true(TRUE)

})
