# Spatial Locally Weighted Regression #

[![Build Status](https://travis-ci.org/XiaosuTong/Spaloess.svg?branch=master)](https://travis-ci.org/XiaosuTong/Spaloess)
[![codecov](https://codecov.io/gh/XiaosuTong/Spaloess/branch/master/graph/badge.svg)](https://codecov.io/gh/XiaosuTong/Spaloess)

This package contains enhancements to the `loess` implementation that comes with base
R. Most of computation in FORTRAN are kept as same in `loess` function except the
distance calculation. Modification is made in spaloess.f file.


## Installation

```r
# from github
devtools::install_github("XiaosuTong/Spaloess")
```

## Adding Feature

Here are some of the added features over `loess`:
- Can predict NA directly without calling `predict.loess` after `loess` to get fitted
value at NA locations
- Add `distance` argument which can be different type of distance function: "Euclid",
"Latlong" for great circle distance
- Add a function which can generate kd-tree from a dataset.
- Includes all locations to build the kd-tree. Argument alltree controls this.

