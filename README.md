# Spatial Locally Weighted Regression #

This package contains enhancements to the `loess` implementation that comes with base
R. 

Here are some of the added features over `loess`:
- can predict NA directly without calling `predict.loess` after `loess` to get fitted
value at NA locations
- add `distance` argument which can be different type of distance function: "Euclid",
"Latlong" for great circle distance, and "Mahal" for Mahalanobis
- add a function which can generate kd-tree from a dataset.


### Fortran code ###
- lowesf: locally weigted regression directly on original data
```
lowesf -> ehg136(does some error checking and then pass arguments to ehg127)
            |-> ehg127(calculates the fitting, called for each observations N)
```
- lowesb: kd-tree construction and fitting
v is a long vector with length lv which contains the whole information about kd-tree and some other
things. Memory is assigned by using Calloc() function in C, `loess_workspace()`.
```
lowesb -> ehg131 
           |-> ehg126(built kd-tree)
           |-> ehg124(not sure what this function is for)
           |-> ehg139(fit at vertices, vval passed into as s(0:od, nv))
                |-> ehg127(calculates the fitting, called for each vertex(nv), s(0:od) is passed into)
                |    |-> ehg106(select q-th smallest by partial sorting)
                |-> ehg137(try to compare the cutting points xi with vertex)
                |-> ehg128(interpolation function is called here based on vval2)
```	
  1. `vval` is vector with length nvmax = max(200, N) in v. It starts at v(iv(13)) in Fortran, which is 
v[iv[12]-1] in C. Length of vval is (d+1)\*nvmax, but useful length is (d+1)\*nv.
  2. `vert` is vector only has max and min of kd-tree vertices for every dimension of predictors. All kd-tree
vertices can be found starting from v(iv(11)) in Fortran, which is v[iv[10]-1] in C, and to 
v(iv(11)+nvmax) in Fortran, which is v[iv[10]-1+nvmax]. Length of vertices is nv which is iv[5].
  3. `xi` is vector of all node points from original predictors. Length of xi is nc which is iv[4]. xi
can be found in v from v(iv(12)) in Fortran, v[iv[11]-1] in C.
  4. In ehg127 function, "b" is the design matrix, b(nf, k). k is iv(29), which equal to 
(d+2)\*(d+1)/2, for example two predictors, degree is 2, then there are 6 terms in local 
regression fit. The maximum of k is 15, which means we only can have 4 predictors at most.
  5. In ehg127, for design matrix, a preliminary factorization X = QR into R and Q with Q'Q = I
followed by SVD of R allows the pseudo-inverse to be computed efficiently.
  6. ehg106, the q nearest are searched by sorting based on the distance of big circle distance.
The weights are also calculated based on big circle distance of longitude and latitude.
  7. Not sure what is vval2?

- lowese: interpolation based on kd-tree
```
lowese -> ehg133 -> ehg128(interpolation, delta is X for each newobs)
```
  1. For each vertex, the fitted value g(hat) and d derivatives of g(hat) which estimated by taking the
slopes of the locally linear or locally quadratic fit are saved and used to do the interpolation.
  2. Each cell boundary consists of four segments that meet at vertices. On each segment, function value
g(hat) are interpolated using the unique *cubic polynomial* determined by the function and derivative 
data at the vertices, this cubic polynomial should be an univariate interpolation since there is only one
dimension at edges of cells; normal derivatives are *interpolated linearly* along the segment.
  3. Finally, blending functions interpolate across the cell by using *cubic polynomial* as well. Certain 
cross derivative terms are neeeded, as described by Barnhill(1977), but we have obtained acceptable 
results by setting those cross derivative to 0.
  4. Cubic spline/Cubic interpolation/Cubic Hermite spline:
On the unit interval (0,1), given a starting point p0 at t=0 and an ending point p1 at t=1 with starting 
tangent m0 at t=0 and ending tangent m1 at t=1, the polynomial can be defined by
  5. In `ehg128`, `z` is the location for interpolation. `z` is a vector with length d. 
loop3 is finding the cell for the `z` location.
loop5 & 6 is about assigning `vval` to `g`.
`lg` is the number of vertices per cell, and `ll` is the lower left of cell, `ur` is the upper right of 
cell.
loop7 calculates the P1P2F(called the tensor product of P1 and P2). `h` is the standardized  value of a particular
edge of cells. output from loop7 is `s`.
Then each section separated by `----` is calculating the projection of values and derivatives on each edge. For
example, `gn` is the blending interpolation to calculate project values on north side of edge, here two 
derivative values used, `g1(1)` and `g0(1)`, are respective to same direction. That is why function value
are interpolated using cubic polynomial using function and derivative data at the vertices.
`gpn` is the linearly interpolation to calculate project derivative on north side of edge, here two derivative 
values used, `g0(2)` and `g0(2)`, are respective to orthogonal direction. That is why derivatives are interpolated
linearly.

### Blending Interpolation ###
In order to understand the blending, first thing we can have a look at is the bilinearly blending on the unit
square S: 0<= u,v <= 1
```
G(u,v) = [1-u u]*[G(0,v)] + [G(u,0) G(u,1)]*[1-v] - [1-u u]*[G(0,0) G(0,1)]*[1-v]
				 [G(1,v)]				    [  v]			[G(1,0) G(1,1)] [  v]
```
All information we used in bilinearly blending are projection of G(u,v) on each edge, G(0,v),G(u,1), etc.,
and the function value at four corners, G(0,0),G(1,0), etc. If projections, G(0,v), etc. are replaced by their
linear interpolants, e.g.,
```
G(0,v) = (1-v)G(0,0) + vG(0,1)
G(1,v) = (1-v)G(1,0) + vG(1,1)
```
the result `G(u,v)` is become only the third term in previous equation, which called bilinear interpolant.
Then based on this, if we replace the linear coefficients, (1-u),u,v,(1-v) by cubic Hermite spline basis
```
c Hermite basis
phi0=(1-h)**2*(1+2*h) --> P_0
phi1=h**2*(3-2*h)     --> P_1
psi0=h*(1-h)**2       --> M_0
psi1=h**2*(h-1)       --> M_1
```
and using cubic Hermite spline(interpolation on a single interval shown as following) 
```
P_t = (2t^3-3t^2+1)P_0 + (t^3-2t^2+t)M_0 + (-2t^3+3t^2)P_1 +(t^3-t^2)M_1
```
to interpolate G(u,v) based on projections, G(0,v), etc.,the interpolation now becomes:
```
G(u,v) = P1F + P2F - P1P2F
```
where
```
P1F = [phi0(u) phi1(u) psi0(u) psi1(u)]*[G(0,v) ]
										[G(1,v) ]
										[Gu(0,v)]
										[Gu(1,v)]

P2F = [G(u,0) G(u,1) Gv(u,0) Gv(u,1)]*[phi0(v)]
                                      [phi1(v)]
                                      [psi0(v)]
                                      [psi1(v)]

P1P2F = [phi0(u) phi1(u) psi0(u) psi1(u)]*B*[phi0(v)]
	                                        [phi1(v)]
      		                                [psi0(v)]
            	                            [psi1(v)]

B = [G(0,0) G(0,1)  | Gv(0,0) Gv(0,1)  ] = [Positions | v-Slopes ]
    [G(1,0) G(1,1)  | Gv(1,0) Gv(1,1)  ]   [----------|----------]
    [---------------|------------------]   [u-Slopes  | Twists   ]
    [Gu(0,0) Gu(0,1)| Gvu(0,0) Gvu(0,1)] 
	[Gu(1,0) Gu(1,1)| Gvu(1,0) Gvu(1,1)]
```
Twists here are second derivatives, they are usually assumed to be independent of the 
order of differentiation. But in our case, we assume these cross-derivatives are zero.
Cubic interpolation is used multiple times. First we used to interpolate the projection of g(u,v) 
on each edge using vertices function values and derivatives. At the same time, we linearly 
interpolated derivatives at projection of g(u,v) on each edge. Then we use four projections and four
linear interpolated derivatives to cubic interpolate the g(u,v). It should be noted that the 
direction of derivatives used in interpolation of projections on each edge are different with the
direction of derivatives first linearly interpolated then used in cubic interpolation with four
projections and four derivatives.

### R code: "kd" element of loess object ###
- kd$xi is the nodes from original data, if loess(z \~ x+y), xi can be either x-coordinate or 
y-coordinate of the cutting point.
- kd$a specifies which dimension the point in kd$xi comes from.
2 means the second dimension
1 means the first dimension
- kd$vert is the max and min vertex coordinate. For instance, d=2, vert is a vector of
(min(x), min(y), max(x), max(y)) 
- kd$vval is the fitted value for kd-tree vertices. Thera are (d+1) values for each vertex, each of
which is result from a dot product. The first value out of (d+1) is the fitted value at kd-tree
vertices, and the rest of d values are the derivatives of locally linear or locally quadratic fit,
which will be used in interpolation

