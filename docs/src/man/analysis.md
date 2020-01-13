# Data Analysis in Space Physics

## Spectral Analysis

### FFT

### Periodogram

### Spectrogram


## Minimum Variance Analysis

A nice introduction is given by Bengt U.Ö.Sonnerup and Maureen Scheible.
 Here is a brief summary of the idea. The implementation of MVA can be found in
 [`MVA.jl`](https://github.com/henry2004y/VisAnaJulia/blob/master/src/space/MVA.jl).

The main purpose of minimum or maximum variance analysis (MVA) is to find, from
single-spacecraft data, an estimator for the direction normal to a
one-dimensional or approximately one-dimensional current layer, wave front, or
other transition layer in a plasma.

For real transition layers observed in space there are usually more or less
pronounced deviations from the ideal 1-D model. The layer is likely to have 2-D
or 3-D internal structures which evolve in time and to have temporal
fluctuations in the orientation of its normal as well.

The minimum variance technique is designed to deal with the situation where some
or all of the non-ideal effects mentioned above, except a systematic temporal
change in the normal direction, ``\widehat{n}``, are present. As the estimate of
 ``\widehat{n}``, the method identifies that direction in space along which the
field-component set {``\mathbf{B}^{(m)}\cdot\widehat{n}``} ``(m = 1, 2, 3...M)``
 has minimum variance. In other words, ``\widehat{n}`` is determined by
 minimisation of
```math
\sigma^2 = \frac{1}{M} \sum_{m=1}^{M}| (\mathbf{B}^{(m)} - \mathbf{B})\cdot\widehat{n} |^2
```

where the average ``<\mathbf{B}>`` is defined by
```math
\langle\mathbf{B}\rangle \equiv \frac{1}{M} \sum_{m=1}^M \mathbf{B}^{(m)}
```
and where the minimisation is subject to the normalisation constraint
``|\widehat{n}|=1``. Using a Lagrange multiplier ``\lambda`` to implement this
constraint, one then seeks the solution of the set of three homogeneous linear
equations
```math
\frac{\partial}{\partial n_x}\Big( \sigma^2 - \lambda (|\widehat{n}|^2 - 1) \Big) = 0 \\
\frac{\partial}{\partial n_y}\Big( \sigma^2 - \lambda (|\widehat{n}|^2 - 1) \Big) = 0 \\
\frac{\partial}{\partial n_z}\Big( \sigma^2 - \lambda (|\widehat{n}|^2 - 1) \Big) = 0
```
where ``\sigma^2`` is given by the equation above and ``\widehat{n}`` is
represented in terms of its three components ``(n_x, n_y, n_z)`` along the
cartesian coordinate system X, Y, Z (e.g., GSE or GSM) in which the field data
``\{\mathbf{B}^{(m)} \}`` are given. When the differentiations in equations
above have been performed, the resulting set of three equations can be written
in matrix form as
```math
\sum_{\nu=1}^{3} M_{\mu\nu}^{B} n_\nu = \lambda n_\mu
```
where the subscripts ``\mu,\nu = 1,2,3`` denote cartesian components along the
X, Y, Z system and
```math
M_{\mu\nu}^{B} \equiv \langle B_\mu B_\nu\rangle - \langle B_\mu\rangle\langle B_\nu\rangle
```
is the magnetic variance matrix. It is seen from the equation that the allowed
``\lambda`` values are the eigenvalues ``\lambda_1,\lambda_2,\lambda_3``
(given here in order of decreasing magnitude) of ``M_{\mu\nu}^{B}``. Since
``M_{\mu\nu}^{B}`` is symmetric, the eigenvalues are all real and the
corresponding eigenvectors, ``x_1``, ``x_2``, and ``x_3``, are orthogonal. The
three eigenvectors represent the directions of maximum, intermediate, and
minimum variance of the field component along each vector.

## Correlation Test Between Two Variables

This part takes the reference from [R](http://www.sthda.com/english/wiki/correlation-test-between-two-variables-in-r).

Correlation test is used to evaluate the association between two or more
variables.

!!! info
If there is no relationship between the two variables (father and son heights), the average height of son should be the same regardless of the height of the fathers and vice versa.

### Methods for correlation analyses

There are different methods to perform correlation analysis:
  * Pearson correlation (r), which measures a linear dependence between two variables (x and y). It’s also known as a parametric correlation test because it depends to the distribution of the data. It can be used only when x and y are from normal distribution. The plot of y = f(x) is named the linear regression curve.
  * Kendall tau and Spearman rho, which are rank-based correlation coefficients (non-parametric).

### Correlation formula

In the formula below,
  * x and y are two vectors of length n
  * ``\bar{x}`` and ``\bar{y}`` corresponds to the means of ``x`` and ``y``, respectively.

Pearson correlation formula
```math
r = \frac{\sum (x-\bar{x})(y-\bar{y}}{\sqrt{\sum(x-\bar{x})^2\sum(y-\bar{y})^2}}
```
The p-value (significance level) of the correlation can be determined :
1. by using the correlation coefficient table for the degrees of freedom : ``df=n−2``, where ``n`` is the number of observation in ``x`` and ``y`` variables.
2 or by calculating the ``t`` value as follows:
```math
t = \frac{r}{\sqrt{1-r^2}}\sqrt{n-2}
```
where the corresponding p-value is determined using t table distribution for ``df=n-2``. If the p-value is ``< 5\%``, then the correlation between ``x`` and ``y`` is significant.
