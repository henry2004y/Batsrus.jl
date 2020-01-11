# Data Analysis in Space Physics

## Spectral Analysis

### FFT

### Periodogram

### Spectrogram


## Minimum Variance Analysis

A nice introduction is given by BENGT U. ``\''{O}``. SONNERUP AND MAUREEN
SCHEIBLE. Here is a brief summary of the idea.

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
<\mathbf{B}> \equiv \frac{1}{M} \sum_{m=1}^M \mathbf{B}^{(m)}
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
M_{\mu\nu}^{B} \equiv <B_\mu B_\nu> - <B_\mu><B_\nu>
```
is the magnetic variance matrix. It is seen from the equation that the allowed
``\lambda`` values are the eigenvalues ``\lambda_1,\lambda_2,\lambda_3``
(given here in order of decreasing magnitude) of ``M_{\mu\nu}^{B}``. Since
``M_{\mu\nu}^{B}`` is symmetric, the eigenvalues are all real and the
corresponding eigenvectors, ``x_1``, ``x_2``, and ``x_3``, are orthogonal. The
three eigenvectors represent the directions of maximum, intermediate, and
minimum variance of the field component along each vector.

The implementation of MVA can be found in [`MVA.jl`](src/space/MVA.jl).
