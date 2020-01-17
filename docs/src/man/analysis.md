# Data Analysis in Space Physics

## Spectral Analysis

### FFT

### Periodogram

This is a group of techniques to determine the periodicity of data. Julia has implementations in the [DSP](https://juliadsp.github.io/DSP.jl/stable/contents/) package. Here we introduce the usage by looking at practical examples.

### Spectrogram

Spectrogram is used a lot in wave analysis. For my purpose, I use it as an approach to visualize time dependent simulation data along a continuous line region.

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

where the average ``\langle\mathbf{B}\rangle`` is defined by
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

!!! note
    In practice, the ratio of intermediate to minimum variance should be larger than 5 to give good fit of LMN.

## ULF Wave Detection

ULF waves are MHD waves: Alfvén wave, fast wave and slow wave. One basic approach to identify waves is to check the correlation of quantity perturbations.

The phase speed of shear Alfvén wave is
```math
v_{pA} = \frac{\omega}{k} = v_A \cos{\theta}
```
where ``v_A`` is the Alfvén speed and ``\theta`` is the angle between wave vector ``\mathbf{k}`` and magnetic field ``\boldsymbol{B}``.

The perturbed quantities of Alfvén waves follow these relations:
```math
\frac{\delta \mathbf{v}}{v_A} = \pm \frac{\delta \mathbf{B}}{B_0} \\
\delta \rho = 0
```
where ``\delta \boldsymbol{v}``, ``\delta \boldsymbol{B}``, and ``\delta ρ`` are perturbed plasma velocity, magnetic fields, and plasma density, respectively, and ``B_0`` is the background magnetic magnitude.

For slow and fast waves, the phase speeds are
```math
v_{p\pm}^2 = \big(\frac{\omega}{k} \big) = \frac{1}{2}(v_s^2 + v_A^2) \pm \frac{1}{2}\Big[ (v_s^2 + v_A^2)^2 - 4v_s^2 v_A^2 \cos^2{\theta}\Big]^{1/2}
```
The "+" is for fast waves and "−" for slow waves, and ``v_S`` is the sound speed. The perturbed quantities for fast and slow waves are
```math
\delta \rho = \frac{\rho_0}{v_p}\frac{v_A^2\sin\theta}{B_0 (v_p - v_s^2/v_p)}\delta B\\
\delta \mathbf{v} = -\frac{v_A^2 \cos{\theta}}{B_0 v_p}\delta\mathbf{B} + \frac{v_A^2 \sin{\theta}\delta B}{B_0 (v_p - v_s^2/v_p)}\frac{\mathbf{k}}{k}
```
Thus generally the Alfvén wave is identified by the correlations between velocity and magnetic field perturbations, and the fast and slow waves are identified by the negative (for slow waves) or positive (for fast waves) correlations between either density and magnetic field perturbation or thermal pressure and magnetic pressure perturbation.

For the magnetosonic waves, consider using ``\delta \mathbf{E}`` and ``\delta \mathbf{B}`` for identifying speed. The slopes of the curves ``\delta E∕\delta B`` correspond to the wave propagation speed in the spacecraft frame.

About the specific names: transverse and shear Alfvén wave refer to actually the same thing.

It is possible for a static satellite to encounter first one wave and then another wave?

A tricky part in practice is how to get the average through smoothing. Note that a real satellite moves both in time and space. Usually people do moving-box-average to get an average state within a short period.

A more careful analysis is called Walen test.

However, always keep in mind that the most reliable way of identifying waves is to calculate the dispersion relation.

## Correlation Test Between Two Variables

This part takes the reference from [R](http://www.sthda.com/english/wiki/correlation-test-between-two-variables-in-r).

Correlation test is used to evaluate the association between two or more
variables.

!!! info
    If there is no relationship between the two variables, the average of ``x`` should be the same regardless of ``y`` and vice versa.

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
t = \frac{|r|}{\sqrt{1-r^2}}\sqrt{n-2}
```
where the corresponding p-value is determined using t table distribution for ``df=n-2``. If the p-value is ``< 5\%``, then the correlation between ``x`` and ``y`` is significant.
