
# Gaussian Basis Functions

This page provides an overview of the centered Gaussian basis functions used for 2‑ and 3‑body quantum calculations in the package. It covers standard real-range Gaussians, an extension to complex-ranged Gaussians, and the methods used for higher angular momenta in three-body systems. A more extensive introduction, however limited to 3D systems can be found in [hiyama2003](@cite).

## Two-Body

The module [GEM2B](@ref GEM2B) provides a solver for two-body systems described by the Schrödinger equation
```math
\left[-\frac{\hbar^2}{2\mu} \nabla_{\vec{r}}^2 + V(|\vec{r}|) \right]\phi(r) = E \phi(r).
```
The Laplacian should be interpreted according to the dimensionality (1D, 2D, 3D).

### Real-range (centered) Gaussian basis functions

The Gaussian expansion method is based on expanding an unknown state into a set
```math
|\phi\rangle = \sum_\alpha c_\alpha |\phi_\alpha\rangle
```
of Gaussian basis functions `` |\phi_\alpha\rangle `` . They are originally defined in 3D as
```math
\phi_\alpha^{(3D)}(r) = N_{n,l}\, r^l e^{-\nu_n r^2} Y_{l,m}(\theta,\phi)
```
with the multi-index `` \alpha = \{n,l\}  ``. This definition can be extended to 1D and 2D systems by considering (the normalization depends on the dimensionality)
```math
\phi_\alpha^{(1D)}(r) = N_{n,l}\, r^l e^{-\nu_n r^2} \\
\phi_\alpha^{(2D)}(r) = N_{n,l}\, r^l e^{-\nu_n r^2} e^{im\phi}.
```

In 3D, the index `` l `` denotes the angular momentum. In 1D, we can associate even and odd values of it to the parity or symmetry of a two-body state. The index `` n `` is associated with the Gaussian ranges `` \nu_n `` which are chosen in a geometric progression
```math
\nu_n = 1/r_n^2, \qquad r_n = r_1 a^{n-1}, \qquad a = \left(\frac{r_{n_{max}}}{r_1} \right)^{1/(n_{max}-1)}
```
defined by the parameters `nmax`, `r1`, `rnmax`. This choice allows for a large function space, while keeping a low number of numerical parameters. The series of ranges is created in the function `buildnu`

Using Gaussian basis functions allows to treat many steps fully analytically, and is especially useful for 3-body calculations. However, they are not orthogonal which implies that finding the solution to a few-body problem requires solving a generalized eigenvalue problem, and not a standard one. For this, we provide the function `eigen2step`, which is similar to `eigen(A,B)`, but it removes eigenvalues of `B` below a certain threshold, to cure a possibly ill-posed problem. Such a situation can arise if a set of too many basis functions with similar ranges is used.

### Complex-ranged Gaussians

For oscillatory states, e.g. highly excited states or metastable resonant states, it can be difficult to capture their details with standard Gaussians basis functions. A simple extension suggested in [hiyama2003](@cite) are complex-ranged Gaussians which are obtained from the standard real-ranged ones by the transformation
```math
\nu \to \nu (1 \pm i\,\omega).
```
This allows reusing the same code structure as for real-ranged Gaussians. Note that to ensure a real wave function, always the pair of both complex shifts is used, and hence the number of basis functions is twice the one defined in `num_params`. Complex-ranged Gaussians can be used with the optional keyword argument `cr_bool = 1` when calling the solver.


## Three-Body

The modules [GEM3B1D](@ref GEM3B1D), and [ISGL](@ref ISGL) provide a solver for three-body systems described by the Schrödinger equation
```math
\left[ -\frac{\hbar^2}{2\mu_{ij}} \nabla_{\vec{r}_{k}}^2 - \frac{\hbar^2}{2\mu_{k}} \nabla_{\vec{R}_{k}}^2 + V_{12}(\vec{r}_{12}) + V_{23}(\vec{r}_{23}) + V_{31}(\vec{r}_{31}) \right] \Psi(\vec{r}_{k},\vec{R}_{k}) = E \Psi(\vec{r}_{k},\vec{R}_{k})
```
with the reduced masses ``\mu_{ij} = (m_i m_j)/(m_i + m_j)) ``, and `` \mu_k = (m_k (m_i+m_j))/(m_i+m_j+m_k)``. The Laplacians should be interpreted accordingly to the dimensionality (1D, 3D). For demonstrational purposes, this equation is written in one specific Jacobi set, defined by assigning ``\{i,j,k\} `` to the particles ``\{1,2,3\}``. The solvers employ a combination of up to three Jacobi sets at once. The required number is determined automatically.

### Jacobi Coordinates and Faddeev Components


![Jacobi coordinate sets](assets/JacobiCoordinates.svg)

**Figure 1**: Three sets of Jacobi coordinates for a three-body system.

For three-body systems we employ Jacobi coordinates. This allows to describe the full system inits center-of-mass frame by two relative coordinates. However, there are three equivalent sets of these coordinates, related to the three different partitions of three particles into a pair of two, and a single one, see Fig. 1. In few-body physics it is therefore common to decompose any given three-body state into a sum
```math
|\Psi\rangle = |\Psi^{(1)}\rangle + |\Psi^{(2)}\rangle + |\Psi^{(3)}\rangle
```
of Faddeev components `` |\Psi^{(i)}\rangle  `` (sometimes called rearrangement channels), each described in a different Jacobi set `` \{\vec{r}_i,\vec{R}_i\} ``. In case some particles do not interact, or two or more are identical, the number of Faddeev components can be reduced. The code does this automatically.

Each component is then expanded into a set
```math
\Psi^{(i)}(\vec{r}_i,\vec{R}_i) = \sum_{\alpha=1}^{\alpha_{max}} A_\alpha \psi_\alpha^{(i)}(\vec{r}_i,\vec{R}_i)
```
of basis functions `` \psi_\alpha^{(i)} `` which itself are composed of products of two functions
```math
\psi_\alpha^{(i)}(\vec{r}_i,\vec{R}_i) = \phi_\alpha(\vec{r}_i) \Phi_\alpha(\vec{R}_i).
```
These functions are each defined as in the [two-body case](#Two-body)
```math
\phi_\alpha(\vec{r}) = N_{l,m} r^l e^{-\nu_n r^2} Y_{l,m}(\hat{r})
```
```math
\Phi_\alpha(\vec{R}) = N_{L,M} R^L e^{-\lambda_N R^2} Y_{L,M}(\hat{R})
```
where now `` \alpha = \{n,l,N,L\}  ``. For 2D and 1D systems, the spherical harmonics are simply replaced by `` e^{i m \varphi} `` or 1, respectively.

The basis is automatically constructed based on the inputs `nmax,r1,rnmax,Nmax,R1,RNmax ` from `gem_params` in `num_params`.

The decomposition into different Faddeev components requires that we have to compute matrix elements of the form
```math
\langle \Psi^{(a)} | \hat{O}^{(c)} |\Psi^{(b)}\rangle
```
with basis functions and operators in possibly different Jacobi sets `` a,b,c ``. Owing to the Gaussian form, this can be done readily for s-wave states, since the Gaussian form is preserved when transforming from one set to another.

### Higher angular momenta via infinitesimally shifted Gaussian lobe (ISGL) functions
If higher angular momenta are required, the transformation between different Jacobi sets quickly becomes laborious due to the presence of multiple spherical harmonics. To avoid this, the code makes use of so-called infinitesimally shifted Gaussian basis functions, to express the spherical harmonics in terms of several Gaussians with infinitesimal shift. This shift is treated analytically before computing the matrix elements. More details on this type of basis functions can be found in the review article [hiyama2003](@cite). We note that this treatment is only required in 3D. In 1D, matrix elements of different Jacobi sets can be computed more easily. Hence, we have two separate modules, [GEM3B1D](#GEM3B1D) for 1D, and [ISGL](#ISGL) for 3D. Three-body systems in 2D are currently not supported.


## Page references

```@bibliography
Pages = ["BasisFunctions.md"]
Canonical = false
```

See also the [full bibliography](@ref References) for further references cited throughout this documentation.