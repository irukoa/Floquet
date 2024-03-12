# Floquet
This is a modern Fortran library companion to Ref. [[1]](#ref1), which implements the computational scheme described in that reference to calculate the quasienergy spectrum of real materials and tight-binding models subjected to time-periodic driving fields.

# Working principle

Suppose that a crystal (real material or tight-binding model) is subjected to a time periodic driving field $E^j(t)$ with period $T$. Provided that the crystal is described by a Hamiltonian $\hat{H}_0$ and by a Berry connection $\hat{\xi}^j$, the interacting dynamics in the dipole approximation are described by the Hamiltonian

$$
\hat{H}(t) = \hat{H}_0 - q \sum_j E^j(t)\cdot r^j,
$$

where $q$ is the particle charge and $r^j$ is the position operator of the system. The intricacies related to this operator in extended systems are discussed in Ref. [[1]](#ref1).

The time evolution operator is defined as

$$
\hat{U}(t, t_0) = \text{exp}\left[\frac{-i}{\hbar}\int_{t_0}^{t_0 + t}\hat{H}(\tau)d\tau\right],
$$

and in practice implemented as

$$
\hat{U}(t, t_0) = \prod_{j=1}^N\text{exp}\left[-\frac{i}{\hbar} \delta t \hat{H}(t_j)  \right]+ \mathcal{O}(\delta t^2),
$$

where

$$
t_j = t_0 + T\frac{j-1}{N_t-1} = t_0 + \delta t (j-1)
$$

The time evolution operator defines the propagation of quantum mechanical states in time. The relevant time-averaged physical properties contained in this operator can be analysed by considering the *effective Floquet Hamiltonian*,

$$
\hat{H}_F[t_0] = i\frac{\hbar}{T}\text{log}\left[\hat{U}(t_0 + T, t_0)\right],
$$

the logarithm of the one-period time evolution operator.

This library computes the eigenvalues $\varepsilon_n(\textbf{k})$ of $\hat{H}_F$ in the Brillouin zone (BZ) points $\textbf{k}$ for user defined driving fields and for the different approximations (methods) used to compute the interaction Hamiltonian $\hat{H}(t)$.

## Providing a crystal

The method described in Ref. [[1]](#ref1) supposes a crystal represented by its "exact tight-binding parameters". These are created by using [WannInt](https://github.com/irukoa/WannInt) crystals. Notice that this way to represent a physical system is quite flexible and can be used to represent real materials, models, and non-extended systems.

## Methods to compute $\hat{H}(t)$

In Ref. [[1]](#ref1) we describe 4 methods to compute $\hat{H}(t)$ from first principles,

1. Extended systems: length gauge, "no intraband" approximation. Non-extended systems: length gauge.

   $$
    \hat{H}(t) = \hat{H}_0 - q \sum_j E^j(t)\cdot \hat{X}^j.
   $$

2. Extended systems: velocity gauge, "no curvature" approximation. Non-extended systems: velocity gauge (exponential form).

   $$
    \hat{H}(t) = \text{exp}\left[\frac{iq}{\hbar}\sum_j A^j(t)\cdot \hat{X}^j\right]\hat{H}_0\;\text{exp}\left[\frac{-iq}{\hbar}\sum_j A^j(t)\cdot \hat{X}^j\right].
   $$

3. Extended systems: velocity gauge, 1st order approximation. Non-extended systems: ditto.

   $$
    \hat{H}(t) = \hat{H}_0 - \frac{q}{M}\sum_j A^j(t)\cdot \hat{p}^j.
   $$

4. Extended systems: velocity gauge, 2nd order approximation. Non-extended systems: ditto.

   $$
    \hat{H}(t) = \hat{H}_0 - \frac{q}{M}\sum_j A^j(t)\cdot \hat{p}^j + \frac{-q^2}{2\hbar^2}\sum_{jl} A^j(t)A^l(t)[\hat{\mathcal{D}}^j, [\hat{\mathcal{D}}^l, \hat{H}_0]].
   $$

Where $M$ is the particle mass, $\hat{\mathcal{D}}^j$ is the covariant derivative described in Ref. [[1]](#ref1), $A^j(t)$ is the vector potential,

$$
A^j(t) = -\int_{t'}^{t}E^j(\tau)d\tau,
$$

and $\hat{X}^j$ is defined in terms of the "matrix elements" of the Berry connection,

$$
X_{nm}^j(\textbf{k}) = (1 - \delta_{nm})\xi_{nm}^j(\textbf{k}).
$$

In non-extended systems, methods 1 and 2 are correct and 3 and 4 incorrect. In extended systems, methods 3 and 4 are increasingly better approximations to the real expression of $\hat{H}(t)$ and methods 1 and 2 are incorrect.

## Specifying $E^j(t)$

Once provided a crystal, all dynamics are solely determined by $E^j(t)$. In Ref. [[1]](#ref1) we consider the most general expression for a $T$-periodic $E^j(t)$ compatible with the Floquet conditions,

$$
E^j(t) = \sum_{n=1}^{N} \hat{E}_n^j\cos\left(n\omega t - \varphi_n^j\right)\;\textbf{e}_j,
$$

Here $N$ is the number of harmonics, $\hat{E}_n^j$ and $\varphi_n^j$ are the amplitude and phase corresponding to harmonic $n$ and component $j$, respectively and $\omega = 2\pi/T$.

# API

The module `Quasienergies` of Floquet defines the object
```fortran
type, extends(task_specifier) :: quasienergies_calc_tsk
```
which is an extension of the [SsTC](https://github.com/irukoa/SsTC_driver) sampling [task](https://github.com/irukoa/SsTC_driver?tab=readme-ov-file#typetask_specifier--tsk).

## `type(quasienergies_calc_tsk) :: tsk`

### Constructor

A Floquet quasienergies task is constructed by calling
```fortran
call tsk%build_floquet_task(crys, &
                            Nharm, &
                            axstart, axend, axsteps, &
                            pxstart, pxend, pxsteps, &
                            aystart, ayend, aysteps, &
                            pystart, pyend, pysteps, &
                            azstart, azend, azsteps, &
                            pzstart, pzend, pzsteps, &
                            omegastart, omegaend, omegasteps, &
                            t0start, t0end, t0steps[, &
                            Nt, htk_calc_method])
```
where

- `type(crystal), intent(in) :: crys` is an initialized [WannInt](https://github.com/irukoa/WannInt) crystal.
- `integer, intent(in) :: Nharm` is a positive integer specifying the number of harmonics of $E^j(t)$.
- `real(dp), intent(in) :: axstart(Nharm)` is an array of real numbers, each element `axstart(l)` containing the starting point of the amplitude $E^x_l$ corresponding to harmonic $l$ to consider in the calculation.
- `real(dp), intent(in) :: axend(Nharm)` is an array of real numbers, each element `axend(l)` containing the ending point of the amplitude $E^x_l$ corresponding to harmonic $l$ to consider in the calculation.
- `integer, intent(in) :: axsteps(Nharm)` is an array of positive integers, each element `axsteps(l)` containing the number of steps in the discretization of the amplitude $E^x_l$ corresponding to harmonic $l$ to consider in the calculation. The variable is discretized according to

   $$
    E^x_l(m) = E^x_l(1) + [E^x_l(M) - E^x_l(1)]\frac{m - 1}{M - 1},
   $$

  where $E^x_l(1)$ = `axstart(l)`, $E^x_l(M)$ = `axend(l)`, $M$ = `axsteps(l)`. If `axsteps(l)` = 1, then `axend(l)` = `axstart(l)`.

- `pxstart, pxend, pxsteps`: same data type and meaning as `axstart, axend, axsteps` except that they describe the phases $\varphi_l^x$.
- `aystart, ayend, aysteps`, `pystart, pyend, pysteps`, `azstart, azend, azsteps`, `pzstart, pzend, pzsteps`: as `axstart, axend, axsteps` and `pxstart, pxend, pxsteps`, respectively except that they describe the variables $E^y_l$, $\varphi_l^y$, $E^z_l$, $\varphi_l^z$, respectively.
- `real(dp), intent(in) :: omegastart` is a real number corresponding to the starting point of the frequency $\omega$ to consider in the calculation.
- `real(dp), intent(in) :: omegaend` is a real number corresponding to the ending point of the frequency $\omega$ to consider in the calculation.
- `integer, intent(in) :: omegasteps` is a positive integer, containing the number of steps in the discretization of the frequency $\omega$ to consider in the calculation. The variable is discretized as all other amplitudes and phases.
- `real(dp), intent(in) :: t0start` is a real number corresponding to the starting point of the initial time $t_0$ to consider in the calculation.
- `real(dp), intent(in) :: t0end` is a real number corresponding to the ending point of the initial time $t_0$ to consider in the calculation.
- `integer, intent(in) :: t0steps` is a positive integer, containing the number of steps in the discretization of the initial time $t_0$ to consider in the calculation. The variable is discretized as all other amplitudes and phases.
- `integer, intent(in), optional :: Nt` is a integer containing the number of points to discretize the time interval $[t_0, t_0 + T]$ in the calculation of the one period time evolution operator. Default is 513.
- `integer, intent(in), optional :: htk_calc_method` is an integer specifying the method to use in the calculation of $\hat{H}(t)$. The possibilities are
  - `-1`: length gauge, "no intraband" approximation.
  - `0`: velocity gauge, "no curvature" approximation.
  - `1`: velocity gauge, 1st order approximation.
  - `2`: velocity gauge, 2nd order approximation.

  Default is `1`.

### Integer and continuous indices

In the language of [SsTC](https://github.com/irukoa/SsTC_driver), a task has a number of integer and continuous indices. Functionally, the quasienergies $\varepsilon_n(\textbf{k})$ depend on a number of external parameters:

- Number of eigenvalues of the crystal $M$. These are passed as integer indices.
- Driving parameters $E^{\{x, y, z\}}_l$, $\varphi_l^{\{x, y, z\}}$ for each harmonic $l\in[1, N]$. Passed as continuous indices.
- Frequency $\omega$. Passed as continuous index.
- Starting time $t_0$. Passed as continuous index.

Each of the $6\times N + 2$ continuous indices can be discretized in a number of steps by providing a suitable `*steps(Nharm)` entry.

### Sampling

An initialized quasienergy calculation task can be sampled in a set of points in the BZ by using the standard SsTC [sampling](https://github.com/irukoa/SsTC_driver?tab=readme-ov-file#sampler) routine.

### Time steps handle

Is called as,
```fortran
Nt = tsk%ntsteps()
```
where `integer :: Nt` is the number of steps $N_t$ in the discretization of $[t_0, t_0 + T]$ in the calculation of the time evolution operator.

### Calculation method handle

Is called as,
```fortran
c_way = tsk%htk_calc_method()
```
where `integer :: c_way` is the method $[-1, 2]$ used to calculate $\hat{H}(t)$.

### Floquet initialization query

Is called as,
```fortran
initialized  = tsk%is_floquet_initialized()
```
where `logical :: initialized ` is `.true.` if the task has been initialized and `.false.` otherwise.

# Build

An automated build is available for [Fortran Package Manager](https://fpm.fortran-lang.org/) users. This is the recommended way to build and use Floquet in your projects. You can add Floquet to your project dependencies by including

```
[dependencies]
Floquet = { git="https://github.com/irukoa/Floquet.git" }
```
to the `fpm.toml` file.

# Usage

The library includes three examples. These programs replicate the data used to create the plots in Ref. [[1]](#ref1).

1. `Particle_in_a_Box`. This example simulates the particle in a box and calculates the quasienergy spectrum for variable driving field.
2. `BC2N_Driving_Scan`. This example simulates the material BC$_2$N and calculates the quasienergy spectrum for variable $E^x$ of a monochromatic, X-linear driving field with frequency $\hbar\omega = 0.5\;\text{eV}$ in the BZ point $\textbf{k} = (1/2, 1/4, 0)$.
3. `BC2N_Kpath`. This example simulates the material BC$_2$N and calculates the quasienergy spectrum for $E^x = 5.2\times10^8\;\text{V/m}$ of a monochromatic, X-linear driving field with frequency $\hbar\omega = 0.5\;\text{eV}$ in the BZ path $\Gamma$ - $X$ - $S$ - $Y$ - $\Gamma$.

These can be run by

```
fpm run --example Particle_in_a_Box
```
```
fpm run --example BC2N_Driving_Scan
```
```
fpm run --example BC2N_Kpath
```

respectively.

# TODO
- Revision.
- Banners.
- Upload to zenodo.
- CITATION/LICENSE.
- Reference main publication.
- Fix readme on github.

<a id="ref1"></a>
[1] A. R. Puente-Uriona *et al*. In preparation.