[![Language](https://img.shields.io/badge/-Fortran-734f96?logo=fortran&logoColor=white)](https://github.com/topics/fortran)
ZENODO
[![Testing suite](https://github.com/irukoa/Floquet/actions/workflows/CI.yml/badge.svg)](https://github.com/irukoa/Floquet/actions/workflows/CI.yml)
# Floquet
This is a modern Fortran library, companion to Ref. [[1]](#ref1), which implements the computational scheme described in that reference to calculate the *effective Floquet Hamiltonian* of real materials and tight-binding models subjected to time-periodic driving fields.

<details>
  <summary>Details: working principle</summary>

## Working principle

Suppose that a crystal (real material or tight-binding model), is subjected to a time periodic driving field $E^j(t)$ with period $T$. Provided that the crystal is described by a Hamiltonian $\hat{H}_0$ and by a Berry connection $\hat{\xi}^j$, the interacting dynamics in the dipole approximation are described by the Hamiltonian

$$
\hat{H}(t) = \hat{H}_0 - q \sum_j E^j(t)\cdot \hat{r}^j,
$$

where $q$ is the particle charge and $\hat{r}^j$ is the position operator of the crystal.

The time evolution operator is defined as

$$
\hat{U}(t, t_0) = \text{exp}\left[\frac{-i}{\hbar}\int_{t_0}^{t_0 + t}\hat{H}(\tau)d\tau\right].
$$

This operator describes the propagation of quantum mechanical states in time. The relevant time-averaged physical properties contained in this operator can be analysed by considering the effective Floquet Hamiltonian,

$$
\hat{H}_F[t_0] = i\frac{\hbar}{T}\text{log}\left[\hat{U}(t_0 + T, t_0)\right],
$$

the logarithm of the one-period time evolution operator.

This library computes $\hat{H}_F$ in the Brillouin zone (BZ) points $\textbf{k}$ for user defined driving fields and for the different approximations (methods) used to compute the interaction Hamiltonian $\hat{H}(t)$.

In practice, the calculation of the one-period time evolution operator is implemented as

$$
\hat{U}(t, t_0) = \prod_{j=1}^{N_t}\text{exp}\left[-\frac{i}{\hbar} \delta t \hat{H}(t_j)  \right]+ \mathcal{O}(\delta t^2),
$$

where

$$
t_j = t_0 + T\frac{j-1}{N_t-1} = t_0 + \delta t (j-1).
$$

### Methods to compute $\hat{H}(t)$

In Ref. [[1]](#ref1) we describe 6 methods to compute $\hat{H}(t)$ from first principles,

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
\hat{H}(t) = \hat{H}_0 - \frac{q}{M}\sum_j A^j(t)\cdot \hat{p}^j + \frac{q^2}{2\hbar^2} \sum_j \sum_l A^j(t)A^l(t)[\hat{\mathcal{D}}^j, [\hat{\mathcal{D}}^l, \hat{H}_0]].
$$

5. Extended systems: velocity gauge, 3rd order approximation. Non-extended systems: ditto.

$$
\hat{H}(t) = \hat{H}_0 - \frac{q}{M}\sum_j A^j(t)\cdot \hat{p}^j + \frac{q^2}{2\hbar^2} \sum_j \sum_l A^j(t)A^l(t)[\hat{\mathcal{D}}^j, [\hat{\mathcal{D}}^l, \hat{H}_0]] +
$$

$$
\frac{-q^3}{6\hbar^3} \sum_p \sum_j \sum_l A^p(t)A^j(t)A^l(t)[\hat{\mathcal{D}}^p, [\hat{\mathcal{D}}^j, [\hat{\mathcal{D}}^l, \hat{H}_0]]]
$$

6. Extended systems: velocity gauge, 4th order approximation. Non-extended systems: ditto.

$$
\hat{H}(t) = \hat{H}_0 - \frac{q}{M}\sum_j A^j(t)\cdot \hat{p}^j + \frac{q^2}{2\hbar^2} \sum_j \sum_l A^j(t)A^l(t)[\hat{\mathcal{D}}^j, [\hat{\mathcal{D}}^l, \hat{H}_0]] +
$$

$$
\frac{-q^3}{6\hbar^3} \sum_p \sum_j \sum_l A^p(t)A^j(t)A^l(t)[\hat{\mathcal{D}}^p, [\hat{\mathcal{D}}^j, [\hat{\mathcal{D}}^l, \hat{H}_0]]] +
$$

$$
\frac{q^4}{24\hbar^4} \sum_o \sum_p \sum_j \sum_l A^o(t)A^p(t)A^j(t)A^l(t)[\hat{\mathcal{D}}^o,[\hat{\mathcal{D}}^p, [\hat{\mathcal{D}}^j, [\hat{\mathcal{D}}^l, \hat{H}_0]]]]
$$

Where $M$ is the particle mass, $\hat{p}^j$ is the momentum operator, $\hat{\mathcal{D}}^j$ is the covariant derivative described in Ref. [[1]](#ref1), $A^j(t)$ is the vector potential,

$$
A^j(t) = -\int_{t'}^{t}E^j(\tau)d\tau,
$$

and $\hat{X}^j$ is defined in terms of the "matrix elements" of the Berry connection,

$$
X_{nm}^j(\textbf{k}) = (1 - \delta_{nm})\xi_{nm}^j(\textbf{k}).
$$

In non-extended systems, methods 1 and 2 are correct and 3, 4, 5 and 6 are increasingly better approximations to the real expression of $\hat{H}(t)$ given by methods 1 and 2. In extended systems, methods 3, 4, 5 and 6 are increasingly better approximations to the real expression of $\hat{H}(t)$ and methods 1 and 2 are unfit for calculations.

</details>

## Specifying a crystal

The method described in Ref. [[1]](#ref1) supposes a crystal represented by its "exact tight-binding parameters", $\hat{H}_0$ and $\hat{\xi}^j$. In practice, crystals are specified employing [WannInt](https://github.com/irukoa/WannInt)'s `crystal` object. Notice that this way to represent a physical system can be used to represent real materials, models, and non-extended systems.

## Specifying $E^j(t)$

Once provided a crystal, all dynamics are solely determined by $E^j(t)$. In Ref. [[1]](#ref1) we consider the most general expression for a $T$-periodic $E^j(t)$ compatible with the Floquet conditions,

$$
E^j(t) = \sum_{n=1}^{N} \hat{E}_n^j\cos\left(n\omega t - \varphi_n^j\right)\textbf{e}_j,
$$

where $N$ is the number of harmonics, $\hat{E}_n^j$ and $\varphi_n^j$ are the amplitude and phase corresponding to harmonic $n$ and component $j$, respectively and $\omega = 2\pi/T$.

# API

The module `Quasienergies` of Floquet defines the object
```fortran
type, extends(task_specifier) :: quasienergies_calc_tsk
```
which is an extension of the [SsTC](https://github.com/irukoa/SsTC_driver) sampling [task](https://github.com/irukoa/SsTC_driver?tab=readme-ov-file#typetask_specifier--tsk). This object is used to calculate the quasienergy spectrum $\varepsilon_n(\textbf{k})$ in the BZ in eV [[1]](#ref1).

The module `Currents`[*](#disc) of Floquet defines the object
```fortran
type, extends(task_specifier) :: currents_calc_tsk
```
which is an extension of the [SsTC](https://github.com/irukoa/SsTC_driver) sampling [task](https://github.com/irukoa/SsTC_driver?tab=readme-ov-file#typetask_specifier--tsk). This object is used to calculate the integrand of the Fourier transform of the current denisty $j^l(\lambda)$ generated by the driving field $\textbf{E}(t)$ in Angstroms [[2]](#ref2).

<details>
  <summary> Quasienergies </summary>

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
- `real(dp), intent(in) :: axstart(Nharm)` is an array of real numbers, each element `axstart(l)` contains the starting point of the amplitude $E^x_l$ corresponding to harmonic $l$ to consider in the calculation.
- `real(dp), intent(in) :: axend(Nharm)` is an array of real numbers, each element `axend(l)` contains the ending point of the amplitude $E^x_l$ corresponding to harmonic $l$ to consider in the calculation.
- `integer, intent(in) :: axsteps(Nharm)` is an array of positive integers, each element `axsteps(l)` contains the number of steps in the discretization of the amplitude $E^x_l$ corresponding to harmonic $l$ to consider in the calculation. The variable is discretized according to

$$
E^x_l(m) = E^x_l(1) + [E^x_l(M) - E^x_l(1)]\frac{m - 1}{M - 1},
$$

where $E^x_l(1)$ = `axstart(l)`, $E^x_l(M)$ = `axend(l)`, $M$ = `axsteps(l)`. If `axsteps(l)` = 1, then `axend(l)` = `axstart(l)`.

- `pxstart, pxend, pxsteps`: same data type and meaning as `axstart, axend, axsteps` except that they describe the phases $\varphi_l^x$.
- `aystart, ayend, aysteps`, `pystart, pyend, pysteps`, `azstart, azend, azsteps`, `pzstart, pzend, pzsteps`: same meaning as `axstart, axend, axsteps` and `pxstart, pxend, pxsteps` pairwise, except that they describe the variables $E^y_l$, $\varphi_l^y$, $E^z_l$, $\varphi_l^z$, respectively.
- `real(dp), intent(in) :: omegastart` is a real number corresponding to the starting point of the frequency $\omega$ to consider in the calculation.
- `real(dp), intent(in) :: omegaend` is a real number corresponding to the ending point of the frequency $\omega$ to consider in the calculation.
- `integer, intent(in) :: omegasteps` is a positive integer, containing the number of steps in the discretization of the frequency $\omega$ to consider in the calculation. The variable is discretized as all other amplitudes and phases.
- `real(dp), intent(in) :: t0start` is a real number corresponding to the starting point of the initial time $t_0$ to consider in the calculation.
- `real(dp), intent(in) :: t0end` is a real number corresponding to the ending point of the initial time $t_0$ to consider in the calculation.
- `integer, intent(in) :: t0steps` is a positive integer, containing the number of steps in the discretization of the initial time $t_0$ to consider in the calculation. The variable is discretized as all other amplitudes and phases.
- `integer, intent(in), optional :: Nt` is a integer containing the number of points in the discretization of $[t_0, t_0 + T]$ for the calculation of the one period time evolution operator. Default is 513 steps.
- `integer, intent(in), optional :: htk_calc_method` is an integer specifying the method to use in the calculation of $\hat{H}(t)$. The possibilities are
  - `-1`: length gauge, "no intraband" approximation.
  - `0`: velocity gauge, "no curvature" approximation.
  - `1`: velocity gauge, 1st order approximation.
  - `2`: velocity gauge, 2nd order approximation.
  - `3`: velocity gauge, 3rd order approximation.
  - `4`: velocity gauge, 4th order approximation.

  Default is `1`.

### Integer and continuous indices

In the language of [SsTC](https://github.com/irukoa/SsTC_driver), a task has a number of integer and continuous indices. Functionally, the quasienergies $\varepsilon_n(\textbf{k})$ depend on a number of external parameters aside of the BZ vector $\textbf{k}$:

- Number of eigenvalues of the crystal $M$. This is passed as an integer index taking $M$ values.
- Driving parameters $E^{\{x, y, z\}}_l$ (3 variables), $\varphi_l^{\{x, y, z\}}$  (3 variables) for each harmonic $l\in[1, N]$. Passed as $6\times N$ continuous indices.
- Frequency $\omega$. Passed as a continuous index.
- Starting time $t_0$. Passed as a continuous index.

Each of the $6\times N + 2$ continuous indices can be discretized in a number of steps by providing a suitable `*steps(Nharm)` entry. The continuous index labeling is the following:

- Labels $[6(l-1)+1, 6(l-1)+6]$ correspond to the variables $E^{x}_l$, $\varphi_l^{x}$, ..., $E^{z}_l$, $\varphi_l^{z}$ of the $l$th harmonic.
- Labels $6\times N + 1$ and $6\times N + 2$ correspond to $t_0$ and $\omega$, respectively.

### Sampling

An initialized quasienergy calculation task can be sampled in a set of points in the BZ by using the standard SsTC [sampling](https://github.com/irukoa/SsTC_driver?tab=readme-ov-file#sampler) routine.

### Time steps handle

Is called as,
```fortran
Nt = tsk%ntsteps()
```
where `integer :: Nt` is the number of steps $N_t$ in the discretization of $[t_0, t_0 + T]$ in the calculation of the one period time evolution operator.

### Calculation method handle

Is called as,
```fortran
c_way = tsk%htk_calc_method()
```
where `integer :: c_way` is the method $[-1, 4]$ used to calculate $\hat{H}(t)$.

### Floquet initialization query

Is called as,
```fortran
initialized  = tsk%is_floquet_initialized()
```
where `logical :: initialized ` is `.true.` if the task has been initialized and `.false.` otherwise.

</details>

<details>
  <summary> Currents[*](#disc) </summary>

## `type(currents_calc_tsk) :: tsk`

This task is simmilar to the previously specified `type(quasienergies_calc_tsk)`, we only document the differences.

A Floquet currents task is constructed by calling
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
                            t0start, t0end, t0steps, &
                            lambdastart, lambdaend, lambdasteps[, &
                            delta_smr, &
                            Nt, Ns, htk_calc_method])
```
where

- `real(dp), intent(in) :: lambdastart` is a real number corresponding to the starting point of the frequency $\lambda$ of $j^l(\lambda)$ to consider in the calculation.
- `real(dp), intent(in) :: lambdaend` is a real number corresponding to the ending point of the frequency $\lambda$ of $j^l(\lambda)$ to consider in the calculation.
- `integer, intent(in) :: lambdasteps` is a positive integer, containing the number of steps in the discretization of the of the frequency $\lambda$ of $j^l(\lambda)$ to consider in the calculation.. The variable is discretized as all other amplitudes and phases.
- `real(dp), intent(in), optional :: delta_smr` is a positive real number $\sigma$ specifying the smearing of the resonant delta function, $\sigma \times \hbar \omega$, in eV. Default is $0.04 \times \hbar \omega$ in eV.
- `integer, intent(in), optional :: Ns` is a integer $N_s$ specifying the number of harmonics $s\in[-N_s, N_s]$ in the calculation of $\hat{Q}_s$ for the calculation of the Fourier components of the time-periodic operator $\hat{P}(t)$. Default is $N_s = 10$ harmonics.

### Integer and continuous indices

In the language of [SsTC](https://github.com/irukoa/SsTC_driver), a task has a number of integer and continuous indices. Functionally, the currents $j^l(\lambda)$ depend on a number of external parameters aside of the BZ vector $\textbf{k}$:

- Cartesinan component $l$. This is passed as an integer index taking 3 values.
- Driving parameters $E^{\{x, y, z\}}_l$ (3 variables), $\varphi_l^{\{x, y, z\}}$  (3 variables) for each harmonic $l\in[1, N]$. Passed as $6\times N$ continuous indices.
- Frequency $\omega$. Passed as a continuous index.
- Starting time $t_0$. Passed as a continuous index.
- Frequency $\lambda$. Passed as a continuous index.

Each of the $6\times N + 3$ continuous indices can be discretized in a number of steps by providing a suitable `*steps(Nharm)` entry. The continuous index labeling is the following:

- Labels $[6(l-1)+1, 6(l-1)+6]$ correspond to the variables $E^{x}_l$, $\varphi_l^{x}$, ..., $E^{z}_l$, $\varphi_l^{z}$ of the $l$th harmonic.
- Labels $6\times N + 1$, $6\times N + 2$, and $6\times N + 3$ correspond to $t_0$, $\omega$ and $\lambda$, respectively.

### Sampling

An initialized current calculation task can be sampled in a set of points in the BZ by using the standard SsTC [sampling](https://github.com/irukoa/SsTC_driver?tab=readme-ov-file#sampler) routine.

### Number of harmonics handle

Is called as,
```fortran
Ns = tsk%nsharms()
```
where `integer :: Ns` is the number of harmonics $N_s$ in the in the calculation of $\hat{Q}_s$.

### Dirac delta smearing handle

Is called as,
```fortran
smr = tsk%smr()
```
where `real(dp) :: smr` is the smearing $\sigma$ of the Dirac delta function, $\sigma \times \hbar \omega$, employed in the calculation in units of eV.

</details>

# Build

An automated build is available for [Fortran Package Manager](https://fpm.fortran-lang.org/) users. This is the recommended way to build and use Floquet in your projects. You can add Floquet to your project dependencies by including

```
[dependencies]
Floquet = { git="https://github.com/irukoa/Floquet.git" }
```
in the `fpm.toml` file.

The repository can be downloaded by running
```
git clone https://github.com/irukoa/Floquet.git
```
in a directory of the users' choice.

# Usage

The repository includes four examples. These programs replicate the data used to create some of the plots in Ref. [[1]](#ref1).

1. `Particle_in_a_Box`. This example simulates the particle in a box and calculates the quasienergy spectrum for variable monochromatic driving field.
2. `BC2N_Driving_Scan`. This example simulates the material BC$_2$N and calculates the quasienergy spectrum for variable $E^x$ of a monochromatic, X-linear driving field with frequency $\hbar\omega = 0.5$ $\text{eV}$ in the BZ point $\textbf{k} = (1/2, 1/4, 0)$.
3. `BC2N_Kpath`. This example simulates the material BC$_2$N and calculates the quasienergy spectrum for $E^x = 5.2\times10^8$ $\text{V/m}$ of a monochromatic, X-linear driving field with frequency $\hbar\omega = 0.5$ $\text{eV}$ in the BZ path $\Gamma$ - $\text{X}$ - $\text{S}$ - $\text{Y}$ - $\Gamma$.
4. `BC2N_Current`[*](#disc). This example simulates current generated in the material BC$_2$N by a monochromatic, X-linear driving field with frequency $\hbar\omega = 2.5$ $\text{eV}$ and amplitude $E^x = 5.2\times10^8$ $\text{V/m}$.

These can be run using

```
fpm run --example Particle_in_a_Box
```
```
fpm run --example BC2N_Driving_Scan
```
```
fpm run --example BC2N_Kpath
```
```
fpm run --example BC2N_Current
```

respectively in the repository directory.

# TODO
- Upload to zenodo.
- Update CITATION
- When done, tag as v1.0.0 in fpm.toml, main.F90 and CITATION.
- Release.

<a id="ref1"></a>
[1] Á. R. Puente-Uriona, M. Modugno, I. Souza, y J. Ibañez-Azpiroz, «Computing Floquet quasienergies in finite and extended systems: Role of electromagnetic and quantum-geometric gauges», [Phys. Rev. B, vol. 110, n.º 2, p. 125203](https://doi.org/10.1103/PhysRevB.110.125203), sep. 2024, DOI: 10.1103/PhysRevB.110.125203.

<a id="ref2"></a>
[2] A. R. Puente-Uriona *et al*. In preparation.

<a id="disc"></a>
[*](#disc) The `Currents` module, and examples employing the source code therein are considered a work in progress related to the completion of Ref. [[2]](#ref2), and thus not ready for production calculations. Use it at your own discretion in compliance with the license agreement.
