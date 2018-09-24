
<a id='Lyapunov-Exponents-1'></a>

# Lyapunov Exponents


Lyapunov exponents measure rates of separation of nearby trajectories in the flow of a dynamical system. The [Wikipedia](https://en.wikipedia.org/wiki/Lyapunov_exponent) and the [Scholarpedia](http://www.scholarpedia.org/article/Lyapunov_exponent) entries have a lot of valuable information about the history and usage of these quantities.


This page treats systems where the equations of motion are known. If instead you have numerical data, see the [nonlinear timeseries analysis page](nlts).


<a id='Lyapunov-Spectrum-1'></a>

## Lyapunov Spectrum


The function `lyapunovs` calculates the entire spectrum of the Lyapunov exponents of a system:

<a id='ChaosTools.lyapunovs' href='#ChaosTools.lyapunovs'>#</a>
**`ChaosTools.lyapunovs`** &mdash; *Function*.



```
lyapunovs(ds::DynamicalSystem, N [, k::Int | Q0]; kwargs...) -> λs
```

Calculate the spectrum of Lyapunov exponents [1] of `ds` by applying a QR-decomposition on the parallelepiped matrix `N` times. Return the spectrum sorted from maximum to minimum.

The third argument `k` is optional, and dictates how many lyapunov exponents to calculate (defaults to `dimension(ds)`). Instead of passing an integer `k` you can pass a pre-initialized matrix `Q0` whose columns are initial deviation vectors (then `k = size(Q0)[2]`).

**Keyword Arguments**

  * `u0 = get_state(ds)` : State to start from.
  * `Ttr = 0` : Extra "transient" time to evolve the system before application of the algorithm. Should be `Int` for discrete systems. Both the system and the deviation vectors are evolved for this time.
  * `dt` : Time of individual evolutions between successive orthonormalization steps. Defaults to `1`. For continuous systems this is approximate.
  * `diffeq...` : Keyword arguments propagated into `init` of DifferentialEquations.jl. See [`trajectory`](../definition/evolve.md#DynamicalSystemsBase.trajectory) for examples. Only valid for continuous systems.

**Description**

The method we employ is "H2" of [2], originally stated in [3]. The deviation vectors defining a `D`-dimensional parallepiped in tangent space are evolved using the tangent dynamics of the system. A QR-decomposition at each step yields the local growth rate for each dimension of the parallepiped. The growth rates are then averaged over `N` successive steps, yielding the lyapunov exponent spectrum (at each step the parallepiped is re-normalized).

**Performance Notes**

This function uses a [`tangent_integrator`](../advanced.md#DynamicalSystemsBase.tangent_integrator). For loops over initial conditions and/or parameter values one should use the low level method that accepts an integrator, and `reinit!` it to new initial conditions. See the "advanced documentation" for info on the integrator object. The low level method is

```
lyapunovs(tinteg, N, dt::Real, Ttr::Real)
```

If you want to obtain the convergence timeseries of the Lyapunov spectrum, use the method

```
ChaosTools.lyapunovs_convergence(tinteg, N, dt, Ttr)
```

(not exported).

**References**

[1] : A. M. Lyapunov, *The General Problem of the Stability of Motion*, Taylor & Francis (1992)

[2] : K. Geist *et al.*, Progr. Theor. Phys. **83**, pp 875 (1990)

[3] : G. Benettin *et al.*, Meccanica **15**, pp 9-20 & 21-30 (1980)


---


As you can see, the documentation string is detailed and self-contained. For example, the lyapunov spectrum of the [folded towel map](http://www.scholarpedia.org/article/Hyperchaos) is calculated as:


```julia
using DynamicalSystems

ds = Systems.towel()
λλ = lyapunovs(ds, 10000)
```

```
3-element Array{Float64,1}:
  0.43246568251305534
  0.3721212450816193
 -3.296959221885987
```


Similarly, for a continuous system, e.g. the Lorenz system, you would do:


```julia
lor = Systems.lorenz(ρ = 32.0) #this is not the original parameter!
λλ = lyapunovs(lor, 10000, dt = 0.1)
```

```
3-element Array{Float64,1}:
   0.9820612905820211
  -0.00018297782315066568
 -14.64854497909885
```


`lyapunovs` is also very fast:


```julia
using BenchmarkTools
ds = Systems.towel()
@btime lyapunovs($ds, 2000);
```


```
  237.226 μs (45 allocations: 4.27 KiB)
```


<a id='Maximum-Lyapunov-Exponent-1'></a>

## Maximum Lyapunov Exponent


It is possible to get only the maximum Lyapunov exponent simply by giving `1` as the third argument of [`lyapunovs`](lyapunovs.md#ChaosTools.lyapunovs). However, there is a second algorithm that allows you to do the same thing, which is offered by the function `lyapunov`:

<a id='ChaosTools.lyapunov' href='#ChaosTools.lyapunov'>#</a>
**`ChaosTools.lyapunov`** &mdash; *Function*.



```
lyapunov(ds::DynamicalSystem, Τ; kwargs...) -> λ
```

Calculate the maximum Lyapunov exponent `λ` using a method due to Benettin [1], which simply evolves two neighboring trajectories (one called "given" and one called "test") while constantly rescaling the test one. `T`  denotes the total time of evolution (should be `Int` for discrete systems).

**Keyword Arguments**

  * `Ttr = 0` : Extra "transient" time to evolve the trajectories before starting to measure the expontent. Should be `Int` for discrete systems.
  * `d0 = 1e-9` : Initial & rescaling distance between the two neighboring trajectories.
  * `upper_threshold = 1e-6` : Upper distance threshold for rescaling.
  * `lower_threshold = 1e-12` : Lower distance threshold for rescaling (in order to  be able to detect negative exponents).
  * `dt = 1` : Time of evolution between each check of distance exceeding the thresholds. For continuous systems this is approximate.
  * `inittest = (u1, d0) -> u1 .+ d0/sqrt(D)` : A function that given `(u1, d0)` initializes the test state with distance `d0` from the given state `u1` (`D` is the dimension of the system). This function can be used when you want to avoid the test state appearing in a region of the phase-space where it would have e.g. different energy or escape to infinity.
  * `diffeq...` : Keyword arguments propagated into `init` of DifferentialEquations.jl. See [`trajectory`](../definition/evolve.md#DynamicalSystemsBase.trajectory) for examples. Only valid for continuous systems.

**Description**

Two neighboring trajectories with initial distance `d0` are evolved in time. At time $t_i$ their distance $d(t_i)$ either exceeds the `upper_threshold`, or is lower than `lower_threshold`, which initializes a rescaling of the test trajectory back to having distance `d0` from the given one, while the rescaling keeps the difference vector along the maximal expansion/contraction direction: $u_2 \to u_1+(u_2−u_1)/(d(t_i)/d_0)$.

The maximum Lyapunov exponent is the average of the time-local Lyapunov exponents

$$
\lambda = \frac{1}{t_{n}}\sum_{i=1}^{n}
\ln\left( a_i \right),\quad a_i = \frac{d(t_{i})}{d_0}.
$$

**Performance Notes**

This function uses a [`parallel_integrator`](../advanced.md#DynamicalSystemsBase.parallel_integrator). For loops over initial conditions and/or parameter values one should use the low level method that accepts an integrator, and `reinit!` it to new initial conditions. See the "advanced documentation" for info on the integrator object. The low level method is

```
lyapunov(pinteg, T, Ttr, dt, d0, ut, lt)
```

**References**

[1] : G. Benettin *et al.*, Phys. Rev. A **14**, pp 2338 (1976)


---


For example:


```julia
using DynamicalSystems
henon = Systems.henon()
λ = lyapunov(henon, 10000, d0 = 1e-7, upper_threshold = 1e-4, Ttr = 100)
```

```
0.42018736282059616
```


The same is done for continuous systems:


```julia
ross = Systems.lorenz(ρ = 32)
λ = lyapunov(ross, 10000.0, dt = 10.0, Ttr = 100.0)
```

```
0.9947201173035588
```

