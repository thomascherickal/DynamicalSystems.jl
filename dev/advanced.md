
<a id='Advanced-documentation-1'></a>

# Advanced documentation


This section overviews the various integrators available from `DynamicalSystemsBase`, as well as gives some insight into the internals, so that other developers that want to use this library can build upon it.


<a id='Integrators-1'></a>

## Integrators

<a id='DynamicalSystemsBase.integrator' href='#DynamicalSystemsBase.integrator'>#</a>
**`DynamicalSystemsBase.integrator`** &mdash; *Function*.



```
integrator(ds::DynamicalSystem [, u0]; diffeq...) -> integ
```

Return an integrator object that can be used to evolve a system interactively using `step!(integ [, Δt])`. Optionally specify an initial state `u0`.

The state of this integrator is a vector.

  * `diffeq...` are keyword arguments propagated into `init` of DifferentialEquations.jl. See [`trajectory`](definition/evolve.md#DynamicalSystemsBase.trajectory) for examples. Only valid for continuous systems.

<a id='DynamicalSystemsBase.parallel_integrator' href='#DynamicalSystemsBase.parallel_integrator'>#</a>
**`DynamicalSystemsBase.parallel_integrator`** &mdash; *Function*.



```
parallel_integrator(ds::DynamicalSystem, states; kwargs...)
```

Return an integrator object that can be used to evolve many `states` of a system in parallel at the *exact same times*, using `step!(integ [, Δt])`.

`states` are expected as vectors of vectors.

**Keyword Arguments**

  * `diffeq...` : Keyword arguments propagated into `init` of DifferentialEquations.jl. See [`trajectory`](definition/evolve.md#DynamicalSystemsBase.trajectory) for examples. Only valid for continuous systems. These keywords can also include `callback` for [event handling](http://docs.juliadiffeq.org/latest/features/callback_functions.html).

It is *heavily* advised to use the functions [`get_state`](advanced.md#DynamicalSystemsBase.get_state) and [`set_state!`](advanced.md#DynamicalSystemsBase.set_state!) to manipulate the integrator. Provide `i` as a second argument to change the `i`-th state.

<a id='DynamicalSystemsBase.tangent_integrator' href='#DynamicalSystemsBase.tangent_integrator'>#</a>
**`DynamicalSystemsBase.tangent_integrator`** &mdash; *Function*.



```
tangent_integrator(ds::DynamicalSystem, Q0 | k::Int; kwargs...)
```

Return an integrator object that evolves in parallel both the system as well as deviation vectors living on the tangent space.

`Q0` is a *matrix* whose columns are initial values for deviation vectors. If instead of a matrix `Q0` an integer `k` is given, then `k` random orthonormal vectors are choosen as initial conditions.

**Keyword Arguments**

  * `u0` : Optional different initial state.
  * `diffeq...` : Keyword arguments propagated into `init` of DifferentialEquations.jl. See [`trajectory`](definition/evolve.md#DynamicalSystemsBase.trajectory) for examples. Only valid for continuous systems. These keywords can also include `callback` for [event handling](http://docs.juliadiffeq.org/latest/features/callback_functions.html).

It is *heavily* advised to use the functions [`get_state`](advanced.md#DynamicalSystemsBase.get_state), [`get_deviations`](advanced.md#DynamicalSystemsBase.get_deviations), [`set_state!`](advanced.md#DynamicalSystemsBase.set_state!), [`set_deviations!`](advanced.md#DynamicalSystemsBase.set_deviations!) to manipulate the integrator.

**Description**

If $J$ is the jacobian of the system then the *tangent dynamics* are the equations that evolve in parallel the system as well as a deviation vector (or matrix) $w$:

$$
\begin{aligned}
\dot{u} &= f(u, p, t) \\
\dot{w} &= J(u, p, t) \times w
\end{aligned}
$$

with $f$ being the equations of motion and $u$ the system state. Similar equations hold for the discrete case.


---


Notice that the state type `integrator.u` of each integrator is quite different and *does change* between the possible versions of a [`DynamicalSystem`](definition/general.md#DynamicalSystemsBase.DynamicalSystem)!


<a id='Integrator-state-functions-1'></a>

## Integrator state functions


There are four functions associated with the integrators that we export:

<a id='DynamicalSystemsBase.get_state' href='#DynamicalSystemsBase.get_state'>#</a>
**`DynamicalSystemsBase.get_state`** &mdash; *Function*.



```
get_state(ds::DynamicalSystem)
```

Return the state of `ds`.

```
get_state(integ [, i::Int = 1])
```

Return the state of the integrator, in the sense of the state of the dynamical system.

If the integrator is a [`parallel_integrator`](advanced.md#DynamicalSystemsBase.parallel_integrator), passing `i` will return the `i`-th state. The function also correctly returns the true state of the system for tangent integrators.

<a id='DynamicalSystemsBase.set_state!' href='#DynamicalSystemsBase.set_state!'>#</a>
**`DynamicalSystemsBase.set_state!`** &mdash; *Function*.



```
set_state!(integ, u [, i::Int = 1])
```

Set the state of the integrator to `u`, in the sense of the state of the dynamical system. Works for any integrator (normal, tangent, parallel).

For parallel integrator, you can choose which state to set (using `i`).

Automatically does `u_modified!(integ, true)`.

<a id='DynamicalSystemsBase.get_deviations' href='#DynamicalSystemsBase.get_deviations'>#</a>
**`DynamicalSystemsBase.get_deviations`** &mdash; *Function*.



```
get_deviations(tang_integ)
```

Return the deviation vectors of the [`tangent_integrator`](advanced.md#DynamicalSystemsBase.tangent_integrator) in a form of a matrix with columns the vectors.

<a id='DynamicalSystemsBase.set_deviations!' href='#DynamicalSystemsBase.set_deviations!'>#</a>
**`DynamicalSystemsBase.set_deviations!`** &mdash; *Function*.



```
set_deviations!(tang_integ, Q)
```

Set the deviation vectors of the [`tangent_integrator`](advanced.md#DynamicalSystemsBase.tangent_integrator) to `Q`, which must be a matrix with each column being a deviation vector.

Automatically does `u_modified!(tang_integ, true)`.


!!! note
    These functions work with *any* possible integrator and it is best to use the to change states robustly!



<a id='Re-initializing-an-integrator-1'></a>

## Re-initializing an integrator


It is more efficient to re-initialize an integrator using `reinit!` than to create a new one. This can be very helpful when looping over initial conditions and/or parameter values.


All high-level functions from `ChaosTools` have a set-up part that creates an integrator, and a low-level part that does the computation. The low level part is your friend! Use it! See the [Using `gali`](chaos/chaos_detection/#using-gali) page for an example as well as the section below.


The `reinit!` call signature is the same for continuous and discrete systems. In the following, `state` is supposed to be a `D` dimensional vector (state of the dynamical system).


1. `reinit!(integ, state)` : to be used with standard [`integrator`](advanced.md#DynamicalSystemsBase.integrator).
2. `reinit!(integ, Vector_of_states)` : to be used with the [`parallel_integrator`](advanced.md#DynamicalSystemsBase.parallel_integrator).
3. `reinit!(integ, state, Q0::AbstractMatrix)` : to be used with the [`tangent_integrator`](advanced.md#DynamicalSystemsBase.tangent_integrator). This three argument version of `reinit!` is exported from `DynamicalSystemsBase`.


<a id='Re-init-of-continuous-tangent-integrator-1'></a>

### Re-init of continuous tangent integrator


Here we compute the [`lyapunovs`](chaos/lyapunovs.md#ChaosTools.lyapunovs) for many different initial conditions.


```julia
ds = Systems.lorenz()
tinteg = tangent_integrator(ds, 2)
ics = [rand(3) for i in 1:100]
for ic in ics
  reinit!(tinteg, ic, orthonormal(3, 2))
  λ = lyapunovs(tinteg, 1000, 0.1, 10.0)
  # reminder: lyapunovs(tinteg, N, dt::Real, Ttr::Real = 0.0)
end
```


<a id='Re-init-of-discrete-parallel-integrator-1'></a>

### Re-init of discrete parallel integrator


Here we compute the [`lyapunov`](chaos/lyapunovs.md#ChaosTools.lyapunov) for many different parameters.


```julia
ds = Systems.henon()
u0 = rand(SVector{2})
ps = 1.2:0.01:1.4
pinteg = parallel_integrator(ds, [u0, u0 + 1e-9rand(SVector{2})])
for p in ps
  set_parameter!(ds, 1, p)
  reinit!(pinteg, [u0, u0 + 1e-9rand(SVector{2})])
  λ = lyapunov(pinteg, 1000, 10, 1, 1e-9, 1e-6, 1e-12)
  # reminder: lyapunov(pinteg, T, Ttr, dt, d0, ut, lt)
end
```


<a id='Using-callbacks-with-integrators-1'></a>

## Using callbacks with integrators


For the case of continuous systems you can add callbacks from the event handling of **DifferentialEquations.jl**. This is done simply as a keyword argument to the initializers.


In this example we use a simple `SavingCallback` to save the distance between the two states of a [`parallel_integrator`](advanced.md#DynamicalSystemsBase.parallel_integrator).


```julia
using DynamicalSystems, DiffEqCallbacks
using LinearAlgebra: norm

kwargs = (abstol=1e-14, reltol=1e-14, maxiters=1e9)
ds = Systems.lorenz()
d0 = 1e-9
T = 100.0

save_func(u, t, integrator) = norm(u[1] - u[2])
saved_values = SavedValues(eltype(ds.t0), eltype(get_state(ds)))
cb = SavingCallback(save_func, saved_values)

u0 = get_state(ds)
pinteg = parallel_integrator(ds, [u0, u0 + rand(SVector{3})*d0*√3];
kwargs..., callback = cb)
step!(pinteg, T)
t = saved_values.t
n = saved_values.saveval
```

```
14119-element Array{Float64,1}:
  1.4981336453946657e-9
  1.5883069805406025e-9
  1.6576550536516179e-9
  1.7419750383332456e-9
  1.8307146213692537e-9
  1.9261408187977782e-9
  2.0265694731200064e-9
  2.1311256367541026e-9
  2.2397692692857393e-9
  2.3524406880868853e-9
  ⋮
 24.658822532770333
 24.629190494297905
 24.639717957836567
 24.694997432972695
 24.79767232922235
 24.952329454386128
 25.15835335949889
 25.42038458643568
 25.735360165327695
```


As expected you can see that the recorded distance between two states is increasing.


<a id='DynamicalSystem-implementation-1'></a>

## `DynamicalSystem` implementation


```julia
abstract type DynamicalSystem{
        IIP,     # is in place , for dispatch purposes and clarity
        S,       # state type
        D,       # dimension
        F,       # equations of motion
        P,       # parameters
        JAC,     # jacobian
        JM,      # jacobian matrix
        IAD}     # is auto-differentiated
    # one-liner: {IIP, S, D, F, P, JAC, JM, IAD}
    # Subtypes of DynamicalSystem have fields:
    # 1. f
    # 2. u0
    # 3. p
    # 4. t0
    # 5. jacobian (function)
    # 6. J (matrix)
end
```


The `DynamicalSystem` stores only the absolutely necessary information. Every other functionality of **DynamicalSystems.jl** initializes an integrator.


The final type-parameter `IAD` is useful when creating the `tangent_integrator`, so that the vector field is not computed twice!

