
<a id='Time-Evolution-of-Systems-1'></a>

# Time Evolution of Systems


!!! info "Trajectory and Timeseries"
    The word "timeseries" can be confusing, because it can mean a univariate (also called scalar or one-dimensional) timeseries or a multivariate (also called multi-dimensional) timeseries. To resolve this confusion, in **DynamicalSystems.jl** we have the following convention: **"timeseries"** always refers to a one-dimensional vector of numbers, which exists with respect to some other one-dimensional vector of numbers that corresponds to a time-vector. On the other hand, the word **"trajectory"** is used to refer to a *multi-dimensional* timeseries, which is of course simply a group/set of one-dimensional timeseries.

    Note that the data representation of a "trajectory" in Julia may vary: from a 2D Matrix to independent Vectors. In our package, a trajectory is always represented using a [`Dataset`](dataset.md#DynamicalSystemsBase.Dataset), which is a `Vector` of `SVector`s, and each `SVector` represents a data-point (the values of the variables at a given time-point).



**DynamicalSystems.jl** provides a convenient function for getting a trajectory of a system at equally spaced time points:

<a id='DynamicalSystemsBase.trajectory' href='#DynamicalSystemsBase.trajectory'>#</a>
**`DynamicalSystemsBase.trajectory`** &mdash; *Function*.



```
trajectory(ds::DynamicalSystem, T [, u]; kwargs...) -> dataset
```

Return a dataset that will contain the trajectory of the system, after evolving it for total time `T`, optionally starting from state `u`. See [`Dataset`](dataset.md#DynamicalSystemsBase.Dataset) for info on how to use this object.

A `W×D` dataset is returned, with `W = length(t0:dt:T)` with `t0:dt:T` representing the time vector (*not* returned) and `D` the system dimension. For discrete systems both `T` and `dt` must be integers.

**Keyword Arguments**

  * `dt` :  Time step of value output during the solving of the continuous system. For discrete systems it must be an integer. Defaults to `0.01` for continuous and `1` for discrete.
  * `Ttr` : Transient time to evolve the initial state before starting saving states.
  * `diffeq...` : Keyword arguments propagated into `init` of DifferentialEquations.jl. For example `abstol = 1e-9`.  Only valid for continuous systems. If you want to specify a solver, do so by using the name `alg`, e.g.: `alg = Tsit5(), maxiters = 1000`. This requires you to have been first `using OrdinaryDiffEq` to access the solvers. See `DynamicalSystemsBase.CDS_KWARGS` for default values. These keywords can also include `callback` for [event handling](http://docs.juliadiffeq.org/latest/features/callback_functions.html). Using a `SavingCallback` with `trajectory` will lead to unexpected behavior!


---


Notice that if you want to do repeated evolutions of different states of a continuous system, you should use the [`integrator`](../advanced.md#DynamicalSystemsBase.integrator) interface instead.


<a id='Solution-precision-for-continuous-systems-1'></a>

## Solution precision for continuous systems


A numerical solution of an ODE is not the "true" solution, uniquely defined by a (well-defined) ODE and an initial condition. Especially for chaotic systems, where deviations are amplified exponentially, one is left worried if the numerical solutions truly are part of the system and can truly give insight in understanding the system.


DifferentialEquations.jl offers a tool, called [Uncertainty Quantification](http://docs.juliadiffeq.org/latest/analysis/uncertainty_quantification.html), which allows users to asses up to what time-scales the numerical solution is close to the "true" solution. For example, using the default solving parameters of **DynamicalSystems.jl**, the Lorenz system is accurate up to time `t = 50.0`.


However, fortunately for us, there is not too much worry about the numerical solution diverging from the true solution. That is because of the [shadowing theorem](http://mathworld.wolfram.com/ShadowingTheorem.html) (or [shadowing lemma](http://www.scholarpedia.org/article/Shadowing_lemma_for_flows)):


!!! quote "Shadowing Theorem"
    Although a numerically computed chaotic trajectory diverges exponentially from the true trajectory with the same initial coordinates, there exists an errorless trajectory with a slightly different initial condition that stays near ("shadows") the numerically computed one.



This simply means that one can always numerically study chaos not only qualitatively but also quantitatively. For more information, see the book *Chaos in Dynamical Systems* by E. Ott, or the [scholarpedia](http://www.scholarpedia.org/article/Shadowing_lemma_for_flows) entry.

