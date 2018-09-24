
<a id='Spatiotemporal-Timeseries-Prediction-1'></a>

# Spatiotemporal Timeseries Prediction


An application and extension of [local modeling](tsprediction/localmodels) to spatiotemporal timeseries.


!!! tip "Examples"
    Several example scripts can be found in `TimeseriesPrediction/examples`. These examples are run in the [examples](stexamples.md) page.



<a id='Spatio-Temporal-Embeddings-1'></a>

## Spatio-Temporal Embeddings


some info here.

<a id='TimeseriesPrediction.AbstractSpatialEmbedding' href='#TimeseriesPrediction.AbstractSpatialEmbedding'>#</a>
**`TimeseriesPrediction.AbstractSpatialEmbedding`** &mdash; *Type*.



```
AbstractSpatialEmbedding <: AbstractEmbedding
```

Super-type of spatiotemporal embedding methods. Valid subtypes:

  * `SpatioTemporalEmbedding`
  * `PCAEmbedding`

<a id='TimeseriesPrediction.SpatioTemporalEmbedding' href='#TimeseriesPrediction.SpatioTemporalEmbedding'>#</a>
**`TimeseriesPrediction.SpatioTemporalEmbedding`** &mdash; *Type*.



```
SpatioTemporalEmbedding{T,Φ,BC,X} → embedding
```

A spatio temporal delay coordinates structure to be used as a functor. Applies to data of `Φ` spatial dimensions and gives an embedding of dimensionality `X`.

```
embedding(rvec, s, t, α)
```

Operates inplace on `rvec` (of length `X`) and reconstructs vector from spatial timeseries `s` at timestep `t` and cartesian index `α`. Note that there are no bounds checks for `t`.

It is assumed that `s` is a `Vector{<:AbstractArray{T,Φ}}`.

**Constructors**

There are some convenience constructors that return intuitive embeddings here:

  * [`cubic_shell_embedding`](spatiotemporal.md#TimeseriesPrediction.cubic_shell_embedding)
  * [`light_cone_embedding`](spatiotemporal.md#TimeseriesPrediction.light_cone_embedding)

The "main" constructor is

```
SpatioTemporalEmbedding{X}(τ, β, bc, fsize)
```

which allows full control over the spatio-temporal embedding.

  * `Χ == length(τ) == length(β)` : dimensionality of resulting reconstructed space.
  * `τ::Vector{Int}` = Vector of temporal delays *for each entry* of the reconstructed space (sorted in ascending order).
  * `β::Vector{CartesianIndex{Φ}}` = vector of *relative* indices of spatial delays *for each entry* of the reconstructed space.
  * `bc::BC` : boundary condition.
  * `fsize::NTuple{Φ, Int}` : Size of each state in the timeseries.

<a id='TimeseriesPrediction.cubic_shell_embedding' href='#TimeseriesPrediction.cubic_shell_embedding'>#</a>
**`TimeseriesPrediction.cubic_shell_embedding`** &mdash; *Function*.



```
cubic_shell_embedding(s, D, τ, B, k, bc) → embedding
```

Create a [`SpatioTemporalEmbedding`](spatiotemporal.md#TimeseriesPrediction.SpatioTemporalEmbedding) instance that includes spatial neighbors in hypercubic *shells*. The embedding is to be used with data from `s`.

**Description**

Points are participating in the embedding by forming hypercubic shells around the current point. The total shells formed are `B`. The points on the shells have spatial distance `k ≥ 1` (distance in indices, like a cityblock metric). `k = 1` means that all points of the shell participate. The points of the hypercubic grid can be separated by `k ≥ 1` points apart (i.e. dropping `k-1` in-between points). In short, in each spatial dimension of the system the cartesian offset indices are `-B*k : k : k*B`.

`D` is the number of temporal steps in the past to be included in the embedding, where each step in the past has additional delay time `τ::Int`. `D=0` corresponds to using only the present. Notice that **all** embedded time frames have the same spatial structure, in contrast to [`light_cone_embedding`](spatiotemporal.md#TimeseriesPrediction.light_cone_embedding).

As an example, consider one of the `D` embedded frames (all are the same) of a system with 2 spatial dimensions (`□` = current point, (included *by definition* in the embedding), `n` = included points in the embedding coming from `n`-th shell, `.` = points not included in the embedding)

```
      B = 2,  k = 1        |        B = 1,  k = 2        |        B = 2,  k = 2
                           |                             |
.  .  .  .  .  .  .  .  .  |  .  .  .  .  .  .  .  .  .  |  2  .  2  .  2  .  2  .  2
.  .  .  .  .  .  .  .  .  |  .  .  .  .  .  .  .  .  .  |  .  .  .  .  .  .  .  .  .
.  .  2  2  2  2  2  .  .  |  .  .  1  .  1  .  1  .  .  |  2  .  1  .  1  .  1  .  2
.  .  2  1  1  1  2  .  .  |  .  .  .  .  .  .  .  .  .  |  .  .  .  .  .  .  .  .  .
.  .  2  1  □  1  2  .  .  |  .  .  1  .  □  .  1  .  .  |  2  .  1  .  □  .  1  .  2
.  .  2  1  1  1  2  .  .  |  .  .  .  .  .  .  .  .  .  |  .  .  .  .  .  .  .  .  .
.  .  2  2  2  2  2  .  .  |  .  .  1  .  1  .  1  .  .  |  2  .  1  .  1  .  1  .  2
.  .  .  .  .  .  .  .  .  |  .  .  .  .  .  .  .  .  .  |  .  .  .  .  .  .  .  .  .
.  .  .  .  .  .  .  .  .  |  .  .  .  .  .  .  .  .  .  |  2  .  2  .  2  .  2  .  2
```

<a id='TimeseriesPrediction.light_cone_embedding' href='#TimeseriesPrediction.light_cone_embedding'>#</a>
**`TimeseriesPrediction.light_cone_embedding`** &mdash; *Function*.



```
light_cone_embedding(s, D, τ, r₀, c, bc) → embedding
```

Create a [`SpatioTemporalEmbedding`](spatiotemporal.md#TimeseriesPrediction.SpatioTemporalEmbedding) instance that includes spatial and temporal neighbors of a point based on the notion of a *light cone*.

The embedding is to be used with data from `s`.

**Description**

Information does not travel instantly but with some finite speed `c ≥ 0.0`. This constructor creates a cone-like embedding including all points in space and time, whose value can influence a prediction based on the information speed `c`. `D` is the number of temporal steps in the past to be included in the embedding, where each step in the past has additional delay time `τ::Int`. `D=0` corresponds to using only the present. `r₀` is the initial radius at the top of the cone, i.e. the radius of influence at the present. `bc` is the boundary condition.

The radius of the light cone evolves as: `r = i*τ*c + r₀` for each step `i ∈ 0:D`.

As an example, in a one-dimensional system with `D = 1, τ = 2, r₀ = 1`, the embedding looks like (`□` = current point (included *by definition* in the embedding), `o` point to be predicted using [`temporalprediction`](spatiotemporal.md#TimeseriesPrediction.temporalprediction), `x` = points included in the embedding, `.` = points not included in the embedding)

```
time  | c = 1.0               | c = 2.0               | c = 0.0

n + 1 | ..........o.......... | ..........o.......... | ..........o..........
n     | .........x□x......... | .........x□x......... | .........x□x.........
n - 1 | ..................... | ..................... | .....................
n - τ | .......xxxxxxx....... | .....xxxxxxxxxx...... | .........xxx.........
```

Besides this example, in the official documentation we show a function `explain_light_cone` which produces a plot of the light cone for 2 spatial dimensions (great for understanding!).

<a id='TimeseriesPrediction.PCAEmbedding' href='#TimeseriesPrediction.PCAEmbedding'>#</a>
**`TimeseriesPrediction.PCAEmbedding`** &mdash; *Type*.



```
PCAEmbedding(s, em::SpatioTemporalEmbedding; kwargs...) → embedding
```

A spatio temporal delay coordinates structure with Principal Component Analysis as a means of dimension reduction, `embedding` can be used as a functor:

```julia
embedding(rvec, s, t, α)
```

which operates inplace on `rvec` and reconstructs vector from spatial time series `s` at timestep `t` and cartesian index `α`.

To instantiate this `embedding`, give the data to be reconstructed `s` as well as an instance of [`SpatioTemporalEmbedding`](spatiotemporal.md#TimeseriesPrediction.SpatioTemporalEmbedding) to `PCAEmbedding`.

**Keyword Arguments**

  * `pratio = 0.99` : Ratio of variances that needs to be preserved in low-dimensional PCA-reconstruction.
  * `maxoutdim = 25`: Upper limit for output dimension. May break `pratio` criterion.
  * `every_t = 1` : Speed up computation by only using every n-th point in time.
  * `every_α = 1` : Speed up computation further by only using every n-th point in space (linear indexing).

To set the output dimension to a certain value `X`, pass `pratio=1, maxoutdim=X`.


---


Boundary conditions

<a id='TimeseriesPrediction.ConstantBoundary' href='#TimeseriesPrediction.ConstantBoundary'>#</a>
**`TimeseriesPrediction.ConstantBoundary`** &mdash; *Type*.



```
ConstantBoundary(c) <: AbstractBoundaryCondition
```

Constant boundary condition type. Enforces constant boundary conditions when passed to [`SpatioTemporalEmbedding`](spatiotemporal.md#TimeseriesPrediction.SpatioTemporalEmbedding) by filling missing out-of-bounds values in the reconstruction with parameter `c`.

<a id='TimeseriesPrediction.PeriodicBoundary' href='#TimeseriesPrediction.PeriodicBoundary'>#</a>
**`TimeseriesPrediction.PeriodicBoundary`** &mdash; *Type*.



```
PeriodicBoundary <: AbstractBoundaryCondition
```

Periodic boundary condition struct. Enforces periodic boundary conditions when passed to [`SpatioTemporalEmbedding`](spatiotemporal.md#TimeseriesPrediction.SpatioTemporalEmbedding) in the reconstruction.


<a id='Prediction-functions-1'></a>

## Prediction functions

<a id='TimeseriesPrediction.temporalprediction' href='#TimeseriesPrediction.temporalprediction'>#</a>
**`TimeseriesPrediction.temporalprediction`** &mdash; *Function*.



```
temporalprediction(U, em::AbstractSpatialEmbedding, tsteps; kwargs...)
```

Perform a spatio-temporal time series prediction for `tsteps` iterations, using local weighted modeling [1] give a time series of the form `U::AbstractVector{<:AbstractArray{T, Φ}}`.

The returned data always contains the final state of `U` as starting point (total returned length is `tsteps+1`). The reconstruction process is defined by `em`. For available methods and interfaces see [`AbstractSpatialEmbedding`](spatiotemporal.md#TimeseriesPrediction.AbstractSpatialEmbedding).

**Keyword Arguments**

  * `ttype = KDTree` : Type/Constructor of tree structure. So far only tested with `KDTree`.
  * `method = AverageLocalModel(ω_safe)` : Subtype of [`AbstractLocalModel`](localmodels.md#TimeseriesPrediction.AbstractLocalModel).
  * `ntype = FixedMassNeighborhood(3)` : Subtype of [`AbstractNeighborhood`](../definition/dataset.md#DynamicalSystemsBase.AbstractNeighborhood).
  * `progress = true` : To print progress done.

**Description**

This method works similarly to [`localmodel_tsp`](localmodels.md#TimeseriesPrediction.localmodel_tsp), by expanding the concept of delay embedding to spatially extended systems. Instead of reconstructing complete states of the system, local states are used. See [`AbstractSpatialEmbedding`](spatiotemporal.md#TimeseriesPrediction.AbstractSpatialEmbedding) for details on the embedding. Predictions are then performed frame by frame and point py point. Once all values for a new frame are found, the frame is added to the end of the timeseries and used to generate new prediction queries for the next time step.

**Performance Notes**

Be careful when choosing embedding parameters as memory usage and computation time depend strongly on the resulting embedding dimension.

**References**

[1] : U. Parlitz & C. Merkwirth, [Phys. Rev. Lett. **84**, pp 1890 (2000)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.84.1890)

<a id='TimeseriesPrediction.crossprediction' href='#TimeseriesPrediction.crossprediction'>#</a>
**`TimeseriesPrediction.crossprediction`** &mdash; *Function*.



```
crossprediction(source_train, target_train, source_pred,
                em::AbstractSpatialEmbedding; kwargs...)
```

Perform a spatio-temporal timeseries cross-prediction for `target` from `source`, using local weighted modeling [1]. This can be used for example when there are coupled spatial fields and one is used to predict the other. It is assumed that `source_train`, `target_train`, `source_pred` are all of the same type, `AbstractVector{<:AbstractArray{T, Φ}}`.

The spatio temporal delay embedding process is defined by `em`. See [`AbstractSpatialEmbedding`](spatiotemporal.md#TimeseriesPrediction.AbstractSpatialEmbedding) for available methods and interfaces.

**Keyword Arguments**

  * `ttype = KDTree` : Type/Constructor of tree structure. So far only tested with `KDTree`.
  * `method = AverageLocalModel(ω_safe)` : Subtype of [`AbstractLocalModel`](localmodels.md#TimeseriesPrediction.AbstractLocalModel).
  * `ntype = FixedMassNeighborhood(3)` : Subtype of [`AbstractNeighborhood`](../definition/dataset.md#DynamicalSystemsBase.AbstractNeighborhood).
  * `progress = true` : To print progress done.

**Description**

The reconstructed state of `source_train[t][i,j,...]` is associated with the output value `target_train[t][i,j,...]`. This establishes a "connection" between `target` and `source`. Taking a reconstructed state of `source_pred` as query point, the function finds its neighbors in the reconstructed space of `source_train` using neighborhood `ntype`. Then, the neighbor *indices* are used to make a prediction for the corresponding value of the `target`, using the established "connection" between fields.

**Additional Interfaces**

To save computation time in the case of repeated predictions with the same training set and embedding parameters we provide an additional interface that allows the user to provide an existing reconstruction and tree structure.

```julia
R = reconstruct(train_in, em)
tree = ttype(R)
params = PredictionParams(em, method, ntype, ttype)
sol = crossprediction(params, train_out, pred_in, R, tree; progress=true)
```

where `params` is an internal container with all relevant parameters.

**Performance Notes**

Be careful when choosing embedding parameters as memory usage and computation time depend strongly on the resulting embedding dimension.

**References**

[1] : U. Parlitz & C. Merkwirth, [Phys. Rev. Lett. **84**, pp 1890 (2000)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.84.1890)

