
<a id='Numerical-Data-1'></a>

# Numerical Data


Numerical data in **DynamicalSystems.jl** is represented by a structure called `Dataset`

<a id='DynamicalSystemsBase.Dataset' href='#DynamicalSystemsBase.Dataset'>#</a>
**`DynamicalSystemsBase.Dataset`** &mdash; *Type*.



```
Dataset{D, T} <: AbstractDataset{D,T}
```

A dedicated interface for datasets. It contains *equally-sized datapoints* of length `D`, represented by `SVector{D, T}`.

When indexed with 1 index, a `dataset` is like a vector of datapoints.

When indexed with 2 indices it behaves like a matrix that has each of the columns be the timeseries of each of the dynamic variables.

**Description of indexing**

In the following let `i, j` be integers,  `typeof(data) <: AbstractDataset` and `v1, v2` be `<: AbstractVector{Int}` (`v1, v2` could also be ranges).

  * `data[i]` gives the `i`th datapoint (returns an `SVector`)
  * `data[v1]` will return a vector of datapoints
  * `data[v1, :]` using a `Colon` as a second index will return a `Dataset` of these points
  * `data[:, j]` gives the `j`th variable timeseries, as `Vector`
  * `data[v1, v2]` returns a `Dataset` with the appropriate entries (first indices being "time"/point index, while second being dynamic variables)
  * `data[i, j]` value of the `j`th variable, at the `i`th timepoint

Use `Matrix(dataset)` or `Dataset(matrix)` to convert. It is assumed that each *column* of the `matrix` is one dynamic variable. If you have various timeseries vectors `x, y, z, ...` pass them like `Dataset(x, y, z, ...)`. You can use `columns(dataset)` to obtain the reverse, i.e. all columns of the dataset in a tuple.


---


In essence a `Dataset` is simply a container for a `Vector` of `SVector`s. However, it is visually represented as a matrix, similarly to how numerical data would be printed on a spreadsheet (with time being the *column* direction). It also offers a lot more functionality than just pretty-printing. Besides the examples in the documentation string, you can also do:


```julia
using DynamicalSystems
hen = Systems.henon()
data = trajectory(hen, 10000) # this returns a dataset
for point in data
# do stuff with each datapoint
# (vector with as many elements as system dimension)
end
```


All functions from **DynamicalSystems.jl** that manipulate and use data are expecting an `AbstractDataset` subtype. This allows us to define efficient methods that coordinate well with other packages, like e.g. [`neighborhood`](dataset.md#DynamicalSystemsBase.neighborhood).


If given a matrix, we first convert to `Dataset`. This means that you should *first convert* your data to a `Dataset` if you want to call functions more than once, to avoid constantly converting.


<a id='Dataset-Functions-1'></a>

## Dataset Functions


Functions that operate on datasets.

<a id='DynamicalSystemsBase.minima' href='#DynamicalSystemsBase.minima'>#</a>
**`DynamicalSystemsBase.minima`** &mdash; *Function*.



```
minima(dataset)
```

Return an `SVector` that contains the minimum elements of each timeseries of the dataset.

<a id='DynamicalSystemsBase.maxima' href='#DynamicalSystemsBase.maxima'>#</a>
**`DynamicalSystemsBase.maxima`** &mdash; *Function*.



```
maxima(dataset)
```

Return an `SVector` that contains the maximum elements of each timeseries of the dataset.

<a id='DynamicalSystemsBase.minmaxima' href='#DynamicalSystemsBase.minmaxima'>#</a>
**`DynamicalSystemsBase.minmaxima`** &mdash; *Function*.



```
minmaxima(dataset)
```

Return `minima(dataset), maxima(dataset)` without doing the computation twice.

<a id='DynamicalSystemsBase.columns' href='#DynamicalSystemsBase.columns'>#</a>
**`DynamicalSystemsBase.columns`** &mdash; *Function*.



```
columns(dataset) -> x, y, z, ...
```

Return the individual columns of the dataset.


---


<a id='Dataset-I/O-1'></a>

## Dataset I/O


Input/output functionality for an `AbstractDataset` is already achieved using base Julia, specifically `writedlm` and `readdlm`.


The thing to note is that all data of an `AbstractDataset` is contained within its field `data`.


To write and read a dataset, simply do:


```julia
using DelimitedFiles

data = Dataset(rand(1000, 2))

# I will write and read using delimiter ','
writedlm("data.txt", data.data, ',')

# Don't forget to convert the matrix to a Dataset when reading
data = Dataset(readdlm("data.txt", ',', Float64))
```


---


<a id='Neighborhoods-in-a-Dataset-1'></a>

## Neighborhoods in a Dataset


Combining the excellent performance of [NearestNeighbors.jl](https://github.com/KristofferC/NearestNeighbors.jl) with the `AbstractDataset` allows us to define a function that calculates a "neighborhood" of a given point, i.e. finds other points near it. The different "types" of the neighborhoods are subtypes of `AbstractNeighborhood`.

<a id='DynamicalSystemsBase.neighborhood' href='#DynamicalSystemsBase.neighborhood'>#</a>
**`DynamicalSystemsBase.neighborhood`** &mdash; *Function*.



```
neighborhood(point, tree, ntype)
neighborhood(point, tree, ntype, n::Int, w::Int = 1)
```

Return a vector of indices which are the neighborhood of `point` in some `data`, where the `tree` was created using `tree = KDTree(data [, metric])`. The `ntype` is the type of neighborhood and can be any subtype of [`AbstractNeighborhood`](dataset.md#DynamicalSystemsBase.AbstractNeighborhood).

Use the second method when the `point` belongs in the data, i.e. `point = data[n]`. Then `w` stands for the Theiler window (positive integer). Only points that have index `abs(i - n) ≥ w` are returned as a neighborhood, to exclude close temporal neighbors. The default `w=1` is the case of excluding the `point` itself.

**References**

`neighborhood` simply interfaces the functions `knn` and `inrange` from [NearestNeighbors.jl](https://github.com/KristofferC/NearestNeighbors.jl) by using the argument `ntype`.

<a id='DynamicalSystemsBase.AbstractNeighborhood' href='#DynamicalSystemsBase.AbstractNeighborhood'>#</a>
**`DynamicalSystemsBase.AbstractNeighborhood`** &mdash; *Type*.



```
AbstractNeighborhood
```

Supertype of methods for deciding the neighborhood of points for a given point.

Concrete subtypes:

  * `FixedMassNeighborhood(K::Int)` : The neighborhood of a point consists of the `K` nearest neighbors of the point.
  * `FixedSizeNeighborhood(ε::Real)` : The neighborhood of a point consists of all neighbors that have distance < `ε` from the point.

See [`neighborhood`](dataset.md#DynamicalSystemsBase.neighborhood) for more.


---

