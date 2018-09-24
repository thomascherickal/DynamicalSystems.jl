
<a id='Entropies-and-Dimensions-1'></a>

# Entropies and Dimensions


<a id='Entropies-1'></a>

## Entropies


In the study of dynamical systems there are many quantities that identify as "entropy". Notice that these quantities are not the more commonly known [thermodynamic ones](https://en.wikipedia.org/wiki/Entropy), used in Statistical Physics. Rather, they are more like the to the entropies of [information theory](https://en.wikipedia.org/wiki/Entropy_(information_theory)), which represents information contained within a dataset, or information about the dimensional scaling of a dataset.


<a id='Generalized-Entropy-1'></a>

### Generalized Entropy

<a id='ChaosTools.genentropy' href='#ChaosTools.genentropy'>#</a>
**`ChaosTools.genentropy`** &mdash; *Function*.



```
genentropy(α, ε, dataset::AbstractDataset; base = e)
```

Compute the `α` order generalized (Rényi) entropy [1] of a dataset, by first partitioning it into boxes of length `ε` using [`non0hist`](entropies.md#ChaosTools.non0hist).

```julia
genentropy(α, p::AbstractArray; base = e)
```

Compute the entropy of an array `p` directly, assuming that `p` is sum-normalized.

Optionally use `base` for the logarithms.

**Description**

Let $p$ be an array of probabilities (summing to 1). Then the Rényi entropy is

$$
R_\alpha(p) = \frac{1}{1-\alpha} \log \left(\sum_i p[i]^\alpha\right)
$$

and generalizes other known entropies, like e.g. the information entropy ($\alpha = 1$, see [2]), the maximum entropy ($\alpha=0$, also known as Hartley entropy), or the correlation entropy ($\alpha = 2$, also known as collision entropy).

**References**

[1] : A. Rényi, *Proceedings of the fourth Berkeley Symposium on Mathematics, Statistics and Probability*, pp 547 (1960)

[2] : C. E. Shannon, Bell Systems Technical Journal **27**, pp 379 (1948)


---


Basically, given a [dataset](system_definition/#numerical-data) you can partition it into boxes to calculate an entropy.


!!! tip "Worried about memory overflow? Don't be!"
    Partitioning the dataset (i.e. doing a *histogram*) is in general a costly operation that depends exponentially on the number of dimensions of the data and algebraically to the box size `ε`. However, in this specific case the partition process has some special aspects that can be taken advantage of, reducing tremendously the memory allocation and spent time!

    In fact, there is an upper bound to the memory allocated by `non0hist`: A constant multiplied by the length of the array, `N = length(p)`. No matter how small `ε` or how many dimensions the data has, the method can at most assign `N` dictionary entries.



The function used internally by `genentropy` is `non0hist`:

<a id='ChaosTools.non0hist' href='#ChaosTools.non0hist'>#</a>
**`ChaosTools.non0hist`** &mdash; *Function*.



```julia
non0hist(ε, dataset::AbstractDataset)
```

Partition a dataset into tabulated intervals (boxes) of size `ε` and return the sum-normalized histogram in an unordered 1D form, discarding all zero elements.

**Performances Notes**

This method is effecient in both memory and speed, because it uses a dictionary to collect the information of bins with elements, while it completely disregards empty bins. This allows computation of entropies of high-dimensional datasets and with small box sizes `ε` without memory overflow.

Use e.g. `fit(Histogram, ...)` from [`StatsBase`](http://juliastats.github.io/StatsBase.jl/stable/) if you wish to keep information about the edges of the binning as well as the zero elements.


---


For example, the Shannon entropy of a coin-flip process should be one bit, [by definition](https://en.wikipedia.org/wiki/Shannon_(unit)). Let's see...


```julia
using DynamicalSystems
y = Float64.(rand(Bool, 1000000)) # just some coin tosses
sh = genentropy(1, 0.1, y)  # this is the shannon entropy
isapprox(sh, log(2),  rtol = 1e-4)
```

```
true
```


Because all entropies are by default calculated on base-$e$, the unit of measurement is "nat" and one bit is $\log(2)\times$nat.


<a id='Permutation-Entropy-1'></a>

### Permutation Entropy


The permutation entropy is introduced by C. Bandt and B. Pompe as a "A Natural Complexity Measure for Timeseries", which directly applies to arbitrary real-world data and is particularly useful in the presence of dynamical or observational noise.

<a id='ChaosTools.permentropy' href='#ChaosTools.permentropy'>#</a>
**`ChaosTools.permentropy`** &mdash; *Function*.



```
permentropy(x::AbstractVector, order [, interval=1]; base = e)
```

Compute the permutation entropy [1] of given `order` from the `x` timeseries.

Optionally, `interval` can be specified to use `x[t0:interval:t1]` when calculating permutation of the sliding windows between `t0` and `t1 = t0 + interval * (order - 1)`.

Optionally use `base` for the logarithms.

**References**

[1] : C. Bandt, & B. Pompe, [Phys. Rev. Lett. **88** (17), pp 174102 (2002)](http://doi.org/10.1103/PhysRevLett.88.174102)


For example, we will compute and compare the [`lyapunov`](lyapunovs.md#ChaosTools.lyapunov) exponent of the logistic map with the order-6 permutation entropy, like in the original paper.


```julia
ds = Systems.logistic()
rs = 3.5:0.001:4
ls = Float64[]; hs = Float64[]
for r in rs
    ds.p[1] = r
    push!(ls, lyapunov(ds, 100000))
    # For 1D systems `trajectory` returns a vector
    push!(hs, permentropy(trajectory(ds, 10000), 6))
end

f = figure(figsize = (10,6))
a1 = subplot(211)
plot(rs, ls); ylim(-2, log(2)); ylabel("\$\\lambda\$")
a1[:axes][:get_xaxis]()[:set_ticklabels]([])
xlim(rs[1], rs[end]);

a2 = subplot(212)
plot(rs, hs; color = "C1"); ylabel("\$h_6\$")
xlim(rs[1], rs[end]); xlabel("\$r\$")
tight_layout()
```


![Permutation Entropy](https://i.imgur.com/tsqSA7a.png)


!!! info "Permutation Entropy performance"
    Even though the current implementation is fine and runs reasonably fast for moderate orders, it can get slow for high orders. Issue [ChaosTools.jl#22](https://github.com/JuliaDynamics/ChaosTools.jl/issues/22) keeps track of this, and contains information on how to improve performance.



<a id='Attractor-Dimension-Estimation-1'></a>

## Attractor Dimension Estimation


There are numerous methods that one can use to calculate a so-called "dimension" of a dataset, like for example the [Fractal dimension](https://en.wikipedia.org/wiki/Fractal_dimension). This real number can offer a lot of information about the object that the dataset represents.


<a id='Generalized-Dimensions-1'></a>

### Generalized Dimensions


Based on the definition of the [generalized entropy](#ChaosTools.genentropy), one can calculate an appropriate dimension, called *generalized dimension*:

<a id='ChaosTools.generalized_dim' href='#ChaosTools.generalized_dim'>#</a>
**`ChaosTools.generalized_dim`** &mdash; *Function*.



```
generalized_dim(α, dataset [, sizes]) -> D_α
```

Return the `α` order generalized dimension of the `dataset`, by calculating the [`genentropy`](entropies.md#ChaosTools.genentropy) for each `ε ∈ sizes`.

**Description**

The returned dimension is approximated by the (inverse) power law exponent of the scaling of the [`genentropy`](entropies.md#ChaosTools.genentropy) versus the box size `ε`, where `ε ∈ sizes`.

Calling this function performs a lot of automated steps:

1. A vector of box sizes is decided by calling `sizes = estimate_boxsizes(dataset)`, if `sizes` is not given.
2. For each element of `sizes` the appropriate entropy is calculated, through `d = genentropy.(α, sizes, dataset)`. Let `x = -log.(sizes)`.
3. The curve `d(x)` is decomposed into linear regions, using [`linear_regions`](entropies.md#ChaosTools.linear_regions)`(x, d)`.
4. The biggest linear region is chosen, and a fit for the slope of that region is performed using the function [`linear_region`](entropies.md#ChaosTools.linear_region). This slope is the return value of `generalized_dim`.

By doing these steps one by one yourself, you can adjust the keyword arguments given to each of these function calls, refining the accuracy of the result.

The following aliases are provided:

  * α = 0 : `boxcounting_dim`, `capacity_dim`
  * α = 1 : `information_dim`


---


!!! danger "Be wary when using `generalized_dim`"
    As stated clearly by the documentation string, calling `generalized_dim` performs a lot of automated steps by calling other functions (see below) with default arguments. It is actually more like a convenient bundle than an actual function and therefore you should be careful when considering the validity of the returned number.


<a id='ChaosTools.estimate_boxsizes' href='#ChaosTools.estimate_boxsizes'>#</a>
**`ChaosTools.estimate_boxsizes`** &mdash; *Function*.



```
estimate_boxsizes(dataset::AbstractDataset; k::Int = 12, z = -1, w = 1)
```

Return `k` exponentially spaced values from 10^`lower + w` to 10^`upper + z`.

`lower` is the magnitude of the minimum pair-wise distance between datapoints while `upper` is the magnitude of the maximum difference between greatest and smallest number among each timeseries.

"Magnitude" here stands for order of magnitude, i.e. `round(log10(x))`.

<a id='ChaosTools.linear_regions' href='#ChaosTools.linear_regions'>#</a>
**`ChaosTools.linear_regions`** &mdash; *Function*.



```
linear_regions(x, y; dxi::Int = 1, tol = 0.2) -> (lrs, tangents)
```

Identify regions where the curve `y(x)` is linear, by scanning the `x`-axis every `dxi` indices (e.g. at `x[1] to x[5], x[5] to x[10], x[10] to x[15]` and so on if `dxi=5`).

If the slope (calculated using `LsqFit`) of a region of width `dxi` is approximatelly equal to that of the previous region, within tolerance `tol`, then these two regions belong to the same linear region.

Return the indices of `x` that correspond to linear regions, `lrs`, and the approximated `tangents` at each region. `lrs` is a vector of `Int`.

A function `plot_linear_regions` visualizes the result of using this `linear_regions` (requires `PyPlot`).

<a id='ChaosTools.linear_region' href='#ChaosTools.linear_region'>#</a>
**`ChaosTools.linear_region`** &mdash; *Function*.



```
linear_region(x, y; dxi::Int = 1, tol = 0.2) -> ([ind1, ind2], slope)
```

Call [`linear_regions`](entropies.md#ChaosTools.linear_regions), identify the largest linear region and approximate the slope of the entire region using `linreg`. Return the indices where the region starts and stops (`x[ind1:ind2]`) as well as the approximated slope.


---


<a id='Example-1'></a>

#### Example


For example, the dimension of the strange attractor of the [Hénon map](system_definition/#DynamicalSystems.Systems.henon) is:


```julia
using DynamicalSystems
hen = Systems.henon()
ts = trajectory(hen, 1000000)
D_hen = information_dim(ts)
```


```
1.2279316105815665
```


As a side note, be sure that you have enough data points, otherwise the values you will get will never be correct, as is demonstrated by J.-P. Eckmann and D. Ruelle (see Physica D **56**, pp 185-187 (1992)).


---


<a id='Kaplan-Yorke-Dimension-1'></a>

### Kaplan-Yorke Dimension

<a id='ChaosTools.kaplanyorke_dim' href='#ChaosTools.kaplanyorke_dim'>#</a>
**`ChaosTools.kaplanyorke_dim`** &mdash; *Function*.



```
kaplanyorke_dim(lyapunovs::AbstractVector)
```

Calculate the Kaplan-Yorke dimension, a.k.a. Lyapunov dimension [1].

**Description**

The Kaplan-Yorke dimension is simply the point where `cumsum(lyapunovs)` becomes zero (interpolated). If the sum of the exponents never becomes negative the function will return the length of the input vector.

Useful in combination with [`lyapunovs`](lyapunovs.md#ChaosTools.lyapunovs).

**References**

[1] :  J. Kaplan & J. Yorke, *Chaotic behavior of multidimensional difference equations*, Lecture Notes in Mathematics vol. **730**, Springer (1979)


---


Notice that calling this function requires you to pass the Lyapunov exponents in an ordered vector form (largest to smallest). Example:


```julia
using DynamicalSystems
hen = Systems.henon()
D_kp = kaplanyorke_dim(lyapunovs(hen, 100000))
```

```
1.258711626829725
```

