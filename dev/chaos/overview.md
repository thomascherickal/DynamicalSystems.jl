
<a id='Features-Overview-1'></a>

# Features Overview


The features offered in this documentation section come from the package [ChaosTools.jl](https://github.com/JuliaDynamics/ChaosTools.jl). If you are encountering an issue with some of the methods, you can report/open a new issue at the GitHub Issues page.


<a id='[Orbit-Diagrams](orbitdiagram)-1'></a>

### [Orbit Diagrams](orbitdiagram)


1. Orbit diagrams (aka bifurcation diagrams) of maps: [`orbitdiagram`](orbitdiagram.md#ChaosTools.orbitdiagram).
2. Poincaré surfaces of section for continuous systems: [`poincaresos`](orbitdiagram.md#ChaosTools.poincaresos).
3. Automated production of orbit diagrams for continuous systems: [`produce_orbitdiagram`](orbitdiagram.md#ChaosTools.produce_orbitdiagram).


<a id='[Lyapunov-Exponents](lyapunovs)-1'></a>

### [Lyapunov Exponents](lyapunovs)


The following treat systems where the equations of motion are known:


1. Maximum Lyapunov exponent for both discrete and continuous systems: [`lyapunov`](lyapunovs.md#ChaosTools.lyapunov).
2. Lyapunov *spectrum* for both discrete and continuous systems: [`lyapunovs`](lyapunovs.md#ChaosTools.lyapunovs).


<a id='[Entropies-and-Dimensions](entropies)-1'></a>

### [Entropies and Dimensions](entropies)


1. Generalized (Renyi) entropy: [`genentropy`](entropies.md#ChaosTools.genentropy).
2. Permutation entropy: [`permentropy`](entropies.md#ChaosTools.permentropy).
3. Fast and cheap (memory-wise) method for computing entropies of large datasets.
4. Generalized dimensions (e.g. capacity dimension, information dimension, etc.): [`generalized_dim`](entropies.md#ChaosTools.generalized_dim).
5. Kaplan-Yorke dimension: [`kaplanyorke_dim`](entropies.md#ChaosTools.kaplanyorke_dim).
6. Automated detection of best algorithmic parameters for calculating attractor dimensions.


And, in order to automatically deduce dimensions, we also offer methods for:


  * Partitioning a function $y(x)$ vs. $x$ into regions where it is approximated by a straight line, using a flexible algorithm with a lot of control over the outcome. See [`linear_regions`](entropies.md#ChaosTools.linear_regions).
  * Detection of largest linear region of a function $y(x)$ vs. $x$ and extraction of the slope of this region.


<a id='[Nonlinear-Timeseries-Analysis](nlts)-1'></a>

### [Nonlinear Timeseries Analysis](nlts)


1. Methods for estimating good [`reconstruct`](../definition/reconstruction.md#DynamicalSystemsBase.reconstruct) parameters.
2. Broomhead-King coordinates: [`broomhead_king`](nlts.md#ChaosTools.broomhead_king).
3. Numerically determining the maximum Lyapunov exponent of a (e.g. experimentally) measured timeseries: [`numericallyapunov`](nlts.md#ChaosTools.numericallyapunov).


<a id='[Periodicity](periodicity)-1'></a>

### [Periodicity](periodicity)


1. Numerical method to find unstable and stable fixed points of *any order* $n$ of a discrete map (of any dimensionality): [`periodicorbits`](periodicity.md#ChaosTools.periodicorbits).

      * Convenience functions for defining and realizing all possible combinations of $\mathbf{\Lambda}_k$ matrices required in the above method.


<a id='[Chaos-Detection](chaos_detection)-1'></a>

### [Chaos Detection](chaos_detection)


1. The Generalized Alignment Index: $\text{GALI}_k$ : [`gali`](chaos_detection.md#ChaosTools.gali).

      * Implemented for both discrete and continuous systems.

