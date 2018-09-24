
<a id='Predefined-Systems-1'></a>

# Predefined Systems


Predefined systems exist in the `Systems` submodule exported by DynamicalSystemsBase.jl, in the form of functions that return a `DynamicalSystem`. They are accessed like:


```julia
using DynamicalSystems
ds = Systems.lorenz(ρ = 32.0)
ts = trajectory(ds, 10.0)
```


So far, the predefined systems that exist in the `Systems` sub-module are:

<a id='DynamicalSystemsBase.Systems.coupledstandardmaps' href='#DynamicalSystemsBase.Systems.coupledstandardmaps'>#</a>
**`DynamicalSystemsBase.Systems.coupledstandardmaps`** &mdash; *Function*.



```julia
coupledstandardmaps(M::Int, u0 = 0.001rand(2M); ks = ones(M), Γ = 1.0)
```

$$
\begin{aligned}
\theta_{i}' &= \theta_i + p_{i}' \\
p_{i}' &= p_i + k_i\sin(\theta_i) - \Gamma \left[
\sin(\theta_{i+1} - \theta_{i}) + \sin(\theta_{i-1} - \theta_{i})
\right]
\end{aligned}
$$

A discrete system of `M` nonlinearly coupled standard maps, first introduced in [1] to study diffusion and chaos thresholds. The *total* dimension of the system is `2M`. The maps are coupled through `Γ` and the `i`-th map has a nonlinear parameter `ks[i]`.

The first `M` entries of the state are the angles, the last `M` are the momenta.

[1] : H. Kantz & P. Grassberger, J. Phys. A **21**, pp 127–133 (1988)

<a id='DynamicalSystemsBase.Systems.double_pendulum' href='#DynamicalSystemsBase.Systems.double_pendulum'>#</a>
**`DynamicalSystemsBase.Systems.double_pendulum`** &mdash; *Function*.



```
double_pendulum(u0 = [π/2, 0, 0, rand()];
                G=10.0, L1 = 1.0, L2 = 1.0, M1 = 1.0, M2 = 1.0)
```

Famous chaotic double pendulum system (also used for our logo!). Keywords are gravity (G), lengths of each rod and mass of each ball (all assumed SI units).

The variables order is [θ1, dθ1/dt, θ2, dθ2/dt].

Jacobian is created automatically (thus methods that use the Jacobian will be slower)!

(please contribute the Jacobian and the e.o.m. in LaTeX :smile:)

The parameter container has the parameters in the same order as stated in this function's documentation string.

<a id='DynamicalSystemsBase.Systems.duffing' href='#DynamicalSystemsBase.Systems.duffing'>#</a>
**`DynamicalSystemsBase.Systems.duffing`** &mdash; *Function*.



```
duffing(u0 = [rand(), rand(), 0]; ω = 2.2, f = 27.0, d = 0.2, β = 1)
```

The (forced) duffing oscillator, that satisfies the equation

$$
\ddot{x} + d\cdot\dot{x} + β*x + x^3 = f\cos(\omega t)
$$

with `f, ω` the forcing strength and frequency and `d` the dampening.

The parameter container has the parameters in the same order as stated in this function's documentation string.

<a id='DynamicalSystemsBase.Systems.gissinger' href='#DynamicalSystemsBase.Systems.gissinger'>#</a>
**`DynamicalSystemsBase.Systems.gissinger`** &mdash; *Function*.



```julia
gissinger(u0 = 3rand(3); μ = 0.119, ν = 0.1, Γ = 0.9)
```

$$
\begin{aligned}
\dot{Q} &= \mu Q - VD \\
\dot{D} &= -\nu D + VQ \\
\dot{V} &= \Gamma -V + QD
\end{aligned}
$$

A continuous system that models chaotic reversals due to Gissinger [1], applied to study the reversals of the magnetic field of the Earth.

The parameter container has the parameters in the same order as stated in this function's documentation string.

[1] : C. Gissinger, Eur. Phys. J. B **85**, 4, pp 1-12 (2012)

<a id='DynamicalSystemsBase.Systems.henon' href='#DynamicalSystemsBase.Systems.henon'>#</a>
**`DynamicalSystemsBase.Systems.henon`** &mdash; *Function*.



```julia
henon(u0=zeros(2); a = 1.4, b = 0.3)
```

$$
\begin{aligned}
x_{n+1} &= 1 - ax^2_n+y_n \\
y_{n+1} & = bx_n
\end{aligned}
$$

The Hénon map is a two-dimensional mapping due to Hénon [1] that can display a strange attractor (at the default parameters). In addition, it also displays many other aspects of chaos, like period doubling or intermittency, for other parameters.

According to the author, it is a system displaying all the properties of the Lorentz system (1963) while being as simple as possible. Default values are the ones used in the original paper.

The parameter container has the parameters in the same order as stated in this function's documentation string.

[1] : M. Hénon, Commun.Math. Phys. **50**, pp 69 (1976)

<a id='DynamicalSystemsBase.Systems.henonheiles' href='#DynamicalSystemsBase.Systems.henonheiles'>#</a>
**`DynamicalSystemsBase.Systems.henonheiles`** &mdash; *Function*.



```
henonheiles(u0=[0, -0.25, 0.42081,0])
```

$$
\begin{aligned}
\dot{x} &= p_x \\
\dot{y} &= p_y \\
\dot{p}_x &= -x -2 xy \\
\dot{p}_y &= -y - (x^2 - y^2)
\end{aligned}
$$

The Hénon–Heiles system [1] was introduced as a simplification of the motion of a star around a galactic center. It was originally intended to study the existence of a "third integral of motion" (which would make this 4D system integrable). In that search, the authors encountered chaos, as the third integral existed for only but a few initial conditions.

The default initial condition is a typical chaotic orbit.

[1] : Hénon, M. & Heiles, C., The Astronomical Journal **69**, pp 73–79 (1964)

<a id='DynamicalSystemsBase.Systems.logistic' href='#DynamicalSystemsBase.Systems.logistic'>#</a>
**`DynamicalSystemsBase.Systems.logistic`** &mdash; *Function*.



```julia
logistic(x0 = rand(); r = 4.0)
```

$$
x_{n+1} = rx_n(1-x_n)
$$

The logistic map is an one dimensional unimodal mapping due to May [1] and is used by many as the archetypal example of how chaos can arise from very simple equations.

Originally intentend to be a discretized model of polulation dynamics, it is now famous for its bifurcation diagram, an immensly complex graph that that was shown be universal by Feigenbaum [2].

The parameter container has the parameters in the same order as stated in this function's documentation string.

[1] : R. M. May, Nature **261**, pp 459 (1976)

[2] : M. J. Feigenbaum, J. Stat. Phys. **19**, pp 25 (1978)

<a id='DynamicalSystemsBase.Systems.lorenz' href='#DynamicalSystemsBase.Systems.lorenz'>#</a>
**`DynamicalSystemsBase.Systems.lorenz`** &mdash; *Function*.



```julia
lorenz(u0=[0.0, 10.0, 0.0]; σ = 10.0, ρ = 28.0, β = 8/3) -> ds
```

$$
\begin{aligned}
\dot{X} &= \sigma(Y-X) \\
\dot{Y} &= -XZ + \rho X -Y \\
\dot{Z} &= XY - \beta Z
\end{aligned}
$$

The famous three dimensional system due to Lorenz [1], shown to exhibit so-called "deterministic nonperiodic flow". It was originally invented to study a simplified form of atmospheric convection.

Currently, it is most famous for its strange attractor (occuring at the default parameters), which resembles a butterfly. For the same reason it is also associated with the term "butterfly effect" (a term which Lorenz himself disliked) even though the effect applies generally to dynamical systems. Default values are the ones used in the original paper.

The parameter container has the parameters in the same order as stated in this function's documentation string.

[1] : E. N. Lorenz, J. atmos. Sci. **20**, pp 130 (1963)

<a id='DynamicalSystemsBase.Systems.lorenz96' href='#DynamicalSystemsBase.Systems.lorenz96'>#</a>
**`DynamicalSystemsBase.Systems.lorenz96`** &mdash; *Function*.



```
lorenz96(N::Int, u0 = rand(M); F=0.01)
```

$$
\frac{dx_i}{dt} = (x_{i+1}-x_{i-2})x_{i-1} - x_i + F
$$

`N` is the chain length, `F` the forcing. Jacobian is created automatically. (parameter container only contains `F`)

<a id='DynamicalSystemsBase.Systems.rikitake' href='#DynamicalSystemsBase.Systems.rikitake'>#</a>
**`DynamicalSystemsBase.Systems.rikitake`** &mdash; *Function*.



```julia
rikitake(u0 = [1, 0, 0.6]; μ = 1.0, α = 1.0)
```

$$
\begin{aligned}
\dot{x} &= -\mu x +yz \\
\dot{y} &= -\mu y +x(z-\alpha) \\
\dot{V} &= 1 - xz
\end{aligned}
$$

Rikitake's dynamo is a system that tries to model the magnetic reversal events by means of a double-disk dynamo system.

[1] : T. Rikitake Math. Proc. Camb. Phil. Soc. **54**, pp 89–105, (1958)

<a id='DynamicalSystemsBase.Systems.roessler' href='#DynamicalSystemsBase.Systems.roessler'>#</a>
**`DynamicalSystemsBase.Systems.roessler`** &mdash; *Function*.



```julia
roessler(u0=rand(3); a = 0.2, b = 0.2, c = 5.7)
```

$$
\begin{aligned}
\dot{x} &= -y-z \\
\dot{y} &= x+ay \\
\dot{z} &= b + z(x-c)
\end{aligned}
$$

This three-dimensional continuous system is due to Rössler [1]. It is a system that by design behaves similarly to the `lorenz` system and displays a (fractal) strange attractor. However, it is easier to analyze qualitatively, as for example the attractor is composed of a single manifold. Default values are the same as the original paper.

The parameter container has the parameters in the same order as stated in this function's documentation string.

[1] : O. E. Rössler, Phys. Lett. **57A**, pp 397 (1976)

<a id='DynamicalSystemsBase.Systems.shinriki' href='#DynamicalSystemsBase.Systems.shinriki'>#</a>
**`DynamicalSystemsBase.Systems.shinriki`** &mdash; *Function*.



```
shinriki(u0 = [-2, 0, 0.2]; R1 = 22.0)
```

Shinriki oscillator with all other parameters (besides `R1`) set to constants. *This is a stiff problem, be careful when choosing solvers and tolerances*.

<a id='DynamicalSystemsBase.Systems.standardmap' href='#DynamicalSystemsBase.Systems.standardmap'>#</a>
**`DynamicalSystemsBase.Systems.standardmap`** &mdash; *Function*.



```julia
standardmap(u0=0.001rand(2); k = 0.971635)
```

$$
\begin{aligned}
\theta_{n+1} &= \theta_n + p_{n+1} \\
p_{n+1} &= p_n + k\sin(\theta_n)
\end{aligned}
$$

The standard map (also known as Chirikov standard map) is a two dimensional, area-preserving chaotic mapping due to Chirikov [1]. It is one of the most studied chaotic systems and by far the most studied Hamiltonian (area-preserving) mapping.

The map corresponds to the  Poincaré's surface of section of the kicked rotor system. Changing the non-linearity parameter `k` transitions the system from completely periodic motion, to quasi-periodic, to local chaos (mixed phase-space) and finally to global chaos.

The default parameter `k` is the critical parameter where the golden-ratio torus is destroyed, as was calculated by Greene [2]. The e.o.m. considers the angle variable `θ` to be the first, and the angular momentum `p` to be the second, while both variables are always taken modulo 2π (the mapping is on the [0,2π)² torus).

The parameter container has the parameters in the same order as stated in this function's documentation string.

[1] : B. V. Chirikov, Preprint N. **267**, Institute of Nuclear Physics, Novosibirsk (1969)

[2] : J. M. Greene, J. Math. Phys. **20**, pp 1183 (1979)

<a id='DynamicalSystemsBase.Systems.towel' href='#DynamicalSystemsBase.Systems.towel'>#</a>
**`DynamicalSystemsBase.Systems.towel`** &mdash; *Function*.



```julia
towel(u0 = [0.085, -0.121, 0.075])
```

$$
\begin{aligned}
x_{n+1} &= a x_n (1-x_n) -0.05 (y_n +0.35) (1-2z_n) \\
y_{n+1} &= 0.1 \left( \left( y_n +0.35 \right)\left( 1+2z_n\right) -1 \right)
\left( 1 -1.9 x_n \right) \\
z_{n+1} &= 3.78 z_n (1-z_n) + b y_n
\end{aligned}
$$

The folded-towel map is a hyperchaotic mapping due to Rössler [1]. It is famous for being a mapping that has the smallest possible dimensions necessary for hyperchaos, having two positive and one negative Lyapunov exponent. The name comes from the fact that when plotted looks like a folded towel, in every projection.

Default values are the ones used in the original paper.

[1] : O. E. Rössler, Phys. Lett. **71A**, pp 155 (1979)

