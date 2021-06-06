# Mixed Symbolic/Numerical Methods for Perturbation Theory

[**Symbolics.jl**](https://github.com/JuliaSymbolics/Symbolics.jl) is a fast and modern Computer Algebra System (CAS) written in the Julia Programming Language. It is an integral part of the [SciML](https://sciml.ai/) ecosystem of differential equation solvers and scientific machine learning packages. While **Symbolics.jl** is primarily designed for modern scientific computing (e.g., auto-differentiation, machine learning), it is a powerful CAS and can also be useful for *classic* scientific computing. One such application is using the *perturbation* theory to solve algebraic and differential equations.

Perturbation methods are a collection of techniques to solve intractable problems that generally don't have a closed solution but depend on a tunable parameter and have closed or easy solutions for some values of the parameter. The main idea is to assume a solution as a power series in the tunable parameter (say ϵ), such that ϵ = 0 corresponds to an easy solution. The hallmark of the perturbation method is the generation of long and convoluted intermediate equations, which are subjected to algorithmic and mechanical manipulations. Therefore, these problems are well suited for CAS. Specifically, the tight coupling between **Symbolics.jl** and **DifferentialEquations.jl** is very useful.

The following tutorials discuss the general steps of using the perturbation theory to solve algebraic and differential equations.

- Perturbation Theory
  - [Mixed Symbolic/Numerical Methods for Perturbation Theory - Algebraic Equations](http://svtsim.com/perturbation/01-perturbation_algebraic.html)
  - [Mixed Symbolic/Numerical Methods for Perturbation Theory - Differential Equations](http://svtsim.com/perturbation/02-perturbation_differential.html)

The code discussed in the tutorials is also available (with some modifications) in `src/perturb.jl`. *Mixed Symbolic/Numerical Methods for Perturbation Theory - Algebraic Equations* examples are `test_quintic` and `test_kepler`, and *Mixed Symbolic/Numerical Methods for Perturbation Theory - Differential Equations* examples are `test_rocket` and `test_oscillator`.

# Power Series (Taylor and Frobenius) Solution for Solving Ordinary Differential Equations

As a related subject to the perturbation method, a standard technique to solve Ordinary Differential Equations (ODE) that have no close solution is to use power series approximations. Generally, Taylor series are used to find the solutions in the neighborhood of a regular (non-singular) point and Frobenius series are used around a regular singular point.

`src/taylor.jl` includes a collection of utility functions to help create, substitute and solve linear ODEs using the Taylor and Frobenius methods.

You can start using these functions by `include("src/taylor.jl")`.

### Harmonic Equation

The harmonic equation (𝑦" + 𝑦 = 0) is a simple and regular second-order linear ODE with a general solution 𝑐₁ sin(𝑥) + 𝑐₂ cos(𝑥). We can find the Taylor series as

```julia
@syms x                   # define a symbolic independent variable
y = sym_taylor(x, 10)     # define an order 10 Taylor series in x
D = Differential(x)       # define the Differentiation operator
solve_taylor(x, y, D(D(y)) + y) # define the problem and solve!
```

The output, which reproduces the power series of *sin* and *cos*, is

```julia
a₀*(1 + (1//2)*(x^2) + (1//24)*(x^4) + (1//720)*(x^6) + (1//40320)*(x^8) + (1//3628800)*(x^10)) +
a₁*(x + (1//6)*(x^3) + (1//120)*(x^5) + (1//5040)*(x^7) + (1//362880)*(x^9))
```

Note that a₀ and a₁ are arbitrary constants.

### Airy Equation

The next example is still regular. We try the Airy equation (𝑦" + 𝑥𝑦 = 0).

```julia
@syms x                   # define a symbolic independent variable
y = sym_taylor(x, 10)     # define an order 10 Taylor series in x
D = Differential(x)       # define the Differentiation operator
solve_taylor(x, y, D(D(y)) + x*y) # define the problem and solve!
```

The solution, which can be recast as Airy functions, is

```julia
a₀*(1 + (1//6)*(x^3) + (1//180)*(x^6) + (1//12960)*(x^9)) +
a₁*(x + (1//12)*(x^4) + (1//504)*(x^7) + (1//45360)*(x^10))
```

### Bessel Equation

For the next example, we have chosen a regular singular problem, the Bessel equation:

𝑥² 𝑦" + 𝑥 𝑦' - (𝑥² + ν²) 𝑦 = 0 .

This equation has a singularity at 𝑥 = 0. For now, we need to replace the parameter ν with a numerical value. We use ν = √2 in this example. First, let's try applying the Taylor method to this problem.

```julia
@syms x                   # define a symbolic independent variable
y = sym_taylor(x, 10)     # define an order 10 Taylor series in x
D = Differential(x)       # define the Differentiation operator
solve_taylor(x, y, x^2 * D(D(y)) + x * D(y) - (x^2 + 2)*y)
```

The result is

```
0 // 1
```

This is technically correct (𝑦 = 0 is always a solution to a homogeneous linear differential equation) but is totally useless! We need to use the Frobenius method for this problem.

```julia
@syms x
y = sym_frobenius(x, 10)
D = Differential(x)
solve_frobenius(x, y, x^2 * D(D(y)) + x * D(y) - (x^2 + 2)*y)
```

Now we get a useful solution

```julia
a₀*(1 + 0.003791260736238831(x^4) + 8.262161907723593e-7(x^8) - (0.10355339059327377(x^2)) - (7.15729744885111e-5(x^6)) - (6.440510459607151e-9(x^10)))
*(x^1.4142135623730951)
```

Note that the series is multiplied by 𝑥^(√2), as expected for a Bessel function. In fact, the resulting
series is s Bessel function. 
