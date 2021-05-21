# Mixed Symbolic/Numerical Methods for Perturbation Theory

[**Symbolics.jl**](https://github.com/JuliaSymbolics/Symbolics.jl) is a fast and modern Computer Algebra System (CAS) written in the Julia Programming Language. It is an integral part of the [SciML](https://sciml.ai/) ecosystem of differential equation solvers and scientific machine learning packages. While **Symbolics.jl** is primarily designed for modern scientific computing (e.g., auto-differentiation, machine learning), it is a powerful CAS and can also be useful for *classic* scientific computing. One such application is using the *perturbation* theory to solve algebraic and differential equations.

Perturbation methods are a collection of techniques to solve intractable problems that generally don't have a closed solution but depend on a tunable parameter and have closed or easy solutions for some values of the parameter. The main idea is to assume a solution as a power series in the tunable parameter (say ϵ), such that ϵ = 0 corresponds to an easy solution. The hallmark of the perturbation method is the generation of long and convoluted intermediate equations, which are subjected to algorithmic and mechanical manipulations. Therefore, these problems are well suited for CAS. Specifically, the tight coupling between **Symbolics.jl** and **DifferentialEquations.jl** is very useful.

The follwoing tutorials discuss the general steps of using the perturbation theory to solve algebraic and differential equations.

- Perturbation Theory
  - [Mixed Symbolic/Numerical Methods for Perturbation Theory - Algebraic Equations](http://svtsim.com/perturbation/01-perturbation_algebraic.html)
  - [Mixed Symbolic/Numerical Methods for Perturbation Theory - Differential Equations](http://svtsim.com/perturbation/02-perturbation_differential.html)

The code discussed in the tutorials is also available (with some modifications) in `src/perturb.jl`. *Mixed Symbolic/Numerical Methods for Perturbation Theory - Algebraic Equations* examples are `test_quintic` and `test_kepler`, and *Mixed Symbolic/Numerical Methods for Perturbation Theory - Differential Equations* examples are `test_rocket` and `test_oscillator`.
