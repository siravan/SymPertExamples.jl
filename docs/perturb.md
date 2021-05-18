# Using SciML Symbolics to Solve Perturbation Problems

## Background

[**Symbolics.jl**](https://github.com/JuliaSymbolics/Symbolics.jl) is a fast and modern Computer Algebra System (CAS) written in Julia Programming Language. It is part of the [SciML](https://sciml.ai/) ecosystem of differential equation solvers and scientific machine learning packages. While **Symbolics.jl** is primarily designed for modern scientific computing (e.g., machine learning), it is a powerful CAS and can be useful in *classic* scientific computing, like *perturbation* problems.

Perturbation methods are a collection of techniques to solve algebraic and differential equations. The target problems generally don't have a closed solution. However, they depend on a tunable parameter and have closed-form or easy solutions for some values of the parameter. The main idea is to assume a solution as a power series in the tunable parameter (say $\epsilon$), such that $\epsilon = 0$ corresponds to a closed solution.

We will discuss the general steps of the perturbation methods in four examples below. One hallmark of the perturbation method is the generation of long and involved intermediate equations, which are subjected to algorithmic and mechanical manipulations. Therefore, these problems are well suited for CAS. In fact, CAS softwares have been used to help with the perturbation calculations since the 1950s.

In this tutorial our goal is to show how to use Julia and **Symbolics.jl** to solve simple perturbation problems.

## Solving the Quintic

We start with the "hello world!" analog of the perturbation problems: solving the quintic (fifth-order) equations. We want to find $x$ such that $x^5 + x = 1$. According to the Abel's theorem, a general quintic equation does not have a closed form solution. Of course, we can easily solve this equation numerically; for example, using the Newton's method. Here, we use the following implementation of the Newton's method:

```Julia
function solve_newton(f, x, x₀; abstol=1e-8, maxiter=50)
    xₙ = Float64(x₀)
    fₙ₊₁ = x - f / Symbolics.derivative(f, x)

    for i = 1:maxiter
        xₙ₊₁ = substitute(fₙ₊₁, Dict(x => xₙ))
        if abs(xₙ₊₁ - xₙ) < abstol
            return xₙ₊₁
        else
            xₙ = xₙ₊₁
        end
    end
    return xₙ₊₁
end
```

In this code, `Symbolics.derivative(eq, x)` does exactly what it names implies: it calculates the symbolic derivative of `eq` (a **Symbolics.jl** expression) with respect to `x` (a **Symbolics.jl** variable). We use `Symbolics.substitute(eq, D)` to evaluate the update formula by substituting variables or sub-expressions (defined as a dictionary `D`) in `eq`. It should be noted that `substitute` is the workhorse of our code and will be used multiple times in the rest of this tutorial. `solve_newton` is written with simplicity and clarity in mind and not performance.

Let's go back to our quintic. We can define a Symbolics variable as `@variables x` and then solve the equation as `solve_newton(x^5 + x - 1, x, 1.0)` (here, `x₀ = 0` is our first guess). The answer is `x = 0.7549`. Now, let's see how we can solve this problem using the perturbation method.

We introduce a tuning parameter $\epsilon$ into our equation: $x^5 + \epsilon x = 1$. If $\epsilon = 1$, we get our original problem. For $\epsilon = 0$, the problem transforms to an easy one: $x^5 = 1$ which has a solution $x = 1$ (and four complex solutions which we ignore here). We expand $x$ as a power series on $\epsilon$:


<img src="https://render.githubusercontent.com/render/math?math=x(\epsilon) = a_0 + a_1 \epsilon + a_2 \epsilon^2 + O(\epsilon^3)">

$a_0$ is the solution of the easy equation, therefore $a_0 = 1$.

Substituting in the original problem,

$$
  (1 + a_1 \epsilon + a_2 \epsilon^2)^5 + \epsilon (1 + a_1 \epsilon + x_a \epsilon^2) - 1 = 0
$$

Expanding the equations, we get

$$
  \epsilon (1 + 5a_1) + \epsilon^2 (a_1 + 5 a_2 + 10 a_1^2) + O(\epsilon^3) = 0
$$

This equation should hold for each power of $\epsilon$,

$$
  1 + 5a_1 = 0
  \,,
$$

and

$$
  a_1 + 5 a_2 + 10 a_1^2 = 0
  \,.
$$

We solve the first equation to get $a_1 = -\frac{1}{5}$. Substituting in the second one and solve for $a_2$:

$$
  a_2 = -\frac{ -\frac{1}{5} + 10 (-\frac{1}{5})^2 }{5} = -\frac{1}{25}
$$

Finally,

$$
  x(\epsilon) = 1 -\frac{1}{5} \epsilon - \frac{1}{25}\epsilon^2 + + O(\epsilon^3)
$$

Solving the original problem, $x(1) = 0.76$, compared to 0.75487767 calculated numerically. We can improve the accuracy by including more terms in the expansion of $x$. However, the calculations, while straightforward, become messy and intractable for do manually very quickly. This is why a CAS is very helpful to solve perturbation problems. So, let's see how we can do these calculations in Julia (`test_quintic` function).

Let `n = 2` be the order of the expansion. We start by defining the variables:

```Julia
  @variables @variables ϵ a[1:n]        
```

Then, we define `x = 1 + a[1]*ϵ + a[2]*ϵ^2`. Note that in `test_quintic` we use the helper function `def_taylor` to define `x` by calling it as `x = def_taylor(ϵ, a, 1)`. The next step is to substitute $x$ is the problem `y = x^5 + ϵ*x - 1`. Now, `y` is

```Julia
  ϵ*(1 + a₁*ϵ + a₂*(ϵ^2)) + (1 + a₁*ϵ + a₂*(ϵ^2))^5 - 1
```

Or in the expanded form (calculated as `expand(y)`):

```Julia
ϵ + a₁*(ϵ^2) + a₂*(ϵ^3) + (a₁^5)*(ϵ^5) + (a₂^5)*(ϵ^10) + 5a₁*ϵ + 5a₂*(ϵ^2) +
10(a₁^2)*(ϵ^2) + 10(a₁^3)*(ϵ^3) + 5(a₁^4)*(ϵ^4) + 10(a₂^2)*(ϵ^4) +
10(a₂^3)*(ϵ^6) + 5(a₂^4)*(ϵ^8) + 20a₁*a₂*(ϵ^3) + 30a₁*(a₂^2)*(ϵ^5) +
20a₁*(a₂^3)*(ϵ^7) + 5a₁*(a₂^4)*(ϵ^9) + 30a₂*(a₁^2)*(ϵ^4) + 20a₂*(a₁^3)*(ϵ^5) +
5a₂*(a₁^4)*(ϵ^6) + 30(a₁^2)*(a₂^2)*(ϵ^6) + 10(a₁^2)*(a₂^3)*(ϵ^8) +
10(a₁^3)*(a₂^2)*(ϵ^7)
```

We need a way to get the coefficients of different powers of $\epsilon$. Function `collect_powers(eq, x, ns)` returns the powers of `x` in expression `eq`. Argument `ns` is the range of the powers.

```Julia
function collect_powers(eq, x, ns; max_power=100)
    eq = substitute(expand(eq), Dict(x^j => 0 for j=last(ns)+1:max_power))

    eqs = []
    for i in ns
        powers = Dict(x^j => (i==j ? 1 : 0) for j=1:last(ns))
        push!(eqs, substitute(expand(eq), powers))
    end
    eqs
end
```

For example, `collect_powers(y, ϵ, 1:2)` returns `eqs = [1 + 5a₁, a₁ + 5a₂ + 10(a₁^2)]`. `collect_powers` uses `substitute` to find the coefficient of a given power of `x` by passing a `Dict` with all powers of `x` set to 0, except the target power which is set to 1. To find the coefficient of `ϵ^2` in `y`, we can write

```julia
  substitute(expand(y), Dict(
    ϵ => 0,
    ϵ^2 => 1,
    ϵ^3 => 0,
    ϵ^4 => 0,
    ϵ^5 => 0,
    ϵ^6 => 0,
    ϵ^7 => 0,
    ϵ^8 => 0)
  )
```

The next step is find the coefficients of `ϵ` in the expansion of `x`. **Symbolics.jl** has a function `Symbolics.solve_for` that can solve systems of linear equations. The system described by `eqs` does not seem linear (note `10(a₁^2)` in `eqs[2]`), but upon closer inspection is found to be in fact linear (this is a feature of the permutation method). We can start by solving `eqs[1]` for `a₁` and then substitute it in `eqs[2]` and solve for `a₂`.  This process is done by function `solve_coef(eqs, ps)`:

```julia
function solve_coef(eqs, ps)
    vals = Dict()

    for i = 1:length(ps)
        eq = substitute(eqs[i], vals)
        vals[ps[i]] = Symbolics.solve_for(eq ~ 0, ps[i])
    end
    vals
end
```

Here, `eqs` is an array of expressions (assumed to be equal to 0) and `ps` is an array of variables. The result is a dictionary of variable => value pairs. For example, `solve_coef(eqs, a)` is `Dict(a₁ => -0.2,
a₂ => -0.04)`. We can use larger values of `n` to improve the accuracy of estimations:

| n | x              |
|---|----------------|
|1  |0.8 |
|2  |0.76|
|3  |0.752|
|4  |0.752|
|5  |0.7533|
|6  |0.7543|
|7  |0.7548|
|8  |0.7550|

Remember the numerical value is 0.7549.

The two functions `collect_powers` and `solve_coef(eqs, a)` are used in all the examples in this tutorial and show the main steps in the permutation method.

## Solving the Kepler's Equation

Historically, the perturbation methods were first invented to solve orbital calculations needed to calculate the orbit of the moon and planets. In homage to this history, our second example has a celestial theme. Our goal is solve the Kepler's equation:

$$
  E - e \sin(E) = M
  \,,
$$

where $e$ is the *eccentricity* of the elliptical orbit, $M$ is the *mean anomaly*, and $E$ (unknown) is the *eccentric anomaly* (the angle between the position of a planet in an elliptical orbit and the point of periapsis). We want to find a function $E(M; e)$. As the first example, it is easy to solve this problem using the Newton's method. For example, let $e = 0.01671$ (the eccentricity of the Earth) and $M=\pi / 2$. We have `solve_newton(x - e*sin(x) - M, x, M)` equals to 1.5875 (compared to $\pi/2 = 1.5708$). Now, we try to solve the same problem using the perturbation techniques (see function `test_kepler`.

For $e = 0$, $E = M$. Therefore, we can use $e$ as our perturbation parameter. For consistency, we rename it to $\epsilon$. We start by defining the variables and $x$ (assuming `n = 3`):

```julia
  @variables ϵ M a[1:n]
  x = def_taylor(ϵ, n, M)  
```

The problem equation is `y = E - \epsilon * sin(E) - M`. We further simplify by substituting $\sin$ with its power series (using `expand_sin` helper function):

$$
  \sin(E) = x - \frac{x^3}{6} + \frac{x^5}{120} - \frac{x^7}{5040} + O(x^9)
  \,.
$$

We follow the same algorithm as before. We collect the coefficients of the powers of $\epsilon$:

```
  eqs = collect_powers(y, ϵ, 1:n)
```

and then solve for `a`:

```
  vals = solve_coef(eqs, a)
```

Finally, we substitute `vals` back in `x` to get a formula to calculate `E`:

```
  sol = substitute(x, vals)
  substitute(sol, Dict(ϵ => 0.01671, M => π/2))
```

The result is 1.5876, compared to the numerical value of 1.5875. It is customary to order `sol` based on the powers of `M` instead of `ϵ`. We can calculate this series as `collect_powers(sol, M, 0:3)
`. The result (after cleanup) is

```julia
  E(M, ϵ) =
    (1 + ϵ + ϵ^2 + ϵ^3)*M
    - (ϵ + 4*ϵ^2 + 10*ϵ^3)*M^3/6
    + (ϵ + 16*ϵ^2 + 91*ϵ^3)*M^5/120
```

Comparing the formula to the one for $E$ in the [Wikipedia article on the Kepler's equation](https://en.wikipedia.org/wiki/Kepler%27s_equation):

$$
  E = \frac{1}{1-\epsilon}M
  - \frac{\epsilon}{(1-\epsilon)^4} \frac{M^3}{3!}
  + \frac{(9\epsilon^2 + \epsilon)}{(1-\epsilon)^7} \frac{M^5}{5!}
  + \cdots
$$

The first deviation is in the coefficient of $\epsilon^3 M^5$.
