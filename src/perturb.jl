using Symbolics, SymbolicUtils

@variables x y z ϵ a[1:4]

def_taylor(x, ps) = sum([a*x^i for (i,a) in enumerate(ps)])
def_taylor(x, ps, p₀) = p₀ + def_taylor(x, ps)
expand_sin(x, n) = sum([(isodd(k) ? -1 : 1)*(-x)^(2k-1)/factorial(2k-1) for k=1:n])

function collect_powers(eq, x, ns; max_power=100)
    eq = substitute(expand(eq), Dict(x^j => 0 for j=last(ns)+1:max_power))

    eqs = []
    for i in ns
        powers = Dict(x^j => (i==j ? 1 : 0) for j=1:last(ns))
        push!(eqs, substitute(expand(eq), powers))
    end
    eqs
end

function solve_coef(eqs, ps)
    vals = Dict()

    for i = 1:length(ps)
        eq = substitute(eqs[i], vals)
        vals[ps[i]] = Symbolics.solve_for(eq ~ 0, ps[i])
    end
    vals
end

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

"""
    x^5 + x = 1
"""
function test_quintic(n=4)
    @variables ϵ z a[1:n]
    x = def_taylor(ϵ, a, 1)
    y = x^5 + ϵ*x - 1
    eqs = collect_powers(y, ϵ, 1:n)
    vals = solve_coef(eqs, a)
    sol = substitute(x, vals)

    xₚ = substitute(sol, Dict(ϵ => 1))
    xₙ = solve_newton(z^5 + z - 1, z, 1.0)
    return xₚ, xₙ
end

"""
    E - e * sin(E) = M
"""
function test_kepler(n=4)
    @variables ϵ z M a[1:n]
    x = def_taylor(ϵ, a, M)
    y = x - ϵ * expand_sin(x, n) - M
    eqs = collect_powers(y, ϵ, 1:n)
    vals = solve_coef(eqs, a)
    sol = substitute(x, vals)
    println(sol)

    𝑒 = 0.01671    # The Earth eccentricity
    M₀ = π/2

    xₚ = substitute(sol, Dict(ϵ => 𝑒, M => M₀))
    xₙ = solve_newton(z - 𝑒*sin(z) - M₀, z, M₀)
    return xₚ, xₙ
end

function test_rocket(n=3)
    @variables ϵ t y[1:n](t) ∂∂y[1:n]
    x = def_taylor(ϵ, y)
    ∂∂x = def_taylor(ϵ, ∂∂y)
    eq = ∂∂x * (1 + ϵ*x)^2 + 1
    eqs = collect_powers(eq, ϵ, 0:n-1)
    vals = solve_coef(eqs, ∂∂y)

    D = Differential(t)
    subs = Dict(∂∂y[i] => D(D(y[i])) for i = 1:n)
    eqs = [substitute(first(v), subs) ~ substitute(last(v), subs) for v in vals]

    sys = ODESystem(eqs, t)
    sys = ode_order_lowering(sys)

    sys
end

function test_oscillator(n=3)
    @variables ϵ t y[1:n](t) ∂y[1:n] ∂∂y[1:n]
    x = def_taylor(ϵ, y)
    ∂x = def_taylor(ϵ, ∂y)
    in∂∂x = def_taylor(ϵ, ∂∂y)
    eq = ∂∂x + 2*ϵ*∂x + x
    eqs = collect_powers(eq, ϵ, 1:n)
    vals = solve_coef(eqs, ∂∂y)

    D = Differential(t)
    subs1 = Dict(∂y[i] => D(y[i]) for i = 1:n)
    subs2 = Dict(∂∂y[i] => D(D(y[i])) for i = 1:n)
    subs = subs1 ∪ subs2
    eqs = [substitute(first(v), subs) ~ substitute(last(v), subs) for v in vals]

    sys = ODESystem(eqs, t)
    sys = ode_order_lowering(sys)

    sys
end
