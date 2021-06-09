using Symbolics, SymbolicUtils
using SymbolicUtils: istree, operation, arguments
using Symbolics: value, get_variables, solve_for, derivative
using SymbolicUtils.Rewriters
using SymbolicUtils.Code

"""
    sym_taylor generates a Taylor series on variable x to order n and at
    point x₀, i.e., Σ aᵢ*(x - x₀)ⁱ
"""
function sym_taylor(x, n, x₀=0)
    x = value(x)
    @variables a[0:n]
    a = value.(a)
    a[1] + sum([a[i+1]*(x-x₀)^i for i = 1:n])
end

"""
    sym_frobenius generates a Frobenius series on variable x to order n and at
    point x₀, i.e., (x - x₀)^α Σ aᵢ*(x - x₀)ⁱ
"""
function sym_frobenius(x, n, x₀=0)
    x = value(x)
    @variables a[0:n]
    a = value.(a)
    @syms α
    (x-x₀)^α * sym_taylor(x, n, x₀)
end

# pox (power-of-x) is a symbolic function to keep track of the powers of x
# pox(k,n) means k*x^n
@syms pox(k, n)

is_pox(x) = istree(x) && operation(x)==pox
is_not_pox(x) = !is_pox(x)

get_coef(p) = is_pox(p) ? arguments(p)[1] : p
get_power(p) = is_pox(p) ? arguments(p)[2] : 0

replace_x(eq, x) = substitute(eq, Dict(x => pox(1,1)))

count_rule1 = @rule ^(pox(~k, ~n1), ~n2) => isequal(~k,1) ? pox(1, ~n1 * ~n2) : pox(^(~k,~n2), ~n1 * ~n2)
count_rule2 = @rule pox(~k1, ~n1) * pox(~k2, ~n2) => pox(~k1 * ~k2, ~n1 + ~n2)
count_rule3 = @acrule pox(~k, ~n) * ~u::is_not_pox => pox(~k * ~u, ~n)

"""
    collect_powers separates the powers of x in eq (a polynomial) and returns
    a dictionary of power => term
"""
function collect_powers(eq, x)
    eq = expand(expand_derivatives(eq))
    eq = replace_x(eq, x)
    #eq = Prewalk(PassThrough(count_rule1))(eq)
    eq = Fixpoint(Prewalk(PassThrough(Chain([count_rule1, count_rule2, count_rule3]))))(eq)

    if !istree(eq)
        return Dict{Any, Any}(0 => eq)
    elseif is_pox(eq)
        return Dict{Any, Any}(get_power(eq) => get_coef(eq))
    else
        eqs = Dict{Any, Any}()

        for term in arguments(eq)
            n = get_power(term)
            if haskey(eqs, n)
                eqs[n] = eqs[n] + get_coef(term)
            else
                eqs[n] = get_coef(term)
            end
        end

        return eqs
    end
end

"""
    strip_powers removes the powers from the output of collect_powers
    and returns a sorted list of terms/sub-equations
"""
function strip_powers(eqs)
    @syms α
    d = Dict(α => 0)
    eqs = Dict(substitute(power, d) => eq for (power, eq) in eqs)
    collect(values(sort(eqs)))
end

"""
    solve_rational gets a triangular system of linear equations and solve
    for the variables one by one starting from the first equation

    eqs is a list of terms as prepared by strip_powers

    solve_rational tries to generate rational results but falls back to reals
    if not possible
"""
function solve_rational(eqs)
    eqs = Num.(eqs)
    vals = Dict()

    for eq in eqs
        vars = get_variables(eq)
        vars = filter(x -> !haskey(vals,x), vars)

        if length(vars) > 0
            dep = vars[end]
            eq = substitute(eq, vals)
            # the last variable of each eq is considered the dependent variable
            A = value(substitute(eq, Dict(dep => 0)))
            B = value(substitute(eq, Dict(dep => 1)))
            try
                vals[dep] = A // (B - A)
            catch
                vals[dep] = A / (B - A)
            end
        end
    end

    vals
end

"""
    finalize_taylor substitutes the calculated coefficients back into the
    original Taylor or Frobenius series
"""
function finalize_taylor(x, y, vals)
    u = substitute(y, vals)
    vars = [v for v in get_variables(y) if !isequal(v,x)]
    factor(u, vars)
end

"""
    A main entry point
    solve_taylor solves a Taylor series expansion of a linear differential equation

    x is the independent variable (a Symbolic Sym)
    y is the dependent variable (a Symbolic Sym)
    diffeq is the differential equation, assuming diffeq ~ 0
"""
function solve_taylor(x, y, diffeq)
    eqs = strip_powers(collect_powers(diffeq, x))
    vals = solve_rational(eqs)
    finalize_taylor(x, y, vals)
end

"""
    A main entry point
    solve_frobenius solves a Frobenius series expansion of a linear differential equation

    x is the independent variable (a Symbolic Sym)
    y is the dependent variable (a Symbolic Sym)
    diffeq is the differential equation, assuming diffeq ~ 0
"""
function solve_frobenius(x, y, diffeq; abstol=1e-8, rationalize_coef=false)
    eqs = strip_powers(collect_powers(diffeq, x))
    eq = eqs[1]
    @syms α
    d = Dict(v => 1 for v in get_variables(eq) if !isequal(v, α))

    rs, zs = find_roots(substitute(eq, d), α)
    isempty(rs) && error("The indicial equation has no real root")
    α₁ = rs[end]
    @info "α is $α₁"

    d = Dict(α => α₁)
    eqs = [substitute(eq,d) for eq in eqs]
    vals = solve_rational(eqs[2:end])
    if rationalize_coef
        vals = rationalize(vals)
    end
    vals[α] = α₁
    finalize_taylor(x, y, vals)
end

###############################################################################

"""
    factor factors a polynomial eq over a list of variables (vars)
    it assumes trivial separatibilty
"""
function factor(eq, vars)
    f = []
    for v in vars
        d = Dict(u => (isequal(v,u) ? 1 : 0) for u in vars)
        push!(f, v * substitute(eq, d))
    end

    +(f...)
end

"""
    solve_newton is a symbolic Newton-Ralphson solver
    f is a symbolic equation to be solved (f ~ 0)
    x is the variable to solve
    x₀ is the initial guess
"""
function solve_newton(f, x, x₀; abstol=1e-8, maxiter=50)
    xₙ = Complex(x₀)
    ∂f = value(derivative(f, x))

    for i = 1:maxiter
        d = Dict(x => xₙ)
        xₙ₊₁ = xₙ - substitute(f, d) / substitute(∂f, d)
        if abs(xₙ₊₁ - xₙ) < abstol
            return xₙ₊₁
        else
            xₙ = xₙ₊₁
        end
    end
    return xₙ
end

"""
    poly_degree returns the degree of a polynomial equation eq
"""
function poly_degree(poly, x)
    eqs = collect_powers(poly, x)
    round(Int, maximum(keys(eqs)))
end

"""
    find_roots returns all the real and complex roots of a polynomial (poly)
    for variable x

    the output is rs, sz, where rs is a list of real roots and zs is a list of
    complex roots

    if poly is not a polynomial, the number of desired roots (n) needs to be
    provided. For example, find_roots(sin(x), x, 10) returns 10 different, but
    not necessarily sequential, multiples of π
"""
function find_roots(poly, x, n=-1; abstol=1e-8)
    n = (n == -1 ? poly_degree(poly, x) : n)

    rs = Float64[]
    zs = Complex[]

    while n > 0
        z = solve_newton(poly, x, Complex(rand(),rand()))
        if abs(imag(z)) < abstol
            r = real(z)
            push!(rs, r)
            poly = poly / (x - r)
            n -= 1
        else
            append!(zs, [z, conj(z)])
            poly = poly / ((x-z)*(x-conj(z)))
            n -= 2
        end
    end
    sort(rs), zs
end

##############################################################################

function continued_fraction(x, n; abstol=1e-8)
    try
        l = zeros(Int, n)
        for i = 1:n
            l[i] = floor(Int, x)
            Δx = x - l[i]
            Δx < abstol && break
            x = 1 / Δx
        end
        y = 0
        for i = n:-1:1
            if l[i] + y > 0
                y = 1 // (l[i] + y)
            end
        end
        return 1 // y
    catch
        return x
    end
end

function float_to_rational(x; n=10, abstol=1e-8)
    s = x < 0
    y = continued_fraction(abs(x), n) * (s ? -1 : +1)
    if abs(x - y) < abstol
        return y
    else
        return x
    end
end

# function rationalize(vals)
#     r = Dict()
#     for p in vals
#         eq = last(p)
#         vars = get_variables(eq)
#         if length(vars) == 1
#             var = vars[1]
#             eq = eq / var
#             if !(eq isa Int) || !(eq isa Rational)
#                 eq = float_to_rational(eq)
#             end
#             r[first(p)] = eq * var
#         else
#             r[first(p)] = eq
#         end
#     end
#     r
# end

################################### Examples! #################################

function test_harmonic()
    @syms x
    y = sym_taylor(x, 10)
    D = Differential(x)
    solve_taylor(x, y, D(D(y)) + y)
end

function test_airy()
    @syms x
    y = sym_taylor(x, 10)
    D = Differential(x)
    solve_taylor(x, y, D(D(y)) + x*y)
end

function test_bessel(ν = sqrt(2); rationalize_coef=false)
    @syms x
    y = sym_frobenius(x, 10)
    D = Differential(x)
    solve_frobenius(x, y, x^2 * D(D(y)) + x * D(y) - (x^2 + ν^2)*y; rationalize_coef=rationalize_coef)
end
