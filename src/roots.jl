eval_for(p::Poly, x₀) = p.a₀ + (p.n > 0 ? sum([coef*x₀^k for (k,coef) in p.terms]) : 0)
eval_for(p::Poly, x, x₀) = isequal(p.x, x) ? eval_for(p, x₀) : eval_for(sym(p), x, x₀)
eval_for(p, x, x₀) = substitute(p, Dict(x => x₀))

(p::Poly)(x₀) = eval_for(p, x₀)

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
        xₙ₊₁ = xₙ - eval_for(f, x, xₙ) / eval_for(∂f, x, xₙ)
        
        if abs(xₙ₊₁ - xₙ) < abstol
            return xₙ₊₁
        else
            xₙ = xₙ₊₁
        end
    end
    return nothing
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
function find_roots(p, x, n=-1; abstol=1e-8)
    n = (n == -1 ? degree(p, x) : n)

    rs = Float64[]
    zs = Complex[]

    while n > 0
        z = solve_newton(p, x, Complex(rand(),rand()))
        if z != nothing  
            if abs(imag(z)) < abstol
                r = real(z)
                push!(rs, r)
                if p isa Poly
                    p = p / gen_poly(x, -r, 1)
                else
                    p = p / (x - r)
                end
                n -= 1
            else
                if abs(real(z)) < abstol
                    z = Complex(0, imag(z))
                end
                append!(zs, [z, conj(z)])
                #poly = poly / ((x-z)*(x-conj(z)))
                if p isa Poly
                    p = p / gen_poly(x, abs2(z), -2*real(z), 1)
                else
                    p = p / (x^2 - 2x*real(z) + abs2(z))
                end
                n -= 2
            end
        end
    end
    sort(rs), zs
end
