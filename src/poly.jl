include("taylor.jl")

mutable struct Poly
    x
    n::Int
    a₀
    terms::Union{Dict,Nothing}

    Poly(x::SymbolicUtils.Sym) = new(x, 0, 0, Dict())
    Poly(val) = new(nothing, 0, val, nothing)
    Poly(p::Poly) = new(p.x, p.n, p.a₀, copy(p.terms))
    Poly(x, n, a₀, terms) = new(value(x), n, a₀, terms)
end

degree(p::Poly) = p.n
degree(x) = 0
terms(p::Poly) = p.terms
var(p::Poly) = p.x
leading(p::Poly) = p[degree(p)]

function update_order!(p::Poly)
    p.n = isempty(p.terms) ? 0 : maximum(keys(p.terms))
end

function poly(eq, x)
    terms = collect_powers(eq, x)
    if haskey(terms, 0)
        a₀ = terms[0]
        delete!(terms, 0)
    else
        a₀ = 0
    end
    n = isempty(terms) ? 0 : maximum(keys(terms))
    Poly(x, n, a₀, terms)
end

function poly(eq)
    xs = get_variables(eq)

    if length(xs) > 1
        error("more than one implicit variable. Pass x to polymerize.")
    elseif length(xs) == 0
        return nothing
    else
        return poly(eq, xs[1])
    end
end

poly(p::Poly) = p

rationalpoly(eq, x) = rationalize(poly(eq, x))
rationalpoly(eq) = rationalize(poly(eq))

function Base.getindex(p::Poly, k::Integer)
    if k == 0
        return p.a₀
    elseif k > p.n || !haskey(p.terms, k)
        return 0
    else
        return p.terms[k]
    end
end

function Base.setindex!(p::Poly, term, k::Number)
    if k == 0
        p.a₀ = term
    elseif isequal(term, 0)
        delete!(p.terms, k)
        if k == p.n
            update_order!(p)
        end
    else
        p.terms[k] = term
        p.n = max(p.n, k)
    end
end

sym(p::Poly) = p.a₀ + (p.n>0 ? sum([p[k]*p.x^k for k=1:p.n]) : 0)
sym(x) = x

Base.show(io::IO, p::Poly) = print(io, sym(p))
Base.iszero(p::Poly) = (p.n == 0) && iszero(p.a₀)

function extract_coef(p::Poly)
    c = [p.a₀]
    if p.n > 0
        append!(c, values(p.terms))
    end
    c
end

function is_int_poly!(p::Poly)
    b = all(isinteger, extract_coef(p))
    if b
        p.a₀ = convert(Int, p.a₀)
        for k = 1:p.n
            p[k] = convert(Int, p[k])
        end
    end
    b
end

function rationalize(p::Poly)
    try
        q = Poly(p.x, p.n, Rational(p.a₀), Dict())
        for k = 1:p.n
            q[k] = Rational(p[k])
        end
        return q
    catch
        return nothing
    end
end

###############################################################################

consensus(u, v) = (isequal(u.x,v.x) || v.x == nothing ? u.x : (u.x == nothing ? v.x : nothing))

function Base.:+(u::Poly, v::Poly)
    x = consensus(u, v)
    x == nothing && return Poly(sym(u) + sym(v))
    q = Poly(x)

    for k = 0:max(u.n, v.n)
        q[k] = u[k] + v[k]
    end
    q
end

Base.:+(u, v::Poly) = poly(u, v.x) + v
Base.:+(u::Poly, v) = u + poly(v, u.x)

function Base.:-(u::Poly, v::Poly)
    x = consensus(u, v)
    x == nothing && return Poly(sym(u) - sym(v))
    q = Poly(x)

    for k = 0:max(u.n, v.n)
        q[k] = u[k] - v[k]
    end
    q
end

Base.:-(u::Poly, v) = u - poly(v, u.x)
Base.:-(u, v::Poly) = poly(u, v.x) - v


function Base.:*(u::Poly, v::Poly)
    x = consensus(u, v)
    x == nothing && return Poly(sym(u) * sym(v))
    q = Poly(x)

    for i = 0:u.n
        for j = 0:v.n
            q[i+j] = u[i] * v[j]
        end
    end
    q
end

Base.:*(u::Poly, v) = u * poly(v, u.x)
Base.:*(u, v::Poly) = poly(u, v.x) * v

function divide_poly(u::Poly, v::Poly; op=/)
    x = consensus(u, v)
    x == nothing && return Poly(sym(u) / sym(v))
    q = Poly(x)
    r = Poly(u)
    m, n = u.n, v.n

    for i = m-n:-1:0
        q[i] = op(r[i+n], v[n])
        for j = n+i:-1:i
            r[j] = r[j] - q[i] * v[j-i]
        end
    end
    q, r
end

function Base.:/(u::Poly, v::Poly)
    q, r = divide_poly(u, v)
    q
end

Base.:/(u::Poly, v) = u / poly(v, u.x)
Base.:/(u, v::Poly) = poly(u, v.x) / v

function Base.:%(u::Poly, v::Poly)
    q, r = divide_poly(u, v)
    r
end

Base.:%(u::Poly, v) = u % poly(v, u.x)
Base.:%(u, v::Poly) = poly(u, v.x) % v

function pseudo_divide_poly(u::Poly, v::Poly)
    m, n = u.n, v.n
    α = v[n]^(m-n+1)
    q, r = divide_poly(u*α, v; op=÷)
    return q, r, α
end

##############################################################################

function gcd_naive(u::Poly, v::Poly)
    while degree(v) > 0
        u, v = v, u % v
    end
    return iszero(v) ? u : 1
end

function gcd_monic(u::Poly, v::Poly)
    u = u / leading(u)
    while degree(v) > 0
        v = v / leading(v)
        u, v = v, u % v
    end
    return iszero(v) ? u : 1
end

function split_poly(p::Poly)
    !is_int_poly!(p) && error("Only int polynomials can be split into cont/pp")
    cont = gcd(extract_coef(p)...)
    pp = p ÷ cont
    cont, pp
end
