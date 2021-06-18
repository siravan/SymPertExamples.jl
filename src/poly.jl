mutable struct Poly
    T::Type
    x
    n::Int
    a₀
    terms::Union{Dict{Int,Any},Nothing}

    Poly(val) = new(typeof(val), nothing, 0, val, nothing)
    Poly(T::Type, x) = new(T, value(x), 0, 0, nothing)
    Poly(T::Type, x, n, a₀, terms) = new(T, value(x), n, a₀, terms)
end

# a series of helper functions to facilitate generation of low-degree polynomials
function gen_poly(T::Type, x, a...)    
    p = Poly(T, x)
    
    for (i, val) in enumerate(a)
        p[i-1] = val 
    end
    
    update_order!(p)
    p
end 

gen_poly(x, a...) = gen_poly(Any, x, a...)

function degree(p::Poly)
    d = p.n
    while d > 0 && iszero(p[d])
        d -= 1
    end
    d
end

degree(p::Poly, x) = isequal(p.x, x) ? degree(p) : degree(sym(p), x)

"""
    poly_degree returns the degree of a polynomial equation eq
"""
function degree(poly, x)
    eqs = collect_powers(poly, x)
    round(Int, maximum(keys(eqs)))
end

terms(p::Poly) = p.terms
var(p::Poly) = p.x
leading(p::Poly) = degree(p) > 0 ? p[degree(p)] : p.a₀
cont(p::Poly) = (p.n > 0 ? gcd(p.a₀, values(p.terms)...) : p.a₀) * sign(leading(p))
prim(p::Poly) = p / cont(p)

function update_order!(p::Poly)
    p.n == 0 && return 0
    while iszero(p[p.n]) && p.n>0
        delete!(p.terms, p.n)
        p.n -= 1
    end
    if p.n == 0
        p.terms = nothing 
    end
    p
end

function poly(T::Type, eq, x)
    x === nothing && return Poly(T, nothing, 0, eq, nothing)
    terms = collect_powers(eq, x)
    terms = Dict{Int,Any}(convert(Int,k) => val for (k,val) in terms)
    
    if haskey(terms, 0)
        a₀ = terms[0]
        delete!(terms, 0)
    else
        a₀ = 0
    end
    n = isempty(terms) ? 0 : maximum(keys(terms))
    Poly(T, x, n, a₀, terms)
end

function poly(T::Type, eq)
    xs = get_variables(eq)

    if length(xs) > 1
        error("more than one implicit variable. Pass x to polymerize.")
    elseif length(xs) == 0
        return nothing
    else
        return poly(T, eq, xs[1])
    end
end

poly(eq, x) = poly(Any, eq, x)
poly(eq) = poly(Any, eq)
poly(T::Type, p::Poly) = convert(T, p)
poly(p::Poly) = p

polyas(p::Poly, val) = poly(p.T, val, p.x)

function Base.getindex(p::Poly, k::Integer)
    if k == 0
        return convert(p.T, p.a₀)
    elseif k > p.n || !haskey(p.terms, k)
        return convert(p.T, 0)
    else
        return convert(p.T, p.terms[k])
    end
end

function Base.setindex!(p::Poly, term, k::Number)
    if k == 0
        p.a₀ = term
    elseif p.terms !== nothing && iszero(term)
        delete!(p.terms, k)
        if k == p.n
            update_order!(p)
        end
    else
        if p.terms === nothing
            p.terms = Dict{Int,Any}()
        end
        p.terms[k] = term
        p.n = max(p.n, k)
    end
end

simplify_rational(x::Rational) = denominator(x) == 1 ? numerator(x) : x
simplify_rational(x) = x

sym(p::Poly) = (p.x === nothing ? p.a₀ : sum([simplify_rational(p[k])*p.x^k for k=0:p.n]))
sym(x) = x

function Base.iszero(p::Poly)
    !iszero(simplify(p.a₀)) && return false
    p.n > 0 && return all(iszero, values(p.terms))
    true
end

Base.show(io::IO, p::Poly) = print(io, sym(p))
Base.copy(p::Poly) = Poly(p.T, p.x, p.n, p.a₀, p.n>0 ? copy(p.terms) : nothing)

function extract_coef(p::Poly)
    c = p.T[p.a₀]
    if p.n > 0
        append!(c, values(p.terms))
    end
    c
end

###############################################################################

# pox (power-of-x) is a symbolic function to keep track of the powers of x
# pox(k,n) means k*x^n
@syms pox(k, n)

is_pox(x) = istree(x) && operation(x)==pox
is_not_pox(x) = !is_pox(x)

get_coef(p) = is_pox(p) ? arguments(p)[1] : p
get_power(p) = is_pox(p) ? arguments(p)[2] : 0

replace_x(eq, x) = substitute(eq, Dict(x => pox(1,1)))

iscomplex(x) = x isa Complex

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

