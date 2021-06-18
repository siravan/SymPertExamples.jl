using Primes

struct FactoredPoly
    factors::Array{Pair{Poly,Int}}

    FactoredPoly() = new(Array{Pair{Poly,Int}}[])
end

function add_factor!(f::FactoredPoly, p, power=1)
    push!(f.factors, p => power)
end

sym(f::FactoredPoly) = prod(map(v -> sym(first(v))^last(v), f.factors); init=1)
poly(f::FactoredPoly) = prod(map(v -> first(v)^last(v), f.factors); init=1)

Base.show(io::IO, f::FactoredPoly) = print(io, sym(f))

###########################################################################################

function split_poly(p::Poly)
    l, p = integer_poly(p)
    cont = gcd(extract_coef(p)...)
    pp = p ÷ cont
    cont * l, pp
end

# Lagrange interpolation algorithm for rational polynomials
function interp_poly(x, xs, ys)
    length(xs) != length(ys) && error("interp_poly needs the same number of x and y values")
    n = length(xs)
    Ps = 0
    for i = 1:n        
        A = *([(x-xs[j]) for j=1:n if j!=i]...)
        B = *([(xs[i]-xs[j]) for j=1:n if j!=i]...)        
        Ps += ys[i] *  A // B
    end
    return poly(Rational, Ps, x)
end

function pre_factor(p::Poly)
    l, p = integer_poly(p)
    x = p.x
    r = copy(p)
    f = Poly[]

    a₀ = p[0]
    aₙ = p[degree(p)]
    
    for root in find_roots(p, x)[1]
        ρ = Rational(root)
        println(ρ)
        if a₀ % numerator(ρ) == 0 && aₙ % denominator(ρ) == 0
            q = gen_poly(x, -numerator(ρ), denominator(ρ))                     
            r = r / q
            push!(f, q)
        end
    end 
    
    return r, f, l
end

function list_divisors(x)
    f = factor(x)
    l = [1]
    
    for (p, k) in f
        r = Int[]
        for i = 0:k
            append!(r, p^i .* l)
        end
        l = r
    end
    l
end

# note k is 1-based
function generate_cross_lists(D)     
    n = *(length.(D)...)
    L = zeros(Int, (length(D), n))

    for k = 1:n
        j = k - 1       
        for i = 1:length(D)
            n = length(D[i])
            L[i,k] = D[i][(j % n) + 1]
            j = j ÷ n
        end
    end
    return L
end

function assemble_factors(f, r, l)
    if !iszero(r - 1)
        # push!(f, r)
        add_factor!(f, r)
    end
    if !isone(l)
        # push!(f, Poly(l))
        add_factor!(f, Poly(l))
    end
    f
end

function factor_schubert_kronecker(p::Poly)
    l, p = integer_poly(p)
    x = p.x
    r = copy(p)
    # f = Poly[]
    f = FactoredPoly()
    A = Int[]
    D = Array{Int}[]

    i = 0
    while length(D) <= degree(r)÷2
        a = ((i+1)÷2) * (isodd(i) ? -1 : +1)
        u = convert(Int, p(a))
        
        if iszero(u)            
            while iszero(r % (x-a))
                r = r / (x-a)
                # push!(f, poly(x-a, x))
                add_factor!(f, poly(x-a, x))
            end
        else
            push!(A, a)
            divs = list_divisors(abs(u))            
            push!(D, [divs; -divs]) 
        end
        i += 1        
    end
    
    for d = 1:degree(r)÷2 
        L = generate_cross_lists(D[1:d+1])
        xs = A[1:d+1]
        for i = 1:size(L,2)
            ys = L[:,i]
            q = interp_poly(x, xs, ys)
            
            while degree(q) > 0 && iszero(r % q)
                r = r / q
                # push!(f, q)
                add_factor!(f, q)
                if degree(r) < 2d
                    return assemble_factors(f, r, l)                    
                end
            end
        end
        if degree(r) < 2*(d+1)
            return assemble_factors(f, r, l)
        end
    end    
    assemble_factors(f, r, l)
end

"""
    square-free decomposition of p
    Yun's algorithm
"""
function decompose(p::Poly)
    p = rationalize(p)
    f = FactoredPoly()
    p′ = derivative(p)
    g = gcd(p, p′)    
    r = p / g
    i = 1
    c = 1
    
    while !iszero(g - 1) && r isa Poly        
        s = gcd(g, r)
        p = r / s
        if !iszero(p - 1)
            c *= leading(p)
            #push!(f, (p / leading(p), i))
            add_factor!(f, p / leading(p), i)
        end
        i += 1
        r = s
        g = g / s    
    end

    if !iszero(r - 1)
        c *= leading(r)
        # push!(f, (r / leading(r), i))
        add_factor!(f, r / leading(r), i)        
    end

    if c != 1
        push!(f, (Poly(c), 0))
    end

    return f
end

function recompose(f)
    y = Poly(1)
    
    for (k, p) in f
        y = y * p^k
    end
    y
end

function Primes.factor(p::Poly; verbose=false)
    _, pᵢ = integer_poly(p)
    c, pₘ = to_monic(pᵢ)
    f = decompose(pₘ)
    verbose && printstyled('\t', f, '\n'; color=:red)
    h = FactoredPoly()
    p₁ = poly(Rational, 1, p.x)

    for w in f.factors
        v, k = first(w), last(w)
        if degree(v) > 0
            for w in factor_schubert_kronecker(v).factors
                w₁ = prim(from_monic(first(w), c))
                # push!(h, (w₁, k))                
                add_factor!(h, w₁, k)
                p₁ = p₁ * w₁^k
            end
        end
    end

    q, r = divide_poly(p, p₁)

    if !iszero(r)        
        @warn "incorrect factorization"
        println(p₁)
    end

    if !iszero(q - 1)        
        add_factor!(h, q, 0)
        # push!(h, (q, 0))        
    end

    return h 
end
