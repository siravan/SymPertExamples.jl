consensus(u, v) = (isequal(u.x,v.x) || v.x == nothing ? u.x : (u.x == nothing ? v.x : nothing))

function Base.:+(u::Poly, v::Poly)
    x = consensus(u, v)
    x === nothing && return poly(sym(u) + sym(v), nothing)
    q = Poly(u.T, x)

    for k = 0:max(u.n, v.n)
        q[k] = u[k] + v[k]
    end
    q
end

Base.:+(u, v::Poly) = polyas(v, u) + v
Base.:+(u::Poly, v) = u + polyas(u, v)

function Base.:-(u::Poly, v::Poly)
    x = consensus(u, v)
    x === nothing && return poly(sym(u) - sym(v), nothing)
    q = Poly(u.T, x)

    for k = 0:max(u.n, v.n)
        q[k] = u[k] - v[k]
    end
    q
end

Base.:-(u::Poly, v) = u - polyas(u, v)
Base.:-(u, v::Poly) = polyas(v, u) - v


function Base.:*(u::Poly, v::Poly)
    x = consensus(u, v)
    x === nothing && return poly(sym(u) * sym(v), nothing)
    q = Poly(u.T, x)

    for i = 0:u.n
        for j = 0:v.n
            q[i+j] += u[i] * v[j]
        end
    end
    q
end

Base.:*(u::Poly, v) = u * polyas(u, v)
Base.:*(u, v::Poly) = polyas(v, u) * v

function divide_poly(u::Poly, v::Poly; op=/)
    x = consensus(u, v)
    x === nothing && return poly(sym(u) / sym(v), nothing)
    q = Poly(u.T, x)
    r = copy(u)
    m, n = u.n, v.n

    for i = m-n:-1:0
        q[i] = op(r[i+n], v[n])
        for j = n+i:-1:i
            r[j] = r[j] - q[i] * v[j-i]
        end
    end

    update_order!(r)

    q, r
end

function Base.:/(u::Poly, v::Poly)
    q, r = divide_poly(u, v)
    q
end

Base.:/(u::Poly, v) = u / polyas(u, v)
Base.:/(u, v::Poly) = polyas(v, u) / v

function Base.:%(u::Poly, v::Poly)
    q, r = divide_poly(u, v)
    r
end

Base.:%(u::Poly, v) = u % polyas(u, v)
Base.:%(u, v::Poly) = polyas(v, u) % v

function pseudo_divide_poly(u::Poly, v::Poly)
    m, n = u.n, v.n
    α = v[n]^(m-n+1)
    q, r = divide_poly(u*α, v; op=÷)
    return q, r, α
end

Base.:^(p::Poly, k) = poly(sym(p)^k, p.x)

##############################################################################

function gcd_naive(u::Poly, v::Poly)
    while degree(v) > 0   
        println(u, '\t', v)     
        u, v = v, u % v
    end    
    return iszero(v) ? u : 1
end

function gcd_monic(u::Poly, v::Poly)
    if degree(u) == 0 || degree(v) == 0
        return gcd(cont(u), cont(v))
    end
    u = u / leading(u)
    while degree(v) > 0        
        v = v / leading(v)
        u, v = v, u % v
    end
    return iszero(v) ? u : 1
end

function gcd_extended(u::Poly, v::Poly; op=/)
    if degree(u) == 0 || degree(v) == 0
        return gcdx(cont(u), cont(v))
    end    
    
    sᵤ = one(u.T)
    tᵤ = zero(u.T)
    sᵥ = zero(v.T)
    tᵥ = one(v.T)
    
    while degree(v) > 0
        q, r = divide_poly(u, v; op)
        s, t = sᵤ - q*sᵥ, tᵤ - q*tᵥ
        u, sᵤ, tᵤ = v, sᵥ, tᵥ
        v, sᵥ, tᵥ = r, s, t          
    end

    l = leading(u)
    if iszero(l)
        return u, sᵤ, tᵤ
    else
        return u / l, sᵤ / l, tᵤ / l
    end
end

function Base.gcd(u::Poly, v::Poly)
    if u.T <: Rational && u.T <: Rational
        return gcd_monic(u, v)
    else
        return gcd_naive(u, v)
    end
end

Base.gcd(u::Poly, v) = gcd(u, polyas(u, v))
Base.gcd(u, v::Poly) = gcd(polyas(v, u), v)

Base.gcdx(u::Poly, v::Poly) = gcd_extended(rationalize(u), rationalize(v))
Base.gcdx(u::Poly, v) = gcdx(u, polyas(u, v))
Base.gcdx(u, v::Poly) = gcdx(polyas(v, u), v)

#############################################################################

Symbolics.derivative(p::Poly, x) = isequal(p.x, x) ? derivative(p) : derivative(sym(p), x)

function Symbolics.derivative(p::Poly)
    q = Poly(p.T, p.x)
    for i = 1:p.n
        q[i-1] = i * p[i]
    end
    q
end

function integrate(p::Poly)
    q = Poly(p.T, p.x)
    for i = 0:p.n
        q[i+1] = p[i] / (i+1)
    end
    q
end

###############################################################################

function rationalize(p::Poly)
    p.T <: Rational && return p
    try
        q = Poly(Rational{BigInt}, p.x, p.n, Rational{BigInt}(p.a₀), Dict())
        for k = 1:p.n
            q[k] = Rational{BigInt}(p[k])
        end
        return q
    catch
        return nothing
    end
end

function to_monic(p::Poly)
    p = rationalize(p)
    c = leading(p)
    n = degree(p)
    q = Poly(p.T, p.x)
    for i = 0:n
        q[i] = p[i] * c ^ (n-1-i)
    end
    c, q
end

function from_monic(p::Poly, c)
    q = Poly(p.T, p.x)
    n = degree(p)
    for i = 0:n
        q[i] = p[i] / c ^ (n-1-i)
    end
    q
end

function integer_poly(p::Poly)
    p = rationalize(p)
    l = lcm(denominator.(extract_coef(p))...)
    return 1//l, p*l    
end 
