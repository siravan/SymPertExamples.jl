@syms x

function generate_rand_poly(x; min_deg=0, max_deg=20, sparcity=0.5, rational=0.25)
    n = rand(min_deg:max_deg)
    p = rand() < rational ? Poly(Rational, x) : Poly(x)
    r = rand(1:10) // rand(1:10)
    for i = 0:n
        if rand() > sparcity || i == n
            s = rand() < 0.5 ? +1 : -1
            p[i] = (p.T == Rational ? s * r * rand(0:10)//rand(1:5) : s * 10 * randn()) 
        end
    end
    p
end

function test_eq(x, f; n=10, abstol=1e-8, min_deg=0, max_deg=20, sparcity=0.5, rational=0.25)
    k = 0
    for i = 1:n
        p = generate_rand_poly(x; max_deg, sparcity, rational)
        q = generate_rand_poly(x; min_deg, max_deg, sparcity, rational)      
        try  
            r = f(p, q)
            if abs(maximum(extract_coef(r))) > abstol            
                @warn "+/- mismatch: $r"  
                println(p, ", ", q)
            else
                k += 1
            end
        catch e
            println(e)
            println(p, ", ", q)
        end        
    end
    @info "$k ok!"
end

test_addsub(x) = test_eq(x, (p,q)->(p+q)-p-q)
test_muldiv(x) = test_eq(x, (p,q)->p-(p/q)*q-(p%q))
test_gcd(x) = test_eq(x, (p,q)->p%gcd_naive(p,q)+q%gcd_naive(p,q); max_deg=5, rational=1.0)

function test_deriv(x; n=10, min_deg=0, max_deg=20, sparcity=0.5, rational=1)
    k = 0
    for i = 1:n
        p = generate_rand_poly(x; min_deg, max_deg, sparcity, rational)
        try 
            q = derivative(p)
            Δp = p - integrate(q)
            if degree(Δp) > 0
                @warn "deriv/integrate mismatch: $Δp"
            else
                k += 1
            end
        catch e
            println(e)
        end
    end 
    @info "$k ok!"
end

function test_shubert(x; n=10, min_deg=1, max_deg=3, sparcity=0.5, rational=1)
    k = 0
    for i = 1:n
        p1 = generate_rand_poly(x; min_deg, max_deg, sparcity, rational)
        p2 = generate_rand_poly(x; min_deg, max_deg, sparcity, rational)
        p3 = generate_rand_poly(x; min_deg, max_deg, sparcity, rational)
        _, p1 = integer_poly(p1)
        _, p2 = integer_poly(p2)
        _, p3 = integer_poly(p3)
        p = p1 * p2 * p3
        try 
            println('(',p1,") * (",p2,") * (",p3,')')
            println("\t=> ", factor_schubert_kronecker(p))
            k += 1
        catch e
            println(e)
        end
    end 
    @info "$k ok!"
end

function test_factor(x)
    ps = [
        gen_poly(x,  64, 56, 14, 1),
        gen_poly(x,  -12, 53, -57, 18),
        gen_poly(x,  273, -86, -73, 6),
        gen_poly(x,  -2, -1, 4, -1, 6),
        gen_poly(x,  -30, -33, 8, -11, 6),
        gen_poly(x,  -2, 7, 20, -24, -6, 5),
        gen_poly(x,  8, -38, 27, 47, -11, 15),
        poly(x^4 - 4, x),
        poly(x^4 -8x^2 - 9, x),
        poly(6x^4 - 7x^3 + 5x^2 - 20x + 17, x),
        poly(x^6 - 1, x),
        poly(x^6 + 1, x),
        poly(x^5 + x + 1, x),
        gen_poly(x,  -35, 11, 6),
        gen_poly(x,  -16, 25),
        gen_poly(x,  4, -13, 24, -19, 6),
        gen_poly(x,  -2, 3, 5, 4, 2, 2, 2),
        gen_poly(x,  -27, 117, -90, 90, -48, 8),
        gen_poly(x,  2, 9, 25, 35, 39, 30),
        gen_poly(x,  -1, 1, 2, -2, -1, 1),
        gen_poly(x,  1, 2, 1, 1, 2, 1, 1, 2, 2, 2, 1),
        poly(x^8 - 4x^6 + 16x^2 - 16, x),
        gen_poly(x,  -1, 0, 3, 2, -3, -6, 6, 3, -2, -3, 1),
    ]

    k = 0
    for p = ps        
        try 
            printstyled(p, '\n'; color=:green)
            f = factor(p; verbose=true)
            println("\t", f)
            k += 1
            # q = expand(*(f...))
            # if iszero(p - q)
            #     @info "OK!"               
            #     k += 1
            # else
            #     @warn "mismatch!"
            # end
        catch e
            println(e)
        end
    end 
    @info "$k ok!"
end