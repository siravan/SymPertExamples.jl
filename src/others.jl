# @syms 𝑛 taylor(x,a,eq) b(n) D(x,eq)

is_taylor(x) = istree(x) && (operation(x)==taylor || any(is_taylor.(arguments(x))))
is_not_taylor(x) = !is_taylor(x)

taylor(x, 𝑛, a₀, eq) = a₀ + sum([x^i * substitute(eq, Dict(𝑛 => Float64(i))) for i=1:5])

tay_sum1 = @acrule taylor(~x,~a1,~eq1) + taylor(~x,~a2,~eq2) => taylor(~x, ~a1 + ~a2,~eq1 + ~eq2)
tay_sum2 = @acrule ~z::is_not_taylor + taylor(~x,~a,~eq) => taylor(~x,~z + ~a,~eq)
tay_sum3 = @acrule taylor(~x,~a1,~eq1) + ~u2::is_not_taylor * taylor(~x,~a2,~eq2) => taylor(~x,~a1 + ~a2*~u2,~eq1 + ~eq2*~u2)
tay_sum4 = @acrule ~u1::is_not_taylor * taylor(~x,~a1,~eq1) + ~u2::is_not_taylor * taylor(~x,~a2,~eq2) => taylor(~x,~a1*~u1 + ~a2*~u2,~eq1*~u1 + ~eq2*~u2)
tay_mul = @acrule ~z::is_not_taylor * taylor(~x,~a,~eq) => taylor(~x, ~a*~z,~eq*~z)
tay_div = @acrule taylor(~x,~a,~eq) / ~z::is_not_taylor => taylor(~x, ~a/~z,~eq/~z)

tay_diff1 = @rule D(~x, taylor(~x,~a,~eq)) => taylor(~x,D(~x,~a) + substitute(~eq, Dict(𝑛 => 1)), substitute(𝑛 * ~eq, Dict(𝑛 => 𝑛+1)))
tay_diff2 = @rule D(~x, b(~n)) => 0
tay_diff3 = @rule D(~x, ~y::is_not_taylor) => Differential(~x)(~y)

tay_rules = Chain([tay_mul, tay_div, tay_sum1, tay_sum2, tay_sum3, tay_sum4])
tay_diffs = Chain([tay_diff1, tay_diff2, tay_diff3])

function simplify_taylor(eq, order)
    w = z
    for i = 1:order
        w = Prewalk(PassThrough(tay_diffs))(w)
        w = expand_derivatives(w)
        w = Fixpoint(tay_rules)(w)
    end
    w
end

function realize_taylor(t, N)
    if istree(t) && operation(t)==taylor
        x, a, eq = arguments(t)

        eqs = [substitute(a, Dict(𝑛 => 0))]
        for i=1:N
            push!(eqs, substitute(eq, Dict(𝑛 => i)))
        end

        subs = Dict(b(i) => c[i+1] for i=0:10)
        eqs = [substitute(eq, subs) for eq in eqs]

        return eqs
    end
end

function apply_initial_values(eqs, u0)
    [substitute(eq,u0) for eq in eqs]
end

###############################################################################

@syms BAdd(x,y) BMul(x,y)

to_BAdd(x) = x
to_BAdd(x, y...) = BAdd(x, to_BAdd(y...))

to_BMul(x) = x
to_BMul(x, y...) = BMul(x, to_BMul(y...))

BAdd_rule = @rule +(~~xs) => to_BAdd(~~xs...)
BMul_rule = @rule *(~~xs) => to_BMul(~~xs...)

binary_op_rules = Chain([BAdd_rule, BMul_rule])

function convert_binary_op(x)
    x = Prewalk(PassThrough(BMul_rule))(x)
    Prewalk(PassThrough(BAdd_rule))(x)
end

iverson(x, y) = exp(-(x-y)^2/1e-8)

sum_rule1 = @rule BAdd(taylor(~x,~a1,~eq1), taylor(~x,~a2,~eq2)) => taylor(~x, ~a1 + ~a2, ~eq1 + ~eq2)
sum_rule2 = @rule BAdd(taylor(~x,~a,~eq), ~u::is_not_taylor) => taylor(~x, ~a + ~u, ~eq)
sum_rule3 = @rule BAdd(~u::is_not_taylor, taylor(~x,~a,~eq)) => taylor(~x, ~a + ~u, ~eq)
mul_rule1 = @rule BMul(taylor(~x,~a,~eq), ~u::is_not_taylor) => taylor(~x, ~a * ~u, ~eq * ~u)
mul_rule2 = @rule BMul(~u::is_not_taylor, taylor(~x,~a,~eq)) => taylor(~x, ~a * ~u, ~eq * ~u)
div_rule = @rule taylor(~x,~a,~eq) / ~u::is_not_taylor => taylor(~x, ~a / ~u, ~eq / ~u)
neg_rule = @rule -taylor(~x,~a,~eq) => taylor(~x, -~a, ~eq)

tay_diff1 = @rule D(~x, taylor(~x,~a,~eq)) => taylor(~x,D(~x,~a) + substitute(~eq, Dict(𝑛 => 1)), substitute(𝑛 * ~eq, Dict(𝑛 => 𝑛+1)))
tay_diff2 = @rule D(~x, b(~n)) => 0
tay_diff3 = @rule D(~x, ~y::is_not_taylor) => Differential(~x)(~y)

tay_rules = Chain([sum_rule1, sum_rule2, sum_rule3, mul_rule1, mul_rule2, div_rule, neg_rule])
tay_diffs = Chain([tay_diff1, tay_diff2, tay_diff3])

function simplify_taylor(eq, order)
    w = z
    for i = 1:order
        w = Prewalk(PassThrough(tay_diffs))(w)
        w = expand_derivatives(w)
    end
    for i = 1:order
        w = Prewalk(PassThrough(tay_rules))(w)
    end
    w
end

############################################################################

diff_vars(n) = diff_vars(n, 0)

function diff_vars(n, x₀)
    y = def_taylor(x - x₀, a[2:n], a[1])
    D = Differential(x)
    x, y, D, n
end

function eqs_linear_diffeq(p, eq)
    eq = expand_derivatives(eq)
    collect_powers(eq, p.x, 0:p.n)
end

function solve_linear_diffeq(p, eq, order)
    eqs = eqs_linear_diffeq(p, eq)
    vars = solve_coef(eqs, a[order+1:p.n])
    return vars
end

#############################################################################

@syms G(x)

g_rule1 = @rule G(~x) + G(~y) => G(~x + ~y)
g_rule2 = @acrule G(~x) + ~y => G(~x + ~y)
g_rule3 = @rule G(~x) * G(~y) => G(~x * ~y)
g_rule4 = @acrule G(~x) * ~y => G(~x * ~y)
g_rule5 = @rule G(~x) / G(~y) => G(~x // ~y)
g_rule6 = @rule G(~x) / ~y => G(~x // ~y)
g_rule7 = @rule ~x / G(~y) => G(~x // ~y)
g_rule8 = @rule ^(G(~x), ~y) => G(^(~x, ~y))
g_rule9 = @rule ^(G(~x), G(~y)) => G(^(~x, ~y))

g_rules = Chain([
    g_rule1,
    g_rule2,
    g_rule3,
    g_rule4,
    g_rule5,
    g_rule6,
    g_rule7,
    g_rule8,
    g_rule9,
])


#############################################################################
