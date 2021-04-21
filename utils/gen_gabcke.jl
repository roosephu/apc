
function rs_find_λ(T, n)
    # https://oeis.org/A087617
    λ = OffsetArray(zeros(T, n + 1), 0:n)
    λ[0] = 1
    eulers = rs_find_eulers(T, n)
    for i in 0:n-1
        c = 0
        for k in 0:i
            c += 2^(4k + 1) * eulers[k + 1] * λ[i - k]
        end
        λ[i + 1] = c / (i + 1)
    end

    λ
end

function rs_find_eulers(T, n)
    # https://oeis.org/A000364
    a = OffsetArray(zeros(T, n + 1), 0:n)
    a[0] = 1
    for i in 1:n
        for j in 0:i-1
            a[i] -= (-1)^(i + j) * a[j] * T(binomial(BigInt(2 * i), BigInt(2 * j)))
        end
    end
    a
end

function rs_find_d(T, m)
    d = OffsetArray(zeros(T, m + 1, m + 1), 0:m, 0:m)
    λ = rs_find_λ(T, m)
    for i in 0:m÷4
        d[i * 4, i * 3] = λ[i]
    end
    for n in 0:m-1
        for k in 0:(3n+2)÷4  # k < 3(n+1)/4  <=> k <= (3n+2)/4
            d[n + 1, k] = (3n + 1 - 4k) * (3n + 2 - 4k) * d[n, k] + (k >= 1 ? d[n, k - 1] : 0)
        end
    end

    d
end

function get_last_nonzero(ps)
    n = length(ps) - 1
    for i in n:-1:0
        if ps[i] != 0
            return i
        end
    end
    0
end

function ps_mul(a, b)
    n = length(a) - 1
    c = OffsetArray(zeros(typeof(a[0]), n+1), 0:n)
    max_i = get_last_nonzero(a)
    max_j = get_last_nonzero(b)
    for i in 0:min(n, max_i)
        for j in 0:min(n-i, max_j)
            c[i + j] += a[i] * b[j]
        end
    end
    c
end

function ps_div(a, b)
    n = length(a) - 1
    c = OffsetArray(zeros(typeof(a[0]), n+1), 0:n)
    for i in 0:n
        v = a[i]
        for j in 0:i-1
            v = v - b[i - j] * c[j]
        end
        c[i] = v / b[0]
    end
    c
end

function ps_compose(coeffs, ps)
    x = ps[0]
    n = length(ps) - 1
    a = copy(ps)
    a[0] = 0

    b = OffsetArray(zeros(typeof(a[0]), n+1), 0:n)
    for i in n:-1:0
        b = ps_mul(b, a)
        b[0] += coeffs[i]
    end
    b
end

function ps_apply(derivatives, ps)
    x = ps[0]
    n = length(ps) - 1
    a = copy(ps)
    a[0] = 0

    b = OffsetArray(zeros(typeof(a[0]), n+1), 0:n)
    factorial = typeof(a[0])(1.0)
    for i in 0:n
        derivatives[i] /= factorial
        factorial *= i + 1
    end
    for i in n:-1:0
        b = ps_mul(b, a)
        b[0] += derivatives[i]
    end
    b
end

function ps_cos(ps)
    n = length(ps) - 1
    x = ps[0]

    cos_x = cos(x)
    sin_x = sin(x)

    derivatives = [-sin_x, -cos_x, sin_x, cos_x]
    fn = OffsetArray(zeros(typeof(x), n+1), 0:n)
    for i in 0:n
        fn[i] = derivatives[mod1(i, 4)]
    end
    fn

    ps_apply(fn, ps)
end

function rs_find_ps(z, n)
    @assert n >= 2

    T = typeof(z)
    numer = OffsetArray(zeros(T, n + 1), 0:n)
    n_pi = T(π)
    numer[0] = z^2 * n_pi / 2 + 3 * n_pi / 8
    numer[1] = n_pi * z
    numer[2] = n_pi / 2
    numer = ps_cos(numer)

    denom = OffsetArray(zeros(T, n + 1), 0:n)
    denom[0] = n_pi * z
    denom[1] = n_pi
    denom = ps_cos(denom)

    ps = ps_div(numer, denom)
#     println(ps_mul(denom, ps) - numer)
    ps
end

function diff_ps(ps)
    n = length(ps) - 1
    ps = copy(ps)
    for i in 1:n
        ps[i - 1] = ps[i] * i
    end
    ps[n] = 0
    ps
end

setprecision(2000)

n = 10
m = 1000

d = rs_find_d(T, n)
ps = rs_find_ps(T(0), m)
derivatives = OffsetArray(zeros(T, 3n+1, m + 1), 0:3n, 0:m)
derivatives[0, :] = ps
for i in 1:3n
    derivatives[i, :] = diff_ps(derivatives[i-1, :])
end

C = OffsetArray(zeros(T, n+1, m + 1), 0:n, 0:m)
for j in 0:n
    for k in 0:(3j ÷ 4)
        C[j, :] = C[j, :] + derivatives[3j - 4k, :] / factorial(T(3j - 4k)) * d[j, k] / T(π)^(2j - 2k)
    end
    C[j, :] /= T(2)^(2j)
end


fac = OffsetArray(zeros(T, 1001), 0:1000)
for i in 0:1000
    fac[i] = factorial(T(i))
end

function binom(n::Int, m::Int)
    fac[n] / fac[m] / fac[n - m]
end



function transform_to_chebyshev(ps) # transform a power series to Chebyshev's
    γ = OffsetArray(zeros(T, m+1), 0:m)

    # even case
    if ps[0] == 0
        ps = OffsetArray([ps[1:m]; 0], 0:m)
    end
    for k in 0:m÷2
        for l in k:m÷2
            γ[2k] += binom(2l, l-k) * ps[2l] / T(2)^(2l - 1)
        end
    end
    γ
end


γ = OffsetArray(zeros(T, n+1, m+1), 0:n, 0:m)
for j in 0:n
    γ[j, :] = transform_to_chebyshev(C[j, :])
end

output_m = 50
println("[")
for i in 0:n
    println("  [")
    for j in 0:2:output_m
        println("    $(repr(j == 0 ? γ[i, 0] / 2 : γ[i, j])),")
    end
    println("  ], ")
end
println("]")
