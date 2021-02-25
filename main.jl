setprecision(100)
T = BigFloat
PI = T(π)
πi = PI * 1im
πi2 = 2πi

function γ(z, ϵ)
    z = z - 1
    half = T(1) / T(2)
    a = floor(max(real(z), 2)) + half
    while true
        ϵ̂ = √(a) / (2 * PI)^(a + half) / (real(z) + a)
        ϵ̂ < ϵ && break
        a += 1
    end

    coef = 1
    k = 1
    fac = T(1)
    while k < a
        c_k = (a - k)^(k - half) * exp(a - k) * (-1)^(k - 1) / √(2 * PI) / fac
        fac *= k
        coef += c_k / (z + k)
        k += 1
    end
    coef * (z + a)^(z + half) * √(2 * PI) / exp(z + a)
end

function brentq(f, xa, xb, xtol, rtol, iter)
    xpre, xcur, xblk, fpre, fcur, fblk, spre, scur = xa, xb, 0, f(xa), f(xb), 0, 0, 0

    fpre * fcur > 0 && return 0
    fpre == 0 && return xpre
    fcur == 0 && return xcur

    for _ in 1:iter
        if fpre * fcur < 0
            xblk = xpre
            fblk = fpre
            spre = xcur - xpre
            scur = spre
        end

        if abs(fblk) < abs(fcur)
            xpre = xcur
            xcur = xblk
            xblk = xpre

            fpre = fcur
            fcur = fblk
            fblk = fpre
        end

        Δ = (xtol + rtol * abs(xcur)) / 2
        sbis = (xblk - xcur) / 2
        (fcur == 0 || abs(sbis) < Δ) && return xcur

        if abs(spre) < Δ && abs(fcur) < abs(fpre)
            if xpre == xblk
                stry = -fcur * (xcur - xpre) / (fcur - fpre)
            else
                dpre = (fpre - fcur) / (xpre - xcur);
                dblk = (fblk - fcur) / (xblk - xcur);
                stry = -fcur * (fblk * dblk - fpre * dpre) / (dblk * dpre * (fblk - fpre))
            end

            if 2 * abs(stry) < min(abs(spre), 3 * abs(sbis) - Δ)
                spre = scur
                scur = stry
            else
                spre, scur = sbis, sbis
            end
        else
            spre, scur = sbis, sbis
        end

        xpre = xcur
        fpre = fcur

        if abs(scur) > Δ
            xcur += scur
        elseif sbis > 0
            xcur += Δ
        else
            xcur -= Δ
        end

        fcur = f(xcur)
    end

    xcur
end

function I0(s, ϵ)
#     println("[I0] s = $s")
    g(z, s) = z^(-s) * exp(z * z * πi)
    H(w) = 1 / (1 - exp(πi2 * w))

    function f(z, s)
        a = exp(z * πi)
        g(z, s) / (a - 1 / a)
    end

    σ = real(s)
    z_s = √(s / πi2)
    n = max(floor(real(z_s) - imag(z_s)), 1)
#     println("z_s = $z_s")

    Δ = max(0, log(n^(-σ) / ϵ) / (2 * PI))
#     println("σ = $σ, n = $n, n^-σ = $(n^-σ)")
    direction = (1 + 1im) / √ T(2)

    criterion(α) = abs(f(n + 0.5 + direction * α * Δ, s)) - ϵ
    α_1 = brentq(criterion, 0,  2, ϵ, ϵ, 20)
    α_2 = brentq(criterion, 0, -2, ϵ, ϵ, 20)

    z_1 = n + 1/2 + direction * α_1 * Δ
    z_2 = n + 1/2 + direction * α_2 * Δ
#     println("α_1 = $α_1, α_2 = $α_1, $(Δ), $σ, $(n^(-σ))")

    z_l = 0
    z_r = 0
    h = 0

    m = 1
    found = false
    for i in 1:1000
        h = (z_2 - z_1) / m
        inv_2h = 1 / (2 * h)
        z_l = sqrt(inv_2h^2 + s / πi2) + inv_2h
        z_r = sqrt(inv_2h^2 + s / πi2) - inv_2h
        err_l = g(z_l, s) * exp((z_1 - z_l) / h * πi2)
        err_r = g(z_r, s) * exp((z_r - z_1) / h * πi2)
        if abs(err_l) ≤ ϵ && abs(err_r) ≤ ϵ
            found = true
            break
        end
        m = ceil(m * 1.1)
#         println("err_l = $err_l, err_r = $err_r, g(z_l, s) = $(g(z_l, s)), z_l = $z_l, z_r = $z_r")
    end
    @assert found

    n_l = ceil(real(z_l) - imag(z_l))
    n_r = floor(real(z_r) - imag(z_r))

    fs2(i) = i^(-s) * H((i - z_1) / h)
    fs3(i) = i^(-s) * H((z_1 - i) / h)
#     println("[ζ] n = $n, m = $m, n_l = $n_l, n_r = $n_r")

    s0 = sum(@. (1:n)^(-s))
    s1 = sum(@. f((0:m) * h + z_1, s))
    s2 = sum(@. fs2(n_l:n))
    s3 = sum(@. fs3((n+1):n_r))

    s0 + s1 * h - s2 + s3
end

χ(s, ϵ) = (2 * PI)^s / cos(PI * s / 2) / γ(s, ϵ) / 2
# ζ(s, ϵ) = I0(s, ϵ / 2) + χ(s, ϵ) * conj(I0(1 - conj(s), ϵ / 2))
function ζ(s, ϵ)
    ret = I0(s, ϵ / 2) + χ(s, ϵ) * conj(I0(1 - conj(s), ϵ / 2))
    println("[ζ] s = $s, ϵ = $ϵ, ret = $ret, $(χ(s, ϵ))")
    ret
end


function erfc(z, ϵ)
    if abs(z) > 2
        h = π / √(log(6 / ϵ))
        K = ceil(√(log(1 / ϵ)) / h)
        ret = 1 / z^2
        PI = T(π)
        for k in 1:K
            ret += 2 * exp(-h^2 * k^2) / (z^2 + h^2 * k^2)
        end
        return ret * exp(-z^2) * h * z / PI + 2 / (1 - exp(2 * PI * z / h))
    else
        if z ≥ 0
            s = 1
        else
            s = -1
            z = -z
        end

        ϵ0 = ϵ / 2
        K = ceil(log2(1 / ϵ) / 2)
        ϵ1 = ϵ0 / 2K
        ϵ1 = 2^floor(log2(ϵ1))

        t = z
        k = T(0)
        S = T(0)
        while abs(t) / (2k + 1) > ϵ0
            S = S + t / (2k + 1)
            t = -z^2 * t / (k + 1) # error in paper
            k += 1
        end
        S += t / (2k + 1)
        1 - s * S * 2 / √π
    end
end


using SpecialFunctions


function find_T(x, λ, σ, ϵ)
#     println("[find T] x = $x, λ = $λ, σ = $σ, ϵ = $ϵ")
    # only to determine bounds, so it's fine to use Float64 here.

    # σ is real and 1 ≤ σ ≤ 2

    E_Σ(T, x, λ, σ) = exp(λ^2 * σ^2 / 2) / 2π * x^σ * log(zeta(σ)) * expint(λ^2 * T^2 / 2)

    # now we solve min T s.t. E_Σ(T, x, λ, σ) < 0.75ϵ
    # T = √(2u) / λ where E_1(u) < 0.75 ϵ / exp(λ^2 * σ^2 / 2) / (2 * π) * log(gamma(σ))
    limit = 0.75ϵ / (exp(λ^2 * σ^2 / 2) / (2 * π) * log(zeta(σ))) / x^σ
#     println("[find T] limit = $limit, expint(0) = $(expint(0))")
    f(u) = expint(u) - limit

    u = brentq(f, λ^2 / 2, -log(limit), ϵ, ϵ, 100)
#     println("[debug] est u = $(-log(limit)), u = $u")
    T = √(2u) / λ
    println("[param] T = $T, residual = $(E_Σ(T, x, λ, σ)) < ϵ = $ϵ")

    T
end

function quad_π_star(x, σ, λ, ϵ)
    ϵ = ϵ / 4
    T = find_T(x, λ, σ, ϵ)

    h_L = 2 * π * (σ - 1) / (log(x / ϵ) + λ^2 / 2 + 1 / x)
    h_R = 2 * π / (log(x / 2) + σ * λ^2 + λ * √(2 * σ * log(x / 2) + σ^2 * λ^2 + 2 * log(3.4 / ϵ)))
    h = min(h_L, h_R)
    println("[param] h = $h, h_L = $h_L, h_R = $h_R")

    ϵ_0 = ϵ / x^σ / log(x)

    Ψ(s, x, λ, ϵ) = exp((λ * s)^2 / 2) * x^s * log(ζ(s, ϵ)) / s

    results = @. real(Ψ(σ + 1im * h * (1:ceil(T/h)+10), x, λ, ϵ_0))
#     println(results)
    s = sum(results)
    h / π * (Ψ(σ, x, λ, ϵ_0) / 2 + s)
end

function get_primes(l, r)
    a = ones(Bool, r - l + 1)
    for i in 2:Int64(floor(√r))
        li = max((l - 1) ÷ i + 1, 2)
        ri = r ÷ i
        for j in li : ri
            a[i * j - l + 1] = 0
        end
    end

    ret = []
    for i in max(l, 2):r
        if a[i - l + 1]
            push!(ret, i)
        end
    end
    ret
end


function Δ(x, λ, ϵ)
    Φ(ρ, ϵ) = erfc(ρ / √2, ϵ) / 2
    ϕ(u, x, λ, ϵ) = Φ(log(u / x) / λ, ϵ)
    E_plus(u, x, λ, ϵ) = x * exp(λ^2 / 2) * Φ(log(u / x) / λ - λ, ϵ / x) - u * Φ(log(u / x) / λ, ϵ / 8 / u)
    E_minus(u, x, λ, ϵ) = u * Φ(-log(u / x) / λ, ϵ / 8 / u) - x * exp(λ^2 / 2) * Φ(λ - log(u / x) / λ, ϵ / x)

    ϵ = ϵ / 2
    x1 = floor(x)
    if E_minus(x1, x, λ, ϵ) > ϵ
        g(u) = E_minus(u, x, λ, ϵ) - 0.75 * ϵ
        u = brentq(g, 2, x, 1e-8, 1e-8, 100)
        x1 = floor(u)
    end

    x2 = ceil(x)
    if E_plus(x2, x, λ, ϵ) > ϵ
        f(u) = E_plus(u, x, λ, ϵ) - 0.75 * ϵ
        u = brentq(f, x * 2, x, 1e-8, 1e-8, 100)
        x2 = ceil(u)
    end
    println("[Δ] x1 = $x1, x2 = $x2")
    println("$(E_minus(900, x, λ, ϵ))")
    println("$(E_plus(x2, x, λ, ϵ))")

    N = (x2 - x1 + 1) + ceil(√x2)
    ϵ_1 = ϵ * 2^ceil(-log2(N))

    ret = 0
    for p in get_primes(Int64(x1), Int64(x2))
        ret -= ϕ(T(p), x, λ, ϵ_1)
        if p ≤ x
            ret += 1
        end
    end
    for p in get_primes(2, Int64(floor(√x2)))
        m = 2
        while p^m < x2
            if p^m < x1
                ret -= T(1 / m)
            else
                ret -= ϕ(p^m, x, λ, ϵ_1) / m
            end
            m += 1
        end
    end

    ret
end


function analytic_π(x)
    μ = 1 + 0
    σ_LR = 2
    σ_0 = (μ + 1/2) / (μ - 1/2)
#     σ = min(σ_LR, σ_0)
    σ = 1.5
    println("[param] σ = $σ")

    λ = 1 / (√x * log(x)^0.25)
    println("[param] λ = $λ")

    ϵ = 1/4
    pi_star = quad_π_star(x, σ, λ, ϵ / 2)
    Delta = Δ(x, λ, ϵ / 2)
    println("pi_star = $pi_star, Delta = $Delta")
    round(Delta + pi_star)
end

print(analytic_π(10^5))
