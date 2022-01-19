use crate::cache_stat::CacheStat;

/// Compute $∑_{u} \phi(u)$ given a list of $u$.
///
///
/// To avoid precision loss of converting `u64` to `f64`, we calculate $\tilde
/// \phi(t) = \phi(x + t)$, with $t = O(\lambda x)$.
///
/// To calculate $\tilde \phi(t)$ fast, we simply use the order-2 Taylor
/// polynomial at $t_0$.  $$ g(t) = \tilde \phi(t_0) + (t - t_0) \tilde
/// \phi'(t_0) + \frac{1}{2}(t - t_0)^2 \tilde \phi''(t_0).  $$ where $$
/// \begin{aligned} \tilde \phi'(t) \& = -\frac{C}{x + t}\lambda^{-1}, \\ \tilde
/// \phi''(t) \& = \frac{C}{(x+t)^2} (\lambda^{-1} + \rho \lambda^{-2}), \\
/// \end{aligned} $$ for $\rho = \lambda^{-1} \ln (1+t/x), C =
/// \frac{\exp(-\rho^2 /2)}{\sqrt{2 \pi}}$. Note that we sometimes find a new
/// $t_0$ and build a new order-2 Taylor polynomial. So we need to determine
/// when the approximation is good. More specifically, we'd like to know what's
/// the maximum $w > 0$, such that, it hold that have $\left|g(t) - \tilde
/// \phi(t)\right| \leq \epsilon$ for all $t$ in $[t_0 - w, t_0 + w]$. By using
/// the Lagrange remainder, we have $$ |g(t) - \tilde \phi(t)| \leq \frac{1}{6}
/// w^3 \max_t |\tilde \phi'''(t)|, $$ where $$ \tilde \phi'''(t) =
/// -\frac{C}{(x+t)^3} (2 \lambda^{-1} + 3 \rho \lambda^{-2} - (1 - \rho^2)
/// \lambda^{-3}).  $$ Note that $\exp(-\rho^2 / 2) \rho \geq -0.63$ and
/// $\exp(-\rho^2 / 2) \rho^2 \geq 0$ for every $\rho$, we can then bound
/// $|\tilde \phi'''(t)|$ (for most $\lambda$): $$ |\tilde \phi'''(t)| \leq
/// \frac{-2 \lambda^{-1} + 1.8 \lambda^{-2} + \lambda^{-3}}{x^3 \sqrt{2 \pi}}.
/// $$
///
/// We limit $w$ to be smaller than $2^{21}$, so that for all primes $p$ in
/// $[t_0, t_0 + w)$, the sum of $(p - t_0)^2$ fits into a 64-bit integer: $$
/// \sum_{t_0 \leq p < t_0 + w} (p - t_0)^2 \leq \frac{1}{3} w^3 < 2^{63}, $$ so
/// that we don't need 128-bit integers. Making most operations in integer type
/// slightly increase the speed of calculating $\phi$ (10%), although the
/// bottleneck is sieving primes. It turns out that it's fine to use
/// high-precision data type to compute $\Delta$, as most operations are integer
/// operations. High-precision data type is only used to compute the
/// coefficients $\tilde \phi(t_0), \tilde \phi'(t_0), \tilde \phi''(t_0)$,
/// which are re-computed once per nearly 25k primes (when $x = 10^{18}$ and
/// $\lambda \approx 4.55 \times 10^{-8}$).
pub struct LittlePhiSum {
    radius: i64,
    eps: f64,
    t: i64,
    x: f64,
    λ: f64,
    coeffs: [f64; 3],
    sum: f64,
    sum0: i64,
    sum1: i64,
    sum2: i64,
    stat: CacheStat,
}

impl LittlePhiSum {
    pub fn stat(&self) -> &CacheStat { &self.stat }

    pub fn sum(&mut self) -> f64 {
        self.flush();
        self.sum
    }

    pub fn new(λ: f64, x: f64, eps: f64) -> Self {
        Self {
            radius: -1,
            eps,
            t: 0,
            x,
            λ,
            coeffs: [0.0, 0.0, 0.0],
            stat: CacheStat::new(),
            sum0: 0,
            sum1: 0,
            sum2: 0,
            sum: 0.0,
        }
    }

    pub fn once(&self, t: i64) -> f64 {
        rgsl::error::erfc((t as f64 / self.x).ln_1p() / self.λ / std::f64::consts::SQRT_2) / 2.0
    }

    #[inline]
    fn fast(&mut self, t: i64) -> Option<()> {
        let δ = t - self.t;
        if δ < -self.radius || δ > self.radius {
            None
        } else {
            self.stat.hit();
            self.sum0 += 1;
            self.sum1 += δ;
            self.sum2 += δ * δ;
            Some(())
        }
    }

    #[inline]
    fn flush(&mut self) {
        let s0 = self.sum0 as f64;
        let s1 = self.sum1 as f64;
        let s2 = self.sum2 as f64;
        self.sum += s0 * self.coeffs[0] + s1 * self.coeffs[1] + s2 * self.coeffs[2];

        self.sum0 = 0;
        self.sum1 = 0;
        self.sum2 = 0;
    }

    fn slow(&mut self, t: i64) {
        self.stat.miss();
        self.flush();

        let λ = self.λ;
        let x = self.x;
        let t_ = t as f64;

        let ρ = (t_ / x).ln_1p() / λ;
        let y0 = rgsl::error::erfc(ρ / std::f64::consts::SQRT_2) / 2.0;

        let coeff = (ρ * ρ / -2.0).exp() / (std::f64::consts::PI * 2.0).sqrt();
        let y1 = coeff / (x + t_).powi(1) * -(1.0 / λ);
        let y2 = coeff / (x + t_).powi(2) * (1.0 / λ + ρ / λ / λ);

        // See Readme.md for derivation.
        let max_y3 = (-2.0 / λ + 1.8 / λ / λ + 1.0 / λ / λ / λ)
            / (std::f64::consts::PI * 2.0).sqrt()
            / x.powi(3);
        let radius = (self.eps / max_y3 * 6.0).cbrt() as i64;
        let radius = std::cmp::min(radius, 1 << 21);

        self.radius = radius;
        self.t = t;
        self.coeffs = [y0, y1, y2 / 2.0];

        self.sum0 = 1;
    }

    #[inline]
    pub fn add(&mut self, t: i64) { self.fast(t).unwrap_or_else(|| self.slow(t)) }
}
