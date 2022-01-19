use std::{
    path::{Path, PathBuf},
    time::Instant,
};

use crate::brentq::brentq;
use crate::platt_integral::*;
use crate::traits::*;
use log::{debug, info};
use num::integer::*;

#[doc(hidden)]
pub struct PlattBuilder {
    hint_λ: f64,
    poly_order: usize,
    segment: u64,
}

/// The algorihtm described in David Platt, “Computing $\pi(x)$ Analytically.”
///
/// Also see [ExpansionIntegrator] to compute $\hat\Phi(s)$ and
/// [HybridPrecIntegrator] for speedup.
#[derive(Default)]
pub struct Platt {
    x: u64,
    λ: f64,
    max_height: f64,
    segment: u64,
    x1: u64,
    x2: u64,
    max_order: usize,
    hint_λ: f64,
}

impl Default for PlattBuilder {
    fn default() -> Self { Self { hint_λ: 10.0, poly_order: 15, segment: 1u64 << 32 } }
}

impl PlattBuilder {
    /// During planning, these hyperparameters (λ, sigma, h, x1, x2,
    /// integral_limits) doesn't need to be very accurate as long as they
    /// satisfy the error bound.
    pub fn build(&self, x: u64) -> Platt {
        let x_ = x as f64;
        let λ = (self.hint_λ * x_.ln() / x_).sqrt();
        // Lemma 4.5
        let ignored = (λ * λ / 2.0).exp() * (12.533141373155 + 2.0 / λ)
            / (2.0 * std::f64::consts::PI * x_ * λ);
        info!("λ = {:.6e}, ignored term = {:.6}", λ, ignored);
        assert!(ignored < 0.1, "Too large ignored term. See [Lemma 4.5, Platt].");

        let (x1, x2) = Platt::plan_Δ_bounds_heuristic(λ, x_, 0.24);
        let max_height = plan_ζ_zeros(λ, x_, 0.1);

        Platt {
            x,
            λ,
            segment: self.segment,
            max_height,
            x1,
            x2,
            max_order: self.poly_order,
            hint_λ: self.hint_λ,
        }
    }

    pub fn hint_λ(&mut self, hint_λ: f64) -> &mut Self {
        self.hint_λ = hint_λ;
        self
    }

    pub fn poly_order(&mut self, poly_order: usize) -> &mut Self {
        self.poly_order = poly_order;
        self
    }

    pub fn sieve_segment(&mut self, segment: u64) -> &mut Self {
        self.segment = segment;
        self
    }
}

fn Φ(r: f64) -> f64 { rgsl::error::erfc(r / std::f64::consts::SQRT_2) / 2.0 }

// change of variable for err_l and err_r: t = ln(u/x) / λ <=> u = x e^(t λ).
#[inline]
fn err_r(t: f64, x: f64, λ: f64) -> f64 {
    x * ((λ * λ / 2.0).exp() * Φ(t - λ) - (t * λ).exp() * Φ(t))
}

#[inline]
fn err_l(t: f64, x: f64, λ: f64) -> f64 {
    x * ((t * λ).exp() * Φ(-t) - (λ * λ / 2.0).exp() * Φ(λ - t))
}

#[inline(never)]
fn calc_Δ1_f64(x: u64, eps: f64, λ: f64, x1: u64, x2: u64, segment: u64) -> f64 {
    let approx_n_primes = (x2 - x1) as f64 / (x1 as f64).ln();
    assert!(f64::EPSILON * approx_n_primes < 0.1, "too many primes for f64 approx");

    let mut ϕ = crate::fast_phi::LittlePhiSum::new(λ, x as f64, eps / (x2 - x1) as f64);
    let mut n_primes_less_than_x = 0u64;

    let mut calc = |p: u64| {
        let t = p as i64 - x as i64;
        if t <= 0 {
            n_primes_less_than_x += 1;
        }
        ϕ.add(t)
    };

    let mut n_primes = 0usize;

    let mut x1 = x1;
    let n_intervals = (x2 - x1) / segment + 1;
    for k in 1..=n_intervals {
        let y1 = std::cmp::min(x2, x1 + segment - 1);
        info!("sieving {}/{}: [{}, {}]", k, n_intervals, x1, y1);
        let sieve_result = crate::sieve::sieve_primesieve(x1, y1);
        for &p in sieve_result.iter() {
            n_primes += 1;
            calc(p);
        }
        x1 = y1 + 1;
    }
    info!("found {} primes in the interval", n_primes);

    ϕ.stat().show("Fast ϕ");
    let Δ_1 = n_primes_less_than_x as f64 - ϕ.sum();

    Δ_1
}

/// See Readme.md for why this can be computed using `f64`
#[inline(never)]
fn calc_Δ_f64(x: u64, eps: f64, λ: f64, x1: u64, x2: u64, segment: u64) -> f64 {
    let primes = crate::sieve::sieve_primesieve(1, x2.sqrt());

    let Δ_1 = calc_Δ1_f64(x, eps, λ, x1, x2, segment);
    let mut Δ_2 = 0.0;

    for &p in primes.iter() {
        let mut m = 1i64;
        let mut power = p;
        while power < x2 / p {
            m += 1;
            power *= p;
            if power < x1 {
                Δ_2 += 1.0 / m as f64;
            } else {
                // only (x_2^1/2 - x_1^1/2) + (x_2^1/3 - x_1^1/3) + ... many
                // The first term dominates, which is still O((x2 - x1)/sqrt(x)) = O(polylog(n)).
                let r = (power as f64 / x as f64).ln() / λ;
                Δ_2 += Φ(r) / m as f64;
            }
        }
    }
    Δ_1 - Δ_2
}

fn cramer_stats(x: f64, λ: f64, d: f64) -> (f64, f64, f64) {
    let integral_l = err_l(-d, x, λ);
    let integral_r = err_r(d, x, λ);
    let p = 1.0 / x.ln();
    let mean = (integral_r - integral_l) * p;
    let max = Φ(d);

    let var = p * (1.0 - p) * (integral_l + integral_r) * max;

    (mean, var, max)
}

// this requires high precision: result ≈ # primes
fn integrate_offline<T: MyReal>(x: u64, λ: f64, max_order: usize) -> T {
    // the result = n / log(n) + o(n)
    let mut integrator =
        PlattIntegrator::<T>::new(T::from_u64(x).unwrap(), T::one(), T::mp(λ), max_order, 0.01);
    let result = integrator.query(T::zero()).im;
    info!("Φ̂(1) = {}", result);

    result
}

fn plan_ζ_zeros(λ: f64, x: f64, eps: f64) -> f64 {
    // let T = (x.ln() + x.ln().ln()).sqrt() / λ;

    let ln_err = |T: f64| -> f64 {
        let PI = std::f64::consts::PI;
        let Tpow = (T / 2.0 / PI) * ((T / 2.0 / PI).ln() - 1.0)
            + 0.875
            + 1.588
            + 0.137 * T.ln()
            + 0.443 * T.ln().ln();
        let a = 2.0
            * (x.sqrt() / T / x.ln() + 1.0 / (λ * λ * T * T * x))
            * (2.0 / (λ * λ * T * T) * Tpow);
        a.ln() + λ * λ * (1.0 - T * T) / 2.0
    };

    let maxT = (x.ln() + x.ln().ln()).sqrt() / λ * 2.0;
    let ln_eps = eps.ln();
    let T = brentq(|t| ln_err(t) - ln_eps, 3.0, maxT, 0.0, 0.0, 30).unwrap();
    let err = ln_err(T).exp();

    info!("[plan ζ zeros] T = {:.6e}, trunc error = {:.6}, maxT = {:.6e}", T, err, maxT);

    T
}

impl Platt {
    pub fn compute<T: MyReal>(&mut self, db: impl ZetaZerosDatabase<T>) -> u64 {
        let max_order = self.max_order;
        let x = self.x;
        let λ = self.λ;

        let t0 = Instant::now();

        let integral_offline = integrate_offline::<T>(x, λ, max_order);
        let integral_critical = Self::integrate_critical::<T>(x, λ, max_order, self.max_height, db);

        let t1 = Instant::now();
        let Δ = calc_Δ_f64(x, 0.5, λ, self.x1, self.x2, self.segment);
        let t2 = Instant::now();

        info!("Δ = {}", Δ);
        let ans = integral_offline - integral_critical * 2.0 - 2.0.ln() + Δ;
        let n_ans = ans.round().to_u64().unwrap();
        info!("ans = {}, rounding error = {}", ans, ans - T::from_u64(n_ans).unwrap());

        let time_integral = (t1 - t0).as_secs_f64();
        let time_sieve = (t2 - t1).as_secs_f64();

        info!(
            "Time: ζ zeros = {:.3} sec, sieve = {:.3} sec, ratio = {:.3}",
            time_integral,
            time_sieve,
            time_integral / time_sieve
        );
        info!("Suggested: --lambda-hint {:.3}", self.hint_λ * time_integral / time_sieve);

        n_ans as u64
    }

    /// We might shorten the interval by applying some heuristics based on
    /// Cramér's random model.
    ///
    /// For brevity, we define $f(u) = \phi(u) - [u \leq x]$, i.e., $f(u) =
    /// \phi(u)$ for $u > x$ and $f(u) = \phi(u) - 1$ for $u \leq x$.
    ///
    /// Given a real $d > 0$,  define $U_L = [1, xe^{-d\lambda}], U_R = [xe^{d
    /// \lambda}, \infty)$ and $U = U_L \cup U_R$. We'd like to sieve all primes
    /// in $(xe^{-d\lambda}, xe^{d \lambda})$ to compute $\Delta$ and estimate
    /// the truncation error $S = \sum_{p \in U}f(p)$. Note that $\lambda =
    /// \tilde O(x^{-1/2})$ and $d = \tilde O(1)$, so $x e^{d\lambda} = x -
    /// \tilde O(x^{1/2})$ and $x e^{d\lambda} = x + \tilde O(x^{1/2})$.
    ///
    /// By Cramér's random model, we assume $$ S = \sum_{p \in U} f(p) \approx
    /// \sum_{u \in U} \frac{f(u)}{\ln u}, $$ One way to improve the
    /// approximation is to add an error bar, which is done by studying the
    /// behavior of a stochastic estimation $\hat S$: $$ \hat S = \sum_{u \in U}
    /// X_u f(u), \quad X_u \sim \text{Bernoulli}(1/\ln u).  $$ We use
    /// Bernstein's inequality to approximate $S = \mathbb{E} [\hat S] + t$ by
    /// solving $t$ for $\epsilon = \exp(-20)$ in $$ \exp\left(
    /// -\frac{\frac{1}{2} t^2}{\text{Var}[\hat S] + \frac{1}{3} M t} \right)
    /// \leq \epsilon, \quad M := \max_{u \in U} |f(u)| = \Phi(d).  $$ Let $q =
    /// 1/\ln x$, so $q \approx 1/\ln u$ for reasonable $u \in U$, i.e., those
    /// with not-too-small $f(u)$.
    ///
    /// Next, we approximate $\mathbb{E}[\hat S]$: $$ \mathbb{E}[\hat S] =
    /// \sum_{u \in U} \frac{f(u)}{\ln u} \approx \sum_{u \in U} q f(u) \approx
    /// q (I_r - I_l), $$ where $$ I_r = x \lambda \int_{d}^{\infty} \Phi(t)
    /// e^{t\lambda} \d t, \quad I_l = x \lambda \int_{d}^\infty \Phi(t)
    /// e^{-t\lambda} \d t.  $$ The last approximation is similar to [Section
    /// 3.2, Galway] and both $I_l$ and $I_r$ have a closed form solution.
    /// Unlike [Galway], we are making use of the fact that $f(u) < 0$ for $u <
    /// x$ and $f(u) > 0$ for $u > x$ and the sum of them might cancel.
    ///
    /// The variance of $\hat S$ can be approximated by: $$ \text{Var}[\hat S] =
    /// \sum_{u \in U} f(u)^2 \frac{1}{\ln u} (1 - 1/\ln u) \approx q(1-q)
    /// \sum_{u \in U} f(u)^2 \leq q(1-q) M \sum_{u \in U} |f(u)| \approx M
    /// q(1-q) (I_l + I_r).  $$ Finally, we minimize $d$ so that $S \leq 0.24$.
    /// In this way, we can get a shorter interval to sieve at a loss of
    /// guarantee.
    ///
    /// In practice, when $x = 10^{11}$ and $\lambda = 3 \times 10^{-5}$,
    /// [Galway] predicts an interval of length 3e7 while this predicts an
    /// interval of length 2e7, while the difference between two sums is around
    /// 0.0112.
    ///
    /// [Büthe] provides a rigorous way to reduce the interval, but I don't
    /// understand... [Theorem 4.3, Büthe2] seems pretty interesting and might
    /// be related to the heuristics I used.
    pub fn plan_Δ_bounds_heuristic(λ: f64, x: f64, eps: f64) -> (u64, u64) {
        let τ = 20.0;
        let err = |d: f64| {
            let (mean, var, max) = cramer_stats(x, λ, d);
            mean + (max * τ) + (max * max * τ * τ + 2.0 * τ * var).sqrt()
        };

        let max_d = (x.ln() + x.ln().ln()).sqrt() * 2.0;
        let d = brentq(|d| err(d) - eps, 0.0, max_d, 0.0, 0.0, 20).unwrap_or(0.0);
        let x1 = (x * (-d * λ).exp()).floor();
        let x2 = (x * (d * λ).exp()).ceil();
        let residue = err(d);
        assert!(residue <= 0.3);

        info!(
            "Δ range = [{:.0}, {:.0}], length = {:.0}, trunc error = {:.6}, d = {:.6}",
            x1,
            x2,
            x2 - x1,
            err(d),
            d,
        );

        (x1 as u64, x2 as u64)
    }

    /// exactly the same as Galway's paper.
    fn plan_Δ_bounds_strict(λ: f64, x: f64, eps: f64) -> (u64, u64) {
        let eps = eps / 2.0;

        let t1 = brentq(|t| err_l(t, x, λ) - eps, 0.0, -8.0, 0.0, 0.0, 100).unwrap_or(x);
        let t2 = brentq(|t| err_r(t, x, λ) - eps, 0.0, 8.0, 0.0, 0.0, 100).unwrap_or(x);

        let x1 = (x * (t1 * λ).exp()).floor();
        let x2 = (x * (t2 * λ).exp()).ceil();

        debug!("t1 = {}, t2 = {}", t1, t2);
        let delta_est = 2.0 * λ * x * (2.0 * (λ * x).ln()).sqrt();
        info!(
            "Δ range = [{}, {}], length = {}, est = {:.0}, residue = ({:.6}, {:.6})",
            x1,
            x2,
            x2 - x1,
            delta_est,
            err_l(t1, x, λ),
            err_r(t2, x, λ)
        );

        (x1.floor() as u64, x2.floor() as u64)
    }

    #[inline(never)]
    fn integrate_critical<T: MyReal>(
        x: u64,
        λ: f64,
        max_order: usize,
        max_height: f64,
        ζ_zeros: impl ZetaZerosDatabase<T>,
    ) -> T {
        let mut integrator = HybridPrecIntegrator::new(
            T::from_u64(x).unwrap(),
            T::one() / 2.0,
            T::mp(λ),
            max_order,
            1e-20,
        );
        let mut result = T::zero();
        let mut last_contribution = T::zero();

        // TODO: check error
        info!("Computing ∑ Φ̂(ρ) for 0 < ℑρ < T.");
        for root in ζ_zeros.map_while(|(x, _)| if x.fp() <= max_height { Some(x) } else { None }) {
            let integral = integrator.query(root).im;
            result += integral;
            last_contribution = integral;
        }

        let max_err = integrator.max_err;
        info!(
            "∑ Φ̂(ρ) = {}, Φ̂(ρ_max) = {:.6e}, max_err/f64::eps = {:.6e}",
            result,
            last_contribution.fp(),
            max_err
        );
        assert!(max_err < 1e13, "possible loss of precision! use PlattIntegrator instead");

        result
    }
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use crate::lmfdb::LMFDB;

    use super::{calc_Δ1_f64, Platt};
    use log::info;
    use F64x2::f64x2;

    #[test]
    fn test_Δ_bounds_heuristic() {
        crate::init();

        let x = 1e11f64;
        let λ = 3e-5;
        let eps = 0.24;
        let segment = 1u64 << 32;

        let (w1, w2) = Platt::plan_Δ_bounds_strict(λ, x, eps);
        let (d1, d2) = Platt::plan_Δ_bounds_heuristic(λ, x, eps);
        info!("   strict bounds = x + [{}, {}]", w1 as i64 - x as i64, w2 - x as u64);
        info!("heuristic bounds = x + [{}, {}]", d1 as i64 - x as i64, d2 - x as u64);

        let mut diff = 0.0;
        if w1 < d1 {
            diff += calc_Δ1_f64(x as u64, eps, λ, w1, d1, segment);
        }
        if d2 < w2 {
            diff += calc_Δ1_f64(x as u64, eps, λ, d2, w2, segment);
        }
        let diff = diff.abs();
        info!("diff = {}", diff);

        assert!(diff < 0.1);
    }

    #[test]
    fn test_integrate_critical() {
        crate::init();

        let x = 1_000_000_000_000_000u64;
        let λ = 3.324516e-7;
        let max_height = 19130220.241455;
        let db = LMFDB::directory("./data/lmfdb");

        let b = Platt::integrate_critical::<f64x2>(x, λ, 15, max_height, db);
        println!("b = {}", b);
        // panic!();
    }
}
