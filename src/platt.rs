use crate::brentq::brentq;
use crate::platt_integral::PlattIntegrator;
use crate::traits::*;
use log::{debug, info};
use num::integer::*;
use num::{Num, One, Zero};
use std::{
    io::{self, Read},
    marker::PhantomData,
};

#[derive(Default)]
pub struct PlattHints {
    pub λ: Option<f64>,
    pub poly_order: Option<usize>,
}

pub struct Platt {
    λ: f64,
    integral_limit: f64,
    x1: u64,
    x2: u64,
}

impl Default for Platt {
    fn default() -> Self { Self { λ: 0.0, integral_limit: 0.0, x1: 0, x2: 0 } }
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

fn segment(mut x1: u64, x2: u64, seg_size: u64) -> Vec<(u64, u64)> {
    assert!(x1 <= x2);
    let mut segments = vec![];
    loop {
        if x2 - x1 > seg_size {
            segments.push((x1, x1 + seg_size));
            x1 += seg_size;
        } else {
            segments.push((x1, x2));
            break;
        }
    }
    segments
}

#[inline(never)]
fn calc_Δ1_f64(x: u64, eps: f64, λ: f64, x1: u64, x2: u64) -> f64 {
    let approx_n_primes = (x2 - x1) as f64 / (x1 as f64).ln();
    assert!(f64::EPSILON * approx_n_primes < 0.1, "too many primes for f64 approx");

    let mut ϕ = crate::fast_phi::LittlePhiFn::new(λ, eps / (x2 - x1) as f64);
    let c = [1.0 / λ, -1.0 / λ / 2.0, 1.0 / λ / 3.0, -1.0 / λ / 4.0];

    let mut calc = |p: u64| -> f64 {
        let t = (p - x) as i64 as f64 / x as f64;
        if p <= x {
            1.0 - ϕ.query(t)
        } else {
            -ϕ.query(t)
        }
    };

    let mut Δ_1 = 0.0;
    let mut n_primes = 0;

    for (l, r) in segment(x1, x2, 1u64 << 34) {
        info!("sieving [{}, {}]", l, r);
        let sieve_result = crate::sieve::sieve_primesieve(l, r);
        for &p in sieve_result.primes {
            n_primes += 1;
            Δ_1 += calc(p);
        }
    }
    info!("found {} primes in the interval", n_primes);

    ϕ.stat().show("Fast ϕ");

    Δ_1
}

/// let delta = 1/2^53 be the relative differencee by converting x from f64x2 to f64.
/// let r = ln(u/x) / λ
/// Φ(r (1 + delta)) - Φ(r) \approx |Φ'(r)| delta = delta exp(-r^2) / sqrt(2pi) <= delta / 2.
/// note that it's fine if we estimate erfc(x) by
#[inline(never)]
fn calc_Δ_f64(x: u64, eps: f64, λ: f64, x1: u64, x2: u64) -> f64 {
    let mut ret = 0.0;
    // let primes = crate::sieve::linear_sieve(x2.sqrt());
    let primes = crate::sieve::sieve_primesieve(1, x2.sqrt());

    // error analysis: we have ~(x2 - x1)/log(x) many p's.
    // for each p: the error by erfc is delta.
    let Δ_1 = calc_Δ1_f64(x, eps, λ, x1, x2);
    ret += Δ_1;

    // error analysis: each has error delta, we have sqrt(x2) many, so in total is $delta sqrt(x2) << 1$.
    for &p in primes.primes {
        let mut m = 1i64;
        let mut power = p;
        while power < x2 / p {
            m += 1;
            power *= p;
            if power < x1 {
                ret -= 1.0 / m as f64;
            } else {
                // only (x_2^1/2 - x_1^1/2) + (x_2^1/3 - x_1^1/3) + ... many
                // The first term dominates, which is still O((x2 - x1)/sqrt(x)) = O(polylog(n)).
                let r = (power as f64 / x as f64).ln() / λ;
                ret -= Φ(r) / m as f64;
            }
        }
    }
    ret
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

/// exactly the same as Galway's paper.
pub fn plan_Δ_bounds_strict(λ: f64, x: f64, eps: f64) -> (u64, u64) {
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

/// Intuition: Galway bounds the difference by considering all numbers, but we only need prime numbers
/// We simply apply Cramer's random model: only O(1/log(n)) numbers are prime numbers.
/// Also see the comments for cramer_stats for .
fn plan_Δ_bounds_heuristic(λ: f64, x: f64, eps: f64) -> (u64, u64) {
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
        "Δ range = [{:.0}, {:.0}], length = {:.0}, residue = {:.6}, d = {:.6}",
        x1,
        x2,
        x2 - x1,
        err(d),
        d,
    );

    (x1 as u64, x2 as u64)
}

/// Problem: How to bound the abs integral here?
fn integrate_critical<T: MyReal>(x: u64, λ: f64, limit: f64, max_order: usize) -> T {
    let mut integrator = PlattIntegrator::new(
        T::from_u64(x).unwrap(),
        T::one() / 2.0,
        T::from_f64(λ).unwrap(),
        max_order,
        1e-20,
    );
    let mut result = T::zero();
    let mut last_contribution = T::zero();

    let roots = crate::lmfdb::LMFDB_reader::<T>(limit).unwrap();

    info!("integrating phi(1/2+it) N(t) for t = 0 to Inf.");
    let mut int_abs = T::zero();
    for i in 1..roots.len() {
        let a = roots[i - 1];
        let b = roots[i];
        let integral = integrator.query(a, b).im;
        result += integral * i as f64;
        last_contribution = integral * i as f64;
        int_abs += integral.abs() * i as f64;
    }
    info!("integral critical = {}, last = {}, int |f| = {}", result, last_contribution, int_abs,);
    integrator.stat.show("IntegratorCritical");
    // assert!(abs_integral.to_f64().unwrap() <= 1e12);

    result
}

fn integrate_offline<T: MyReal>(x: u64, λ: f64, limit: f64, max_order: usize) -> T {
    // the result = n / log(n) + o(n)
    let mut integrator_offline = PlattIntegrator::<T>::new(
        T::from_u64(x).unwrap(),
        T::one(),
        T::from_f64(λ).unwrap(),
        max_order,
        0.01,
    );
    let result = integrator_offline.query(T::zero(), T::from_f64(limit).unwrap()).im;
    info!("offline integral = {}", result);
    integrator_offline.stat.show("IntegratorOffline");

    result
}

fn plan_integral(λ: f64, x: f64, eps: f64) -> f64 { (x.ln() + x.ln().ln()).sqrt() / λ }

impl Platt {
    pub fn new() -> Self { Self::default() }

    /// During planning, these hyperparameters (λ, sigma, h, x1, x2, integral_limits)
    /// doesn't need to be very accurate
    /// as long as they satisfy the error bound.
    pub fn plan(&mut self, x: f64, hints: PlattHints) {
        let λ = (hints.λ.unwrap_or(50.0) * x.ln() / x).sqrt();
        info!("λ = {:.6e}", λ);

        let (x1, x2) = plan_Δ_bounds_heuristic(λ, x, 0.24);

        let integral_limit = plan_integral(λ, x, 0.1);
        info!("integral limit = {:.6}", integral_limit);

        self.λ = λ;
        self.integral_limit = integral_limit;
        self.x1 = x1;
        self.x2 = x2;
    }

    pub fn compute<T: MyReal>(&mut self, x: u64, hints: PlattHints) -> u64 {
        let max_order = hints.poly_order.unwrap_or(15);
        self.plan(x as f64, hints);

        // this requires high precision: result ≈ # primes
        let integral_offline = integrate_offline::<T>(x, self.λ, self.integral_limit, max_order);

        // maybe this only requires low precision, e.g., f64?
        let integral_critical = integrate_critical::<T>(x, self.λ, self.integral_limit, max_order);

        let Δ = calc_Δ_f64(x, 0.5, self.λ, self.x1, self.x2);
        info!("Δ = {}", Δ);
        let ans = integral_offline - integral_critical * 2.0 - 2.0.ln() + Δ;
        let n_ans = ans.round().to_u64().unwrap();
        info!("ans = {}, rounding error = {}", ans, ans - T::from_u64(n_ans).unwrap());

        n_ans as u64
    }
}

#[cfg(test)]
mod tests {
    use super::integrate_critical;
    use super::{calc_Δ1_f64, cramer_stats, plan_Δ_bounds_heuristic, plan_Δ_bounds_strict};
    use log::info;
    use F64x2::f64x2;

    #[test]
    fn test_Δ_bounds_heuristic() {
        env_logger::init();

        let x = 1e11f64;
        let λ = 3e-5;
        let eps = 0.24;

        let (w1, w2) = plan_Δ_bounds_strict(λ, x, eps);
        let (d1, d2) = plan_Δ_bounds_heuristic(λ, x, eps);
        info!("   strict bounds = x + [{}, {}]", w1 as i64 - x as i64, w2 - x as u64);
        info!("heuristic bounds = x + [{}, {}]", d1 as i64 - x as i64, d2 - x as u64);

        let mut diff = 0.0;
        if w1 < d1 {
            diff += calc_Δ1_f64(x as u64, eps, λ, w1, d1);
        }
        if d2 < w2 {
            diff += calc_Δ1_f64(x as u64, eps, λ, d2, w2);
        }
        let diff = diff.abs();
        info!("diff = {}", diff);

        assert!(diff < 0.1);
    }

    #[test]
    fn test_integrate_critical() {
        env_logger::init();

        let x = 1_000_000_000_000_000u64;
        let λ = 3.324516e-7;
        let limit = 19130220.241455;
        let b = integrate_critical::<f64x2>(x, λ, limit, 15);
        println!("b = {}", b);
        panic!();
    }
}
