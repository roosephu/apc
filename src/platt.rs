use crate::brentq::brentq;
use crate::{f64x2, platt_integral::PlattIntegrator};
use crate::{traits::*, unchecked_cast::UncheckedCast};
use log::{debug, info};
use num::integer::*;
use num::{Num, One, Zero};
use std::io::{self, Read};

type T = f64x2;

#[derive(Default)]
pub struct PlattHints {
    pub lambda: Option<f64>,
}

pub struct Platt<T> {
    lambda: T,
    integral_limit: T,
    x1: u64,
    x2: u64,
}

impl<T: MyReal> Default for Platt<T> {
    fn default() -> Self { Self { lambda: T::zero(), integral_limit: T::zero(), x1: 0, x2: 0 } }
}

fn Phi(r: f64) -> f64 { rgsl::error::erfc(r / std::f64::consts::SQRT_2) / 2.0 }

/// let delta = 1/2^53 be the relative differencee by converting x from f64x2 to f64.
/// let r = ln(u/x) / lambda
/// Phi(r (1 + delta)) - Phi(r) \approx |Phi'(r)| delta = delta exp(-r^2) / sqrt(2pi) <= delta / 2.
/// note that it's fine if we estimate erfc(x) by
#[inline(never)]
fn calc_delta_f64(x: u64, eps: f64, lambda: f64, x1: u64, x2: u64) -> f64 {
    let mut ret = 0.0;
    let primes = crate::sieve::linear_sieve(x2.sqrt());

    let c1 = -1.0 / lambda;
    let c2 = -1.0 / lambda / 2.0;
    let c3 = -1.0 / lambda / 3.0;

    let mid1 = (x + x1) / 2;
    let mid2 = (x + x2) / 2;

    // error analysis: we have ~(x2 - x1)/log(x) many p's.
    // for each p: the error by erfc is delta.
    for p in crate::sieve::sieve(&primes, x1, x2) {
        // here we approximate r = ln(u / x) / lambda = ln(1 - (x-u)/x) / lambda.
        // Expad ln(1 - (x-u)/x) at u = x.
        let t = (x as i64 - p as i64) as f64 / x as f64;
        let r = c1 * t + c2 * t * t + c3 * t * t * t;
        let f;
        if p <= x {
            f = 1.0 - Phi(r);
        } else {
            f = -Phi(r);
        }
        ret += f;
    }

    // error analysis: each has error delta, we have sqrt(x2) many, so in total is $delta sqrt(x2) << 1$.
    for p in primes {
        let mut m = 1i64;
        let mut power = p;
        while power < x2 / p {
            m += 1;
            power *= p;
            if power < x1 {
                ret -= 1.0 / m as f64;
            } else {
                // only (x2^1/2 - x1^1/2) + (x2^1/3 - x1^1/3) + ... many
                // The first term dominates, which is still O((x2 - x1)/sqrt(x)) = O(polylog(n)).
                let r = (power as f64 / x as f64).ln() / lambda;
                ret -= Phi(r) / m as f64;
            }
        }
    }
    ret
}

impl<T: MyReal> Platt<T> {
    pub fn new() -> Self { Self::default() }

    #[inline]
    fn Phi(&self, p: T, eps: f64) -> T { (p / T::SQRT_2()).erfc(eps) / 2.0f64 }

    #[inline]
    fn phi(&self, u: T, x: T, eps: f64) -> T { self.Phi((u / x).ln() / self.lambda, eps) }

    /// During planning, these hyperparameters (lambda, sigma, h, x1, x2, integral_limits)
    /// doesn't need to be very accurate
    /// as long as they satisfy the error bound.
    pub fn plan(&mut self, x: f64, hints: PlattHints) {
        let lambda = hints.lambda.unwrap_or(1.0) / x.sqrt();

        let (x1, x2) = self.plan_delta_bounds(lambda, x, 0.24);
        let integral_limit = self.plan_integral(lambda, x, 0.1);
        info!("lambda = {:.6}", lambda);
        let delta_est = 2.0 * lambda * x * (2.0 * (lambda * x).ln()).sqrt();
        info!("delta range = [{}, {}], length = {}, est = {:.0}", x1, x2, x2 - x1, delta_est,);
        info!("integral limit = {:.6}", integral_limit);

        self.lambda = lambda.unchecked_cast();
        self.integral_limit = integral_limit.unchecked_cast();
        self.x1 = x1;
        self.x2 = x2;
    }

    /// we are integrating \hat_\phi(s), which is approximately x^sigma (-\lambda^2 h^2 / 2) with sigma = 0.5 or 1.
    fn plan_integral(&mut self, lambda: f64, x: f64, eps: f64) -> f64 { 6.0 / lambda }

    fn plan_delta_bounds(&mut self, lambda: f64, x: f64, eps: f64) -> (u64, u64) {
        let eps = eps / 2.0;
        let Ep = |u: f64| {
            x * (lambda * lambda / 2.0).exp() * Phi((u / x).ln() / lambda - lambda)
                - u * Phi((u / x).ln() / lambda)
        };
        let Em = |u: f64| {
            u * Phi(-(u / x).ln() / lambda)
                - x * (lambda * lambda / 2.0).exp() * Phi(lambda - (u / x).ln() / lambda)
        };

        let x1 = brentq(|u| Em(u) - eps, 2.0, x, 0.0, 0.0, 100).unwrap_or(x);
        let x2 = brentq(|u| Ep(u) - eps, x * 2.0, x, 0.0, 0.0, 100).unwrap_or(x);

        // let x1 = x - (x - x1) * 0.7;
        // let x2 = x + (x2 - x) * 0.7;
        info!("x1 residue = {}, x2 residue = {}", Em(x1), Ep(x2));

        (x1.floor() as u64, x2.floor() as u64)
    }

    pub fn compute(&mut self, n: u64, hints: PlattHints) -> u64 {
        self.plan(n as f64, hints);

        // the result = n / log(n) + o(n)
        let integral_offline =
            PlattIntegrator::<T>::new(T::from_u64(n).unwrap(), T::one(), self.lambda, 20, 0.01)
                .query(T::zero(), self.integral_limit)
                .im;
        info!("offline integral = {}", integral_offline);

        let mut integrator_critical = PlattIntegrator::<T>::new(
            T::from_u64(n).unwrap(),
            T::one() / 2.0,
            self.lambda,
            20,
            1e-20,
        );
        let mut integral_critical = T::zero();
        let mut last_contribution = T::zero();

        let roots = crate::lmfdb::LMFDB_reader(self.integral_limit).unwrap();

        let mut abs_integral = T::zero();
        for i in 0..roots.len() - 1 {
            let a = roots[i];
            let b = roots[i + 1];
            let integral = integrator_critical.query(a, b).im;
            integral_critical += integral * (i + 1) as f64;
            abs_integral += integral.abs() / (n as f64).sqrt() * (i + 1) as f64;
            last_contribution = integral * (i + 1) as f64;
            // debug!("current zero: {}, integral = {}, est = {}", roots[i], integral, (-self.lambda * self.lambda * a * a / 2.0).exp() * (b - a));
        }
        info!(
            "integral critical = {}, last = {}, abs = {}",
            integral_critical, last_contribution, abs_integral
        );

        let delta = calc_delta_f64(n, 0.5, self.lambda.unchecked_cast(), self.x1, self.x2);
        info!("delta = {}", delta);
        let ans = integral_offline - integral_critical * 2.0 - 2.0.ln() + delta;
        let n_ans = ans.round().unchecked_cast::<i64>();
        info!("ans = {}, residue = {}", ans, ans - n_ans.unchecked_cast::<T>());

        n_ans as u64
    }
}
