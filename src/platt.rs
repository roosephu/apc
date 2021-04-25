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

impl<T: MyReal> Platt<T> {
    pub fn new() -> Self { Self::default() }

    #[inline]
    fn Phi(&self, p: T, eps: f64) -> T { (p / T::SQRT_2()).erfc(eps) / 2.0f64 }

    #[inline]
    fn phi(&self, u: T, x: T, eps: f64) -> T { self.Phi((u / x).ln() / self.lambda, eps) }

    #[inline(never)]
    fn linear_sieve(n: u64) -> Vec<u64> {
        let mut mark = bit_vec::BitVec::from_elem(n as usize + 1, false);
        let mut primes = vec![];

        for i in 2..=n {
            if !mark[i as usize] {
                primes.push(i);
            }
            for &p in &primes {
                let t = p * i;
                if t > n {
                    break;
                }
                mark.set(t as usize, true);
                if i as u64 % p == 0 {
                    break;
                }
            }
        }

        primes
    }

    #[inline(never)]
    fn sieve(primes: &[u64], l: u64, r: u64) -> Vec<u64> {
        let ub = (r - l) as usize;
        let mut mark = bit_vec::BitVec::from_elem(ub + 1, false);
        for &p in primes {
            let x = std::cmp::max((l - 1) / p + 1, 2) * p;
            let mut y = (x - l) as usize;
            while y <= ub {
                mark.set(y, true);
                y += p as usize;
            }
        }
        let mut ret = vec![];
        for (idx, &s) in mark.storage().iter().enumerate() {
            let mut s = !s;
            while s != 0 {
                let w = s & (s - 1);
                let offset = (s ^ w).trailing_zeros();
                let x = l + ((idx as u64) << 5) + offset as u64;
                ret.push(x);
                s = w;
            }
        }

        ret
    }

    #[inline(never)]
    fn calc_delta(&self, x: u64, eps: f64) -> T {
        let mut ret = T::zero();
        let fx = (x as i64).unchecked_cast::<T>();
        let (x1, x2) = (self.x1, self.x2);
        let eps = eps / ((x2 - x1 + 1) + x2.sqrt() + 1) as f64;

        let primes = Self::linear_sieve(x2.sqrt());
        for p in Self::sieve(&primes, x1, x2) {
            ret -= self.phi((p as i64).unchecked_cast(), fx, eps);
            if p <= x {
                ret += T::one();
            }
        }

        for p in primes {
            let mut m = 1i64;
            let mut power = p;
            while power < x2 / p {
                m += 1;
                power *= p;
                if power < x1 {
                    ret -= m.unchecked_cast::<T>();
                } else {
                    ret -= self.phi((power as i64).unchecked_cast(), fx, eps)
                        / m.unchecked_cast::<T>();
                }
            }
        }
        ret
    }

    /// let delta = 1/2^53 be the relative differencee by converting x from f64x2 to f64.
    /// let r = ln(u/x) / lambda
    /// Phi(r (1 + delta)) - Phi(r) \approx |Phi'(r)| delta = delta exp(-r^2) / sqrt(2pi) <= delta / 2.
    /// note that it's fine if we estimate erfc(x) by
    #[inline(never)]
    fn calc_delta_f64(&self, x: u64, eps: f64) -> T {
        let mut ret = 0.0;
        let fx = x as f64;
        let (x1, x2) = (self.x1, self.x2);
        let eps = eps / ((x2 - x1 + 1) + x2.sqrt() + 1) as f64;
        info!("delta eps = {:.e}", eps);
        let lambda: f64 = self.lambda.unchecked_cast();

        fn Phi(r: f64) -> f64 { rgsl::error::erfc(r / std::f64::consts::SQRT_2) / 2.0 }

        let w = lambda * x as f64;
        let c1 = -1.0 / (lambda * w);
        let c2 = -0.5 / (lambda * w * w);
        let c3 = -1.0 / 3.0 / (lambda * w * w * w);

        // error analysis: we have ~(x2 - x1)/log(x) many p's.
        // for each p: the error by erfc is delta.
        let primes = Self::linear_sieve(x2.sqrt());
        for p in Self::sieve(&primes, x1, x2) {
            // here we approximate r = ln(u / x) / lambda.
            // Define t = lambda (x - u), and we have r = ln(1 - t/(x lambda)) / lambda.
            // typically t is not large, while x lambda is very large
            // Expad ln(1 - t / (x lambda)) at t = 0.
            // since x - p = O(sqrt(x)) so there is no precision loss.
            let t = (x as i64 - p as i64) as f64 * lambda;
            let r = c1 * t + c2 * t * t + c3 * t * t * t;
            ret -= Phi(r);
            if p <= x {
                ret += 1.0;
            }
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
                    // The first term dominates, which is still O((x2 - x1)/sqrt(x)) = O(1).
                    let r = ((power as f64) / fx).ln() / lambda;
                    ret -= Phi(r) / m as f64;
                }
            }
        }
        ret.unchecked_cast()
    }

    /// During planning, these hyperparameters (lambda, sigma, h, x1, x2, integral_limits)
    /// doesn't need to be very accurate
    /// as long as they satisfy the error bound.
    pub fn plan(&mut self, x: f64, hints: PlattHints) {
        let lambda = hints.lambda.unwrap_or(1.0) / x.sqrt();

        let (x1, x2) = self.plan_delta_bounds(lambda, x, 0.24);
        let integral_limit = self.plan_integral(lambda, x, 0.1);
        info!("lambda = {:.6}", lambda);
        info!(
            "delta range = [{}, {}], length = {}, est = {:.0}",
            x1,
            x2,
            x2 - x1,
            2.0 * lambda * x * (2.0 * (lambda * x).ln()).sqrt()
        );
        info!("integral limit = {:.6}", integral_limit,);

        self.lambda = lambda.unchecked_cast();
        self.integral_limit = integral_limit.unchecked_cast();
        self.x1 = x1;
        self.x2 = x2;
    }

    /// we are integrating \hat_\phi(s), which is approximately x^sigma (-\lambda^2 h^2 / 2) with sigma = 0.5 or 1.
    fn plan_integral(&mut self, lambda: f64, x: f64, eps: f64) -> f64 { 6.0 / lambda }

    fn plan_delta_bounds(&mut self, lambda: f64, x: f64, eps: f64) -> (u64, u64) {
        let eps = eps / 2.0;
        let Phi = |p| rgsl::error::erfc(p / std::f64::consts::SQRT_2) / 2.0;
        let Ep = |u: f64| {
            x * (lambda * lambda / 2.0).exp() * Phi((u / x).ln() / lambda - lambda)
                - u * Phi((u / x).ln() / lambda)
        };
        let Em = |u: f64| {
            u * Phi(-(u / x).ln() / lambda)
                - x * (lambda * lambda / 2.0).exp() * Phi(lambda - (u / x).ln() / lambda)
        };

        let x1 = brentq(|u| Em(u) - 0.75 * eps, 2.0, x, 0.0, 0.0, 100).unwrap_or(x);
        let x2 = brentq(|u| Ep(u) - 0.75 * eps, x * 2.0, x, 0.0, 0.0, 100).unwrap_or(x);

        (x1.floor() as u64, x2.floor() as u64)
    }

    pub fn compute(&mut self, n: u64, hints: PlattHints) -> u64 {
        self.plan(n as f64, hints);

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
        info!("largest zeta roots {}, # zeros = {}", roots.last().unwrap(), roots.len());

        for i in 0..roots.len() - 1 {
            let a = roots[i];
            let b = roots[i + 1];
            let integral = integrator_critical.query(a, b).im;
            integral_critical += integral * (i + 1) as f64;
            last_contribution = integral * (i + 1) as f64;
            // debug!("current zero: {}, integral = {}, est = {}", roots[i], integral, (-self.lambda * self.lambda * a * a / 2.0).exp() * (b - a));
        }
        info!("integral critical = {}, last = {}", integral_critical, last_contribution);

        let delta = self.calc_delta_f64(n, 0.5);
        info!("delta = {}", delta);
        let ans =
            integral_offline - integral_critical * 2.0 - 2.0.unchecked_cast::<T>().ln() + delta;
        info!("ans = {}", ans);

        ans.round().unchecked_cast::<i64>() as u64
    }
}
