use std::f64::consts::PI;

use crate::{brentq::brentq, unchecked_from::UncheckedCast, zeta::FnZeta};
use crate::{context::Context, traits::MyFloat};
use log::{debug, info};
use num::integer::*;
use num::Complex;

pub struct Galway<'a, T, Z: FnZeta<T>> {
    ctx: &'a Context<T>,
    fn_zeta: &'a mut Z,
    lambda: T,
    sigma: T,
    integral_limit: T,
    h: T,
    x1: u64,
    x2: u64,
}

impl<T: MyFloat, Z: FnZeta<T>> Galway<'_, T, Z> {
    fn Phi(&self, p: T, eps: f64) -> T {
        (p / T::SQRT_2()).erfc(eps) / 2.0f64.unchecked_cast::<T>()
    }

    fn phi(&self, u: T, x: T, eps: f64) -> T { self.Phi((u / x).ln() / self.lambda, eps) }

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

    fn sieve(primes: &[u64], l: u64, r: u64) -> Vec<u64> {
        let mut mark = bit_vec::BitVec::from_elem((r - l + 1) as usize, false);
        for &p in primes {
            let mut x = std::cmp::max((l - 1) / p + 1, 2) * p;
            while x <= r {
                mark.set((x - l) as usize, true);
                // might overflow!
                x += p;
            }
        }
        let mut ret = vec![];
        for x in l..=r {
            if !mark[(x - l) as usize] {
                ret.push(x);
            }
        }

        ret
    }

    fn calc_delta(&self, x: u64, eps: f64) -> T {
        let mut ret = T::zero();
        let fx = (x as i64).unchecked_cast::<T>();
        let (x1, x2) = (self.x1, self.x2);
        let eps = eps / ((x2 - x1 + 1) + x2.sqrt() + 1) as f64;

        let primes = Galway::<T, Z>::linear_sieve(x2.sqrt());
        for p in Galway::<T, Z>::sieve(&primes, x1, x2) {
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
                    ret -= m.unchecked_cast::<T>().recip();
                } else {
                    ret -= self.phi(m.unchecked_cast(), fx, eps) / m.unchecked_cast::<T>();
                }
            }
        }
        ret
    }
}

impl<T: MyFloat, Z: FnZeta<T>> Galway<'_, T, Z> {
    fn init_F_taylor(&mut self, N: usize) {
        // [de Reyna]
        let pi = T::PI();
        let ctx = self.ctx;
        for n in 0..=(N / 2) {
            let mut s1 = Complex::zero();
            for k in 0..=n {
                s1 += 4.0f64.unchecked_cast::<T>().pow((n - k) as i32) * ctx.euler(n - k)
                    / ctx.factorial(2 * k)
                    / ctx.factorial(2 * (n - k));
            }
            let mut s2 = Complex::<T>::zero();
            for j in 0..=n {
                s2 += (if j % 2 == 0 { T::one() } else { -T::one() }) * ctx.euler(2 * j)
                    / ctx.factorial(2 * j)
                    * Complex::i().powu((n - j) as u32)
                    * pi.pow((n + j) as i32)
                    / ctx.factorial(n - j)
                    * 2.0f64.unchecked_cast::<T>().pow((n - j + 1) as i32);
            }
        }
    }
}

impl<T: MyFloat, Z: FnZeta<T>> Galway<'_, T, Z> {
    /// a little bit different from Galway's definition
    /// I divied it by x^sigma.
    fn Psi(&mut self, s: Complex<T>, ln_x: T, eps: f64) -> Complex<T> {
        ((self.lambda * s).powi(2) / 2.0f64.unchecked_cast::<T>()
            + Complex::new(T::zero(), s.im * ln_x))
        .exp()
            * self.fn_zeta.zeta(s, eps).ln()
            / s
    }

    fn calc_pi_star(&mut self, x: T, eps: f64) -> T {
        let eps = eps
            / 4.0
            / x.unchecked_cast::<f64>().powf(self.sigma.unchecked_cast())
            / x.unchecked_cast::<f64>().ln();
        let ln_x = x.ln();

        let mut ans = Complex::<T>::zero();

        let n_total_evals: i64 = (self.integral_limit / self.h).ceil().unchecked_cast();
        for t in 1..=n_total_evals {
            let s = Complex::new(self.sigma, self.h * t.unchecked_cast::<T>());
            ans += self.Psi(s, ln_x, eps);
            if t % (n_total_evals / 100).max(1) == 0 || t == n_total_evals {
                info!(
                    "n total evals = {}, progress = {}, height = {:.6}, ans = {:.16}, Psi = {:.6e}",
                    n_total_evals,
                    t,
                    self.h * t.unchecked_cast::<T>(),
                    ans,
                    self.Psi(s, ln_x, eps)
                );
            }
        }
        // multiply the result by x^sigma, as noted in Psi.
        self.h / T::PI()
            * (self.Psi(Complex::new(self.sigma, T::zero()), ln_x, eps)
                / 2.0f64.unchecked_cast::<T>()
                + ans)
                .re
            * x.powf(self.sigma)
    }
}

impl<'a, T: MyFloat, Z: FnZeta<T>> Galway<'a, T, Z> {
    pub fn new(ctx: &'a Context<T>, fn_zeta: &'a mut Z) -> Self {
        Self {
            ctx,
            fn_zeta,
            lambda: T::zero(),
            sigma: T::zero(),
            x1: 0,
            x2: 0,
            h: T::zero(),
            integral_limit: T::zero(),
        }
    }

    pub fn compute(&mut self, x: u64) -> u64 {
        self.plan(x as f64);

        let eps = 0.4;
        let pi_star = self.calc_pi_star((x as i64).unchecked_cast(), eps / 2.0);
        debug!("pi^* = {:.6}", pi_star);
        let delta = self.calc_delta(x, eps / 2.0);
        debug!("delta = {:.6}", delta);
        (pi_star + delta).round().unchecked_cast::<i64>() as u64
    }

    /// During planning, these hyperparameters (lambda, sigma, h, x1, x2, integral_limits)
    /// doesn't need to be very accurate
    /// as long as they satisfy the error bound.
    fn plan(&mut self, x: f64) {
        let sigma = 1.5;
        let lambda = 30.0 / x.sqrt();

        let (x1, x2) = self.plan_delta_bounds(lambda, x, 0.24);
        let h = self.plan_h(sigma, lambda, x, 0.2);
        let integral_limit = self.plan_integral(sigma, lambda, x, 0.1);
        info!("sigma = {:.6}", sigma);
        info!("lambda = {:.6}", lambda);
        info!("h = {:.6}", self.h);
        info!("delta range = [{}, {}], length = {}", x1, x2, x2 - x1);
        info!(
            "integral limit = {:.6}, # zeta evals = {}",
            integral_limit,
            (integral_limit / h).ceil()
        );

        self.sigma = sigma.unchecked_cast();
        self.lambda = lambda.unchecked_cast();
        self.h = h.unchecked_cast();
        self.integral_limit = integral_limit.unchecked_cast();
        self.x1 = x1;
        self.x2 = x2;
        self.fn_zeta.prepare_multi_eval(self.h, 0.0);
    }

    fn plan_integral(&mut self, sigma: f64, lambda: f64, x: f64, eps: f64) -> f64 {
        let limit = 0.75 * eps / (lambda * lambda * sigma * sigma / 2.0) / (2.0 * PI)
            * rgsl::zeta::riemann::zeta(sigma)
            / x.powf(sigma);
        let u = brentq(
            |x| rgsl::exponential_integrals::E1(x) - limit,
            lambda * lambda / 2.0,
            -limit.ln(),
            eps,
            eps,
            100,
        )
        .unwrap();

        (2.0 * u).sqrt() / lambda
    }

    fn plan_h(&mut self, sigma: f64, lambda: f64, x: f64, eps: f64) -> f64 {
        let h1 = 2.0 * PI * (sigma - 1.0) / ((x / eps).ln() + lambda * lambda / 2.0 + 1.0 / x);
        let h2 = 2.0 * PI
            / ((x / 2.0).ln()
                + sigma * lambda * lambda
                + lambda
                    * (2.0 * sigma * (x / 2.0).ln()
                        + sigma * sigma * lambda * lambda
                        + 2.0 * (3.4 / eps).ln())
                    .sqrt());
        if h1 < h2 {
            h1
        } else {
            h2
        }
    }

    fn plan_delta_bounds(&mut self, lambda: f64, x: f64, eps: f64) -> (u64, u64) {
        let eps = eps / 2.0;
        let Phi = |p| rgsl::error::erfc(p / f64::SQRT_2()) / 2.0;
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
}
