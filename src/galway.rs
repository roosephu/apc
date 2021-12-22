use std::f64::consts::PI;

use crate::{brentq::brentq, zeta::FnZeta};
use crate::{context::Context, traits::MyReal};
use log::{debug, info};
use num::integer::*;
use num::Complex;

pub struct Galway<'a, T, Z> {
    ctx: &'a Context<T>,
    fn_zeta: &'a mut Z,
    lambda: T,
    sigma: T,
    integral_limit: T,
    h: T,
    x1: u64,
    x2: u64,
}

impl<T: MyReal, Z> Galway<'_, T, Z> {
    fn Phi(&self, p: T) -> T { (p / T::SQRT_2()) / T::mp(2.0) }

    fn phi(&self, u: T, x: T) -> T { self.Phi((u / x).ln() / self.lambda) }

    fn calc_delta(&self, x: u64, _: f64) -> T {
        let mut ret = T::zero();
        let fx = T::from_u64(x).unwrap();
        let (x1, x2) = (self.x1, self.x2);
        // let eps = eps / ((x2 - x1 + 1) + x2.sqrt() + 1) as f64;

        for &p in crate::sieve::sieve_primesieve(x1, x2).primes {
            ret -= self.phi(T::from_u64(p).unwrap(), fx);
            if p <= x {
                ret += T::one();
            }
        }

        for &p in crate::sieve::sieve_primesieve(1, x2.sqrt()).primes {
            let mut m = 1i64;
            let mut power = p;
            while power < x2 / p {
                m += 1;
                power *= p;
                if power < x1 {
                    ret -= T::from_i64(m).unwrap().recip();
                } else {
                    ret -= self.phi(T::from_u64(power).unwrap(), fx) / T::from_i64(m).unwrap();
                }
            }
        }
        ret
    }
}

impl<T: MyReal, Z> Galway<'_, T, Z> {
    fn init_F_taylor(&mut self, N: usize) {
        // [de Reyna]
        let pi = T::PI();
        let ctx = self.ctx;
        for n in 0..=(N / 2) {
            let mut s1 = Complex::zero();
            for k in 0..=n {
                s1 += T::mp(4.0).pow((n - k) as i32) * ctx.euler(n - k)
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
                    * T::mp(2.0).pow((n - j + 1) as i32);
            }
        }
    }
}

impl<T: MyReal, Z: FnZeta<T>> Galway<'_, T, Z> {
    /// a little bit different from Galway's definition
    /// I divied it by x^sigma.
    fn Psi(&mut self, s: Complex<T>, ln_x: T, eps: f64) -> Complex<T> {
        ((self.lambda * s).powi(2) / T::mp(2.0) + Complex::new(T::zero(), s.im * ln_x)).exp()
            * self.fn_zeta.zeta(s, eps).ln()
            / s
    }

    fn calc_pi_star(&mut self, x: T, eps: f64) -> T {
        let n_total_evals = (self.integral_limit / self.h).ceil().to_i64().unwrap();

        let eps = eps / 4.0 / x.mp().powf(self.sigma.mp()) / x.mp().ln() / n_total_evals as f64;
        let ln_x = x.ln();

        let mut ans = Complex::<T>::zero();

        for t in 1..=n_total_evals {
            let s = Complex::new(self.sigma, self.h * T::from_i64(t).unwrap());
            ans += self.Psi(s, ln_x, eps);
            // println!("s = {}, Psi = {}, eps = {}", s, self.Psi(s, ln_x, eps), eps);
            if t % (n_total_evals / 100).max(1) == 0 || t == n_total_evals {
                debug!(
                    "n total evals = {}, progress = {}, height = {:.6}, ans = {:.16}, Psi = {:.6e}",
                    n_total_evals,
                    t,
                    self.h * T::from_i64(t).unwrap(),
                    ans,
                    self.Psi(s, ln_x, eps)
                );
            }
        }
        // multiply the result by x^sigma, as noted in Psi.
        self.h / T::PI()
            * (self.Psi(Complex::new(self.sigma, T::zero()), ln_x, eps) / T::mp(2.0) + ans).re
            * x.powf(self.sigma)
    }
}

#[derive(Default)]
pub struct GalwayHints {
    pub lambda: Option<f64>,
}

impl<'a, T: MyReal, Z: FnZeta<T>> Galway<'a, T, Z> {
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

    pub fn compute(&mut self, x: u64, hints: GalwayHints) -> u64 {
        self.plan(x as f64, hints);

        let eps = 0.4;
        let pi_star = self.calc_pi_star(T::from_u64(x).unwrap(), eps / 2.0);
        info!("pi^* = {:.6}", pi_star);
        let delta = self.calc_delta(x, eps / 2.0);
        info!("delta = {:.6}", delta);
        info!("sum = {:.6}", pi_star + delta);
        (pi_star + delta).round().to_u64().unwrap()
    }

    /// During planning, these hyperparameters (lambda, sigma, h, x1, x2, integral_limits)
    /// doesn't need to be very accurate
    /// as long as they satisfy the error bound.
    pub fn plan(&mut self, x: f64, hints: GalwayHints) {
        let sigma = 1.5;
        let lambda = hints.lambda.unwrap_or(1.0) / x.sqrt();

        let (x1, x2) = self.plan_delta_bounds(lambda, x, 0.24);
        let h = self.plan_h(sigma, lambda, x, 0.2);
        let integral_limit = self.plan_integral(sigma, lambda, x, 0.1);
        info!("sigma = {:.6}", sigma);
        info!("lambda = {:.6}", lambda);
        info!("h = {:.6}", h);
        info!("delta range = [{}, {}], length = {}", x1, x2, x2 - x1);
        info!(
            "integral limit = {:.6}, # zeta evals = {}",
            integral_limit,
            (integral_limit / h).ceil()
        );

        self.sigma = T::mp(sigma);
        self.lambda = T::mp(lambda);
        self.h = T::mp(h);
        self.integral_limit = T::mp(integral_limit);
        self.x1 = x1;
        self.x2 = x2;
        self.fn_zeta.prepare_multi_eval(self.h, 0.0);
    }

    fn plan_integral(&mut self, sigma: f64, lambda: f64, x: f64, eps: f64) -> f64 {
        let limit = 0.75 * eps
            / ((lambda * lambda * sigma * sigma / 2.0).exp() / (2.0 * PI)
                * rgsl::zeta::riemann::zeta(sigma).ln()
                * x.powf(sigma));
        let u = brentq(
            |x| rgsl::exponential_integrals::E1(x).ln() - limit.ln(),
            lambda * lambda / 2.0,
            -limit.ln(),
            0.0,
            0.0,
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

#[cfg(test)]
mod tests {
    use crate::zeta::FnZeta;
    use crate::*;
    use num::Complex;
    use F64x2::f64x2;

    #[test]
    fn test_tmp() {
        let ctx1 = Context::<f64>::new(100);
        let mut zeta_galway1 = ZetaGalway::new(&ctx1);
        zeta_galway1.prepare_multi_eval(0.23939822958279525, 1e-10);

        let ctx2 = Context::<f64x2>::new(100);
        let mut zeta_galway2 = ZetaGalway::new(&ctx2);
        zeta_galway2.prepare_multi_eval(f64x2 { hi: 0.23939822958279525, lo: 0.0 }, 1e-10);

        let s1 = Complex::<f64>::new(1.5, 0.23939822958279525);
        let s2 = Complex::<f64x2>::new(
            f64x2 { hi: 1.5, lo: 0.0 },
            f64x2 { hi: 0.23939822958279525, lo: 0.0 },
        );
        let result1 = zeta_galway1.zeta(s1, 1e-10);
        let result2 = zeta_galway2.zeta(s2, 1e-10);
        println!("{:?}", result1);
        println!("{:?}", result2);

        // let ctx = Context::<T>::new(100);
        // let mut zeta_galway = ZetaGalway::new(&ctx);
        // let mut galway = Galway::new(&ctx, &mut zeta_galway);
        // let hints = GalwayHints { lambda: opts.lambda };
        // let ans = galway.compute(n, hints);
        // println!("[Galway] ans = {}", ans);
    }
}
