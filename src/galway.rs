use std::f64::consts::PI;

use crate::context::Context;
use crate::{brentq::brentq, zeta::FnZeta};
use log::{debug, info};
use num::integer::*;
use num::traits::{FloatConst, Pow, Zero};
use rgsl::error::erfc;
use rustfft::FftPlanner;

type Float = f64;
type Complex = num::Complex<Float>;

pub struct Galway<'a, Z: FnZeta<f64>> {
    ctx: &'a Context<f64>,
    fn_zeta: &'a mut Z,
    lambda: Float,
    sigma: Float,
    integral_limit: Float,
    h: Float,
    x1: u64,
    x2: u64,
}

impl<Z: FnZeta<f64>> Galway<'_, Z> {
    const EXP_APPROX: [Complex; 18] = [
        Complex {
            re: 1.0,
            im: 1.0812144107799266e-23,
        },
        Complex {
            re: -4.479610786345225e-21,
            im: 1.0,
        },
        Complex {
            re: -0.5,
            im: 3.027415987199093e-19,
        },
        Complex {
            re: -8.413828845781633e-18,
            im: -0.16666666666666669,
        },
        Complex {
            re: 0.04166666666666679,
            im: 1.1656663448809618e-16,
        },
        Complex {
            re: -1.0603511152022404e-15,
            im: 0.008333333333332324,
        },
        Complex {
            re: -0.0013888888888827402,
            im: 5.789264486508273e-15,
        },
        Complex {
            re: -2.491995923872859e-14,
            im: -0.00019841269843586228,
        },
        Complex {
            re: 2.4801587374768556e-5,
            im: 6.704175576866034e-14,
        },
        Complex {
            re: -1.594987515102099e-13,
            im: 2.7557317787217356e-6,
        },
        Complex {
            re: -2.755729303110001e-7,
            im: 2.3127460502103687e-13,
        },
        Complex {
            re: -3.2663668749921504e-13,
            im: -2.5052389834713885e-8,
        },
        Complex {
            re: 2.087985316554709e-9,
            im: 2.5867211760028217e-13,
        },
        Complex {
            re: -2.2167241850689593e-13,
            im: 1.6041263496425594e-10,
        },
        Complex {
            re: -1.1352710114429515e-11,
            im: 8.943908448871146e-14,
        },
        Complex {
            re: -4.542339711641447e-14,
            im: -7.962911435347713e-13,
        },
        Complex {
            re: 5.979573239083729e-14,
            im: 7.185782517642856e-15,
        },
        Complex {
            re: -1.970149077208406e-15,
            im: 1.9701490772084063e-15,
        },
    ];

    fn ifft(&self, f: &mut Vec<Complex>) {
        let mut planner = FftPlanner::new();
        let fft = planner.plan_fft_inverse(f.len());
        fft.process(f)
    }

    // compute \sum_{n=1}^N n^{-(s + i n t)} for every 0 <= t < M
    pub fn sum_trunc_dirichlet(
        &self,
        N: usize,
        s: Complex,
        M: usize,
        delta: Float,
    ) -> Vec<Complex> {
        let M2 = (M + M % 2) / 2;
        let s = s + Complex::new(0.0, M2 as f64 * delta);
        let a: Vec<_> = (1..=N)
            .map(|x| Complex::new(x as Float, 0.0).powc(-s))
            .collect();
        let g: Vec<_> = (1..=N).map(|x| (x as Float).ln() * -delta).collect();
        let R = (M + 1).next_power_of_two() * 2;
        let div = 2.0 * PI / R as Float;
        let w: Vec<_> = g.iter().map(|x| ((x / div).round() as i64)).collect();
        let d: Vec<_> = g
            .iter()
            .zip(w.iter())
            .map(|(x, &y)| x - y as f64 * div)
            .collect();

        let mut ret = vec![Complex::zero(); M + 1];
        let mut f = vec![Complex::zero(); R];
        for (e, c) in Self::EXP_APPROX.iter().enumerate() {
            for x in f.iter_mut() {
                x.set_zero()
            }
            for j in 0..N {
                f[w[j].rem_euclid(R as i64) as usize] += a[j] * (d[j] * M2 as f64).pow(e as i32);
            }
            self.ifft(&mut f);
            for x in 0..=M {
                ret[x] += c * f[(x + R - M2) % R] * (x as f64 / M2 as f64 - 1.0).pow(e as i32);
            }
        }

        ret
    }
}

impl<Z: FnZeta<f64>> Galway<'_, Z> {
    fn Phi(&self, p: Float, eps: f64) -> Float {
        erfc(p / Float::SQRT_2()) / 2.0
    }

    fn phi(&self, u: Float, x: Float, eps: f64) -> Float {
        self.Phi((u / x).ln() / self.lambda, eps)
    }

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

    fn calc_delta(&self, x: u64, eps: f64) -> Float {
        let mut ret = Float::zero();
        let (x1, x2) = (self.x1, self.x2);
        let eps = eps / ((x2 - x1 + 1) + x2.sqrt() + 1) as f64;

        let primes = Galway::<Z>::linear_sieve(x2.sqrt());
        for p in Galway::<Z>::sieve(&primes, x1, x2) {
            ret -= self.phi(p as f64, x as f64, eps);
            if p <= x {
                ret += 1.0;
            }
        }

        for p in primes {
            let mut m = 1i32;
            let mut power = p;
            while power < x2 / p {
                m += 1;
                power *= p;
                if power < x1 {
                    ret -= 1. / Float::from(m);
                } else {
                    ret -= self.phi(power as f64, x as f64, eps) / Float::from(m);
                }
            }
        }
        ret
    }
}

impl<Z: FnZeta<f64>> Galway<'_, Z> {
    fn init_F_taylor(&mut self, N: usize) {
        let pi = Float::PI();
        let ctx = self.ctx;
        for n in 0..=(N / 2) {
            let mut s1 = Complex::zero();
            for k in 0..=n {
                s1 += 4.0.pow((n - k) as i32) * ctx.euler(n - k)
                    / ctx.factorial(2 * k)
                    / ctx.factorial(2 * (n - k));
            }
            let mut s2 = Complex::zero();
            for j in 0..=n {
                s2 += (if j % 2 == 0 { 1.0 } else { -1.0 }) * ctx.euler(2 * j)
                    / ctx.factorial(2 * j)
                    * Complex::i().powu((n - j) as u32)
                    * pi.pow((n + j) as i32)
                    / ctx.factorial(n - j)
                    * 2.0.pow((n - j + 1) as i32);
            }
        }
    }
}

impl<Z: FnZeta<f64>> Galway<'_, Z> {
    fn Psi(&mut self, s: Complex, ln_x: f64, eps: f64) -> Complex {
        ((self.lambda * s).powi(2) / 2.0 + s * ln_x).exp() * self.fn_zeta.zeta(s, eps).ln() / s
    }

    fn calc_pi_star(&mut self, x: f64, eps: f64) -> Float {
        let eps = eps / 4.0 / x.powf(self.sigma) / x.ln();
        let ln_x = (x as f64).ln();

        let mut ans = Complex::zero();

        let n_total_evals = (self.integral_limit / self.h).ceil() as u64;
        for t in 1..=n_total_evals {
            let s = Complex::new(self.sigma, self.h * t as f64);
            ans += self.Psi(s, ln_x, eps);
            if t % (n_total_evals / 100).max(1) == 0 || t == n_total_evals {
                info!(
                    "n total evals = {}, progress = {}, height = {:.6}, ans = {}, Psi = {}",
                    n_total_evals,
                    t,
                    self.h * t as f64,
                    ans,
                    self.Psi(s, ln_x, eps)
                );
            }
            assert!(
                ans.is_normal(),
                "ans is not normal! s = {}, Psi = {}, eps = {}",
                s,
                self.Psi(s, ln_x, eps),
                eps
            );
        }
        self.h / PI * (self.Psi(Complex::new(self.sigma, 0.0), ln_x, eps) / 2.0 + ans).re
    }
}

impl<'a, Z: FnZeta<f64>> Galway<'a, Z> {
    pub fn new(ctx: &'a Context<f64>, fn_zeta: &'a mut Z) -> Self {
        Self {
            ctx,
            fn_zeta,
            lambda: 0.0,
            sigma: 0.0,
            x1: 0,
            x2: 0,
            h: 0.0,
            integral_limit: 0.0,
        }
    }

    pub fn compute(&mut self, x: u64) -> i64 {
        let x = x as f64;
        self.plan(x);

        let eps = 0.4;
        let pi_star = self.calc_pi_star(x, eps / 2.0);
        debug!("pi^* = {}", pi_star);
        let delta = self.calc_delta(x as u64, eps / 2.0);
        debug!("delta = {}", delta);
        (pi_star + delta).round() as i64
    }

    fn plan(&mut self, x: Float) {
        self.sigma = 1.5;
        info!("sigma = {:.6}", self.sigma);

        self.lambda = 30.0 / (x as f64).sqrt();
        info!("lambda = {:.6}", self.lambda);

        self.plan_delta_bounds(x, 0.24);
        self.plan_h(x, 0.2);
        self.plan_integral(x, 0.1);
    }

    fn plan_integral(&mut self, x: f64, eps: f64) {
        let sigma = self.sigma;
        let lambda = self.lambda;
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
        self.integral_limit = (2.0 * u).sqrt() / self.lambda;

        info!(
            "integral limit = {:.6}, # zeta evals = {}",
            self.integral_limit,
            (self.integral_limit / self.h).ceil()
        );
    }

    fn plan_h(&mut self, x: f64, eps: Float) {
        let lambda = self.lambda;
        let sigma = self.sigma;
        let h1 = 2.0 * PI * (self.sigma - 1.0) / ((x / eps).ln() + lambda * lambda / 2.0 + 1.0 / x);
        let h2 = 2.0 * PI
            / ((x / 2.0).ln()
                + sigma * lambda * lambda
                + lambda
                    * (2.0 * sigma * (x / 2.0).ln()
                        + sigma * sigma * lambda * lambda
                        + 2.0 * (3.4 / eps).ln())
                    .sqrt());
        self.h = if h1 < h2 { h1 } else { h2 };
        info!("h = {:.6}", self.h);
    }

    fn plan_delta_bounds(&mut self, x: f64, eps: f64) {
        let eps = eps / 2.0;
        let lambda = self.lambda;
        let Ep = |u: f64| {
            x * (lambda * lambda / 2.0).exp() * self.Phi((u / x).ln() / lambda - lambda, eps)
                - u * self.Phi((u / x).ln() / lambda, eps)
        };
        let Em = |u: f64| {
            u * self.Phi(-(u / x).ln() / lambda, eps)
                - x * (lambda * lambda / 2.0).exp() * self.Phi(lambda - (u / x).ln() / lambda, eps)
        };

        let x1;
        if Em(x) > eps {
            x1 = brentq(|u| Em(u) - 0.75 * eps, 2.0, x, 1e-8, 1e-8, 100)
                .unwrap()
                .floor() as u64;
        } else {
            x1 = x as u64;
        }

        let x2;
        if Ep(x) > eps {
            x2 = brentq(|u| Ep(u) - 0.75 * eps, x * 2.0, x, 1e-8, 1e-8, 100)
                .unwrap()
                .floor() as u64;
        } else {
            x2 = x as u64;
        }

        self.x1 = x1;
        self.x2 = x2;
        info!("delta range = [{}, {}], length = {}", x1, x2, x2 - x1);
    }
}

#[cfg(test)]
mod tests {
    use crate::zeta::ZetaGalway;

    use super::*;

    #[test]
    fn test_sum_trunc_dirichlet() {
        let ctx = Context::new(100);
        let mut zeta_galway = ZetaGalway::new(&ctx);
        let galway = Galway::new(&ctx, &mut zeta_galway);
        let N = 20;
        let s = Complex::new(1.3, 2.6);
        let M = 9;
        let delta = 1.1;
        let result = galway.sum_trunc_dirichlet(N, s, M, delta);
        for i in 0..=M {
            let mut sum = Complex::zero();
            let z = s + Complex::new(0.0, i as f64 * delta);
            for j in 1..=N {
                sum += Complex::new(j as f64, 0.0).powc(-z);
            }
            assert!((sum - result[i]).norm() < 2e-15);
        }
    }
}
