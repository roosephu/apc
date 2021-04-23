use crate::brentq::brentq;
use crate::{f64x2, platt_integral::PlattIntegrator};
use crate::{traits::*, unchecked_cast::UncheckedCast};
use byteorder::{LittleEndian, ReadBytesExt};
use log::{debug, info};
use num::integer::*;
use num::{Num, One, Zero};
use std::io::{self, Read};

fn LMFDB_reader<T: MyReal>(limit: T) -> Result<Vec<T>, std::io::Error> {
    let data_files = [
        "./data/zeros/zeros_14.dat",
        "./data/zeros/zeros_5000.dat",
        "./data/zeros/zeros_26000.dat",
        "./data/zeros/zeros_236000.dat",
        "./data/zeros/zeros_446000.dat",
        "./data/zeros/zeros_2546000.dat",
    ];
    info!("Loading zeta zeros up to {}", limit);

    let eps = 2.0.unchecked_cast::<T>().powi(-101);

    let mut ret = vec![];
    for &file_name in data_files.iter() {
        let mut file = std::fs::File::open(file_name)?;
        let n_blocks = file.read_u64::<LittleEndian>()?;

        for b in 0..n_blocks {
            let t0 = file.read_f64::<LittleEndian>()?;
            let t1 = file.read_f64::<LittleEndian>()?;
            let n0 = file.read_u64::<LittleEndian>()?;
            let n1 = file.read_u64::<LittleEndian>()?;
            debug!(
                "[LMFDB] loading {} block {}, from N({}) = {} to N({}) = {}",
                file_name, b, t0, n0, t1, n1
            );

            let t0 = T::from_f64(t0).unwrap();
            let mut z = 0u128;

            for i in n0..n1 {
                let z1 = file.read_u64::<LittleEndian>()? as u128;
                let z2 = file.read_u32::<LittleEndian>()? as u128;
                let z3 = file.read_u8()? as u128;
                z = z + z1 + (z2 << 64) + (z3 << 96);

                let zz = t0 + T::from_u128(z).unwrap() * eps;

                if zz > limit {
                    return Ok(ret);
                }
                // debug!("read zero: {}", z);
                ret.push(zz);
            }
        }
    }
    panic!("Insufficient zeta zeros data")
}

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

impl<T: MyReal> Platt<T> {
    pub fn new() -> Self { Self { lambda: T::zero(), integral_limit: T::zero(), x1: 0, x2: 0 } }

    #[inline]
    fn Phi(&self, p: T, eps: f64) -> T {
        (p / T::SQRT_2()).erfc(eps) / 2.0f64.unchecked_cast::<T>()
    }

    #[inline]
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

        let primes = Self::linear_sieve(x2.sqrt());
        for p in Self::sieve(&primes, x1, x2) {
            ret -= self.phi((p as f64).unchecked_cast(), fx, eps);
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
                    ret -= self.phi((power as f64).unchecked_cast(), fx, eps)
                        / m.unchecked_cast::<T>();
                }
            }
        }
        ret
    }

    /// During planning, these hyperparameters (lambda, sigma, h, x1, x2, integral_limits)
    /// doesn't need to be very accurate
    /// as long as they satisfy the error bound.
    pub fn plan(&mut self, x: f64, hints: PlattHints) {
        let lambda = hints.lambda.unwrap_or(1.0) / x.sqrt();

        let (x1, x2) = self.plan_delta_bounds(lambda, x, 0.24);
        let integral_limit = self.plan_integral(lambda, x, 0.1);
        info!("lambda = {:.6}", lambda);
        info!("delta range = [{}, {}], length = {}", x1, x2, x2 - x1);
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

        let roots = LMFDB_reader(self.integral_limit).unwrap();
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

        let delta = self.calc_delta(n, 0.5);
        info!("delta = {}", delta);
        let ans =
            integral_offline - integral_critical * 2.0 - 2.0.unchecked_cast::<T>().ln() + delta;
        info!("ans = {}", ans);

        ans.round().unchecked_cast::<i64>() as u64
    }
}
