use num::Complex;

use crate::sum_trunc_dirichlet::sum_trunc_dirichlet;
use crate::{
    traits::{MyReal, Sinc},
    unchecked_cast::UncheckedCast,
};

type HParams<T> = (T, T);

/// a data structure for querying \sum_{t=1}^k n^{-sigma - i t}
struct BandwidthInterp<T> {
    k: usize,
    tau: T,
    sigma: T,
    t0: T,

    alpha: T,
    beta: T,
    h_params: HParams<T>,
    data: Vec<Complex<T>>,
}

impl<T: MyReal + Sinc> BandwidthInterp<T> {
    pub fn new(k: usize, sigma: T) -> Self {
        let k0 = T::one();
        let k1 = (k as f64).unchecked_cast::<T>();
        let tau = (k1 / k0).ln();
        let beta = tau + 20.0;
        let eps = (beta - tau) / 2.0;
        let alpha = (k0 * k1).ln() / 2.0;
        let c = 80.0f64.unchecked_cast::<T>();
        let min_t = T::PI() * 2.0 * k1 * k1;
        let max_t = T::PI() * 2.0 * (k1 + 1.0) * (k1 + 1.0);
        let t0 = min_t - c / eps;
        let t1 = max_t + c / eps;
        let delta = T::PI() / beta;
        let m = ((t1 - t0) / delta).unchecked_cast::<i64>() as usize;
        let data = sum_trunc_dirichlet(Complex::new(sigma, t0), k, m, delta);
        let data = data
            .iter()
            .enumerate()
            .map(|(i, &x)| Complex::new(T::zero(), (t0 + delta * i as f64) * alpha).exp() * x)
            .collect();
        let coeff = c / c.sinh();
        println!("coef = {}, sinh = {}", coeff, c.sinh());

        println!("precompute {} terms", m);

        Self { k, tau, sigma, alpha, beta, data, h_params: (c, eps), t0 }
    }

    fn h(&self, t: T) -> T {
        let (c, eps) = self.h_params;
        let w = (c * c - eps * eps * t * t).sqrt();
        let coeff = c / c.sinh();
        if w.is_zero() {
            coeff
        } else {
            coeff * w.sinh() / w
        }
    }

    /// we have all precomputed data of
    /// G1(dt) = G(dt + t0) = e^{i alpha (t0 + dt)} \sum_{n=k0}^{k1} n^{-sigma - i(t0 + dt)}
    /// for dt multipliers of delta = pi/beta
    pub fn query(&self, t: T) -> Complex<T> {
        let dt = t - self.t0;
        let delta = T::PI() / self.beta;
        let (c, eps) = self.h_params;
        let r = ((dt + c / eps) / delta).floor();
        let l = ((dt - c / eps) / delta).ceil();

        let mut ret = Complex::<T>::zero();
        for a in l.unchecked_cast::<i32>()..=r.unchecked_cast::<i32>() {
            ret += self.data[a as usize]
                * self.h(dt - delta * (a as f64))
                * (self.beta * dt - T::PI() * (a as f64)).sinc();
            // println!("{} {}", self.h(dt - delta * (a as f64)), (self.beta * dt - T::PI() * (a as f64)).sinc());
        }
        println!(
            "# interp terms = {:}, l = {}, r = {}",
            (r - l + 1.0).unchecked_cast::<i32>(),
            l.unchecked_cast::<i32>(),
            r.unchecked_cast::<i32>()
        );

        ret * Complex::new(T::zero(), -self.alpha * t).exp()
    }
}

#[cfg(test)]
mod tests {
    use super::BandwidthInterp;
    use crate::f64x2;
    use num::traits::FloatConst;
    use num::{Complex, Zero};

    type T = f64x2;

    #[test]
    fn test_bandwidth_limited_interp() {
        let k = 100;
        let sigma = T::from(0.5);
        let ds = BandwidthInterp::<T>::new(k, sigma);
        let t = T::PI() * 2.0 * 10100.0;
        // let t = ds.t0;
        let mut gt = Complex::<T>::zero();
        let s = Complex::new(sigma, t);
        for i in 1..=k {
            gt += Complex::new(T::from(i as f64), T::zero()).powc(-s);
        }
        println!("ds params: alpha = {}, beta = {}, tau = {}", ds.alpha, ds.beta, ds.tau);

        let output = ds.query(t);
        // let output = ds.data[0];
        // gt *= Complex::new(T::zero(), ds.alpha * t).exp();
        println!("gt = {}, output = {}, diff = {:.e}", gt, output, gt - output);
        panic!();
    }
}
