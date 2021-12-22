use num::{Complex, Signed};

use crate::contexts::{ExpPolyApprox, Sinc};
use crate::sum_trunc_dirichlet::sum_trunc_dirichlet;
use crate::traits::MyReal;
use log::debug;

/// a data structure for querying \sum_{t=1}^k n^{-sigma - i t}
pub struct BandwidthInterp<T> {
    k0: usize,
    k1: usize,
    tau: T,
    sigma: T,
    t0: T,
    gap: T,

    alpha: T,
    beta: T,
    delta: T,
    data: Vec<Complex<T>>,
}

/// determine c: c / sinh(c) = eps
fn find_c(eps: f64) -> f64 {
    let c = 1.0;
    let c = (2.0 / eps * c).ln();
    let c = (2.0 / eps * c).ln();
    c
}

impl<T: MyReal + Sinc + ExpPolyApprox + Signed> BandwidthInterp<T> {
    /// tau = ln(k1/k0)
    /// precompute = O((c + k1 eps)(tau/gap + 2))
    /// query = O(k0 + c (tau/gap + 1))
    pub fn new(k: usize, sigma: T) -> Self {
        // let k0_int = std::cmp::min(4, k);
        let k0_int = 1;
        let k0 = T::mp(k0_int as f64);
        let k1 = T::mp(k as f64);
        let tau = (k1 / k0).ln();
        let gap = tau * 0.5;
        let beta = tau + gap * 2.0;
        let alpha = (k0 * k1).ln() / 2.0;
        let max_c = T::mp(80.0);
        let min_t = T::PI() * 2.0 * k1 * k1;
        let max_t = T::PI() * 2.0 * (k1 + 1.0) * (k1 + 1.0);
        let t0 = min_t - max_c / gap;
        let t1 = max_t + max_c / gap;
        let delta = T::PI() / beta;
        let m = ((t1 - t0) / delta).to_usize().unwrap();
        let data = sum_trunc_dirichlet(Complex::new(sigma, t0), k0_int, k, m, delta);
        let data = data
            .iter()
            .enumerate()
            .map(|(i, &x)| Complex::new(T::zero(), (t0 + delta * i as f64) * alpha).exp() * x)
            .collect();
        debug!("precompute {} terms", m);

        Self { k0: k0_int, k1: k, tau, sigma, alpha, beta, data, gap, t0, delta }
    }

    #[inline]
    fn h(&self, c: T, t: T) -> T {
        let w = (c * c - self.gap * self.gap * t * t).sqrt();
        if w.is_zero() {
            T::one()
        } else {
            w.sinh() / w
        }
    }

    /// we have all precomputed data of
    /// G1(dt) = G(dt + t0) = e^{i alpha (t0 + dt)} \sum_{n=k0}^{k1} n^{-sigma - i(t0 + dt)}
    /// for dt multipliers of delta = pi/beta
    pub fn query(&self, t: T, eps: f64) -> Complex<T> {
        let c = T::mp(find_c(eps / 10.0));
        let c_over_c_sinh = c / c.sinh();

        let dt = t - self.t0;
        let r = ((dt + c / self.gap) / self.delta).floor().fp() as usize;
        let l = ((dt - c / self.gap) / self.delta).ceil().fp() as usize;

        let mut ret = Complex::<T>::zero();
        for a in l..=r {
            ret += self.data[a]
                * self.h(c, dt - self.delta * (a as f64))
                * (self.beta * dt * T::FRAC_1_PI() - a as f64).sinc();
            // println!("{} {}", self.h(c, dt - delta * (a as f64)), (self.beta * dt - T::PI() * (a as f64)).sinc());
        }
        ret *= Complex::new(T::zero(), -self.alpha * t).exp() * c_over_c_sinh;

        let s = Complex::new(self.sigma, t);
        for n in 1..self.k0 {
            ret += (-s * T::from_i32(n as i32).unwrap().ln()).exp();
        }
        // println!(
        //     "c = {}, coeff = {}, # interp terms = {:}, l = {}, r = {}",
        //     c,
        //     c_over_c_sinh,
        //     (r - l + 1.0).to_i32().unwrap(),
        //     l.to_i32().unwrap(),
        //     r.to_i32().unwrap(),
        // );

        ret
    }
}

#[cfg(test)]
mod tests {
    use super::BandwidthInterp;
    use crate::traits::MyReal;
    use num::Complex;
    use F64x2::f64x2;
    use F64x2::test_utils::*;

    #[test]
    fn test_bandwidth_limited_interp() {
        type T = f64x2;

        let k = 100;
        let sigma = T::mp(0.5);
        let eps = 1e-26;
        let ds = BandwidthInterp::<T>::new(k, sigma);
        let t = T::PI() * 2.0 * 10100.0;
        // let t = ds.t0;
        let mut gt = Complex::<T>::zero();
        let s = Complex::new(sigma, t);
        for i in 1..=k {
            gt += (-s * T::mp(i as f64).ln()).exp();
        }
        println!("ds params: alpha = {}, beta = {}, tau = {}", ds.alpha, ds.beta, ds.tau);

        let output = ds.query(t, eps);
        // let output = ds.data[0];
        // gt *= Complex::new(T::zero(), ds.alpha * t).exp();
        println!("gt = {}, output = {}, diff = {:.e}", gt, output, gt - output);
        assert_complex_close(output, gt, eps);
    }

    #[test]
    fn test_bandwidth_limited_interp_f64() {
        type T = f64;

        let k = 100;
        let sigma = T::mp(0.5);
        let eps = 1e-10;
        let ds = BandwidthInterp::<T>::new(k, sigma);
        let t = T::PI() * 2.0 * 10100.0;
        let mut gt = Complex::<T>::zero();
        let s = Complex::new(sigma, t);
        for i in 1..=k {
            gt += (-s * T::mp(i as f64).ln()).exp();
        }
        println!("ds params: alpha = {}, beta = {}, tau = {}", ds.alpha, ds.beta, ds.tau);

        let output = ds.query(t, eps);
        println!("gt = {}, output = {}, diff = {:.e}", gt, output, gt - output);
        // (a - b).norm().approx() / b.norm().approx();
        let diff = (output - gt).norm() / gt.norm();
        assert!(diff <= eps);
        // assert_complex_close(output, gt, eps);
    }
}
