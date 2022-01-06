use super::utils::{impl_from_uninit_cell, mpf_to_f64x2, read_data};
use crate::traits::MyReal;
use F64x2::f64x2;

pub trait RiemannSiegelTheta {
    fn rs_theta(&self, eps: f64) -> Self;
}

pub trait RiemannSiegelThetaCoeffs {
    fn rs_theta_coeff(n: usize) -> Self;
}

impl<T: MyReal + RiemannSiegelThetaCoeffs> RiemannSiegelTheta for T {
    /// See [Sec 3.11, Pugh].
    /// TODO: add test
    #[inline]
    fn rs_theta(&self, eps: f64) -> T {
        let t = *self;

        // as it's typically used with RiemannSiegelZ, we hope it's not too small.
        assert!(t.fp() >= 200.0 && eps > 1e-33);
        const K: usize = 7;

        // needs high precision base computation here.
        let mut ret = t / 2.0 * (t / 2.0 / T::PI() / T::E()).ln() - T::FRAC_PI_8();
        let mut tpow_inv = t.recip();
        let tsqr_inv = tpow_inv * tpow_inv;
        for i in 1..=K {
            ret += T::rs_theta_coeff(i) * tpow_inv;
            tpow_inv *= tsqr_inv;
        }
        ret
    }
}

impl_from_uninit_cell!(RiemannSiegelThetaCoeffs, rs_theta_coeff, f64);
impl_from_uninit_cell!(RiemannSiegelThetaCoeffs, rs_theta_coeff, f64x2);

pub fn init() {
    let data = read_data("data/rs_theta_coeffs.txt", 1000)
        .expect("can't load Bernoulli numbers from `data/rs_theta_coeffs.txt`");

    TABLE_f64.set(data.iter().map(|x| x.to_f64()).collect());
    TABLE_f64x2.set(data.iter().map(mpf_to_f64x2).collect());
}
