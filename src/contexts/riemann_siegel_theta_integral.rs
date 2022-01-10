// use super::utils::{impl_from_uninit_cell, mpf_to_f64x2, read_data};
use crate::traits::MyReal;
use F64x2::f64x2;

pub trait RiemannSiegelThetaIntegral {
    fn rs_theta_integral(&self, eps: f64) -> Self;
}

pub trait RiemannSiegelThetaIntegralCoeffs {
    fn rs_theta_integral_coeff(n: usize) -> Self;
}

impl<T: MyReal + RiemannSiegelThetaIntegralCoeffs> RiemannSiegelThetaIntegral for T {
    /// Use the expansion: $\theta(t) = \frac{1}{2} t \ln t - \frac{\ln 2 \pi e}{2} t - \frac{\pi}{8} + ...$
    /// So $\int \theta(t) \mathrm{d} t = \frac{1}{2} t^2 \ln t - \frac{1}{4} (2 + \ln 2 \pi) t^2 - \frac{\pi}{8} t + ...$
    #[inline]
    fn rs_theta_integral(&self, eps: f64) -> T {
        let t = *self;
        assert!(t.fp() >= 200.0 && eps > 1e-33);
        const K: usize = 7;
        let quadratic_coeff = ((T::PI() * 2.0).ln() + 2.0) * 0.25; // TODO: precompute a constant

        let ln_t = t.ln();
        let mut ret = t * t * (ln_t * 0.5 - quadratic_coeff) - T::FRAC_PI_8() * t;
        let t_inv = t.recip();
        let tsqr_inv = t_inv * t_inv;
        let mut tpow_inv = tsqr_inv;
        for i in 1..=K {
            let coeff = T::rs_theta_integral_coeff(i);
            if i == 1 {
                ret += coeff * ln_t;
            } else {
                ret += coeff * tpow_inv;
                tpow_inv *= tsqr_inv;
            }
        }
        ret
    }
}

impl_from_uninit_cell!(RiemannSiegelThetaIntegralCoeffs, rs_theta_integral_coeff, f64);
impl_from_uninit_cell!(RiemannSiegelThetaIntegralCoeffs, rs_theta_integral_coeff, f64x2);

pub fn init() {
    let data = read_data("data/rs_theta_integral_coeffs.txt", 1000)
        .expect("can't load Bernoulli numbers from `data/rs_theta_integral_coeffs.txt`");

    TABLE_f64.set(data.iter().map(|x| x.to_f64()).collect());
    TABLE_f64x2.set(data.iter().map(mpf_to_f64x2).collect());
}
