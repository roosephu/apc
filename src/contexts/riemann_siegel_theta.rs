// use super::utils::{impl_from_uninit_cell, mpf_to_f64x2, read_data};
use crate::traits::MyReal;
use F64x2::f64x2;

pub trait RiemannSiegelTheta {
    fn rs_theta(&self, atol: f64) -> Self;
}

pub trait RiemannSiegelThetaCoeffs {
    fn rs_theta_coeff(n: usize) -> Self;
}

const RIEMANN_SIEGEL_THETA_COEFFS_F64: [f64; 31] =
    include!("../../const_table/rs_theta_coeffs_f64.rs");
const RIEMANN_SIEGEL_THETA_COEFFS_F64X2: [f64x2; 31] =
    include!("../../const_table/rs_theta_coeffs_f64x2.rs");

impl RiemannSiegelThetaCoeffs for f64 {
    fn rs_theta_coeff(n: usize) -> Self { RIEMANN_SIEGEL_THETA_COEFFS_F64[n] }
}

impl RiemannSiegelThetaCoeffs for f64x2 {
    fn rs_theta_coeff(n: usize) -> Self { RIEMANN_SIEGEL_THETA_COEFFS_F64X2[n] }
}

impl<T: MyReal + RiemannSiegelThetaCoeffs> RiemannSiegelTheta for T {
    /// See On asymptotic approximations to the log-Gamma and Riemann-Siegel
    /// theta functions, by Richard Brent.
    ///
    /// Cutting off at k-th term has error of $exp(-\pi t) / 2 + \eta_k
    /// \sqrt{\pi t} T_k$. For k <= 20, we have $\eta_k \sqrt{\pi k} ≤ 7.93$.
    #[inline]
    fn rs_theta(&self, atol: f64) -> T {
        let t = *self;

        let eps = atol / 7.93;

        // needs high precision base computation here.
        let mut ret = t / 2.0 * (t / 2.0 / T::PI() / T::E()).ln() - T::FRAC_PI_8();
        if t.fp() <= 25.0 || atol < 1e-33 {
            // When both fails, we simply ignore the $arctan(exp(-\pi t)) / 2$ term,
            // because it's too small ( $< 8 × 10^{-35}$) compared to atol.
            ret += (-t * T::PI()).exp().atan() / 2.0;
        }
        let mut tpow_inv = t.recip();
        let tsqr_inv = tpow_inv * tpow_inv;
        for k in 1..=30 {
            let term = T::rs_theta_coeff(k) * tpow_inv;
            ret += term;

            if term.fp() < eps {
                break;
            }
            tpow_inv *= tsqr_inv;
        }
        ret
    }
}

#[cfg(test)]
mod test {
    use F64x2::test_utils::assert_close;

    use super::*;

    #[test]
    fn rs_theta() {
        crate::init();

        let x = f64x2::mp(101.0);
        let atol = 1e-28;
        let y = x.rs_theta(atol);
        let gt = f64x2::new(89.35830143691956, 6.162759849478704e-15);
        assert_close(y, gt, atol, 0.0);

        let atol = 1e-12;
        let y = x.fp().rs_theta(atol);
        assert_close(y, gt.fp(), atol, 0.0);

        let t = f64x2::new(1000.0, 0.0);
        let eps = 1e-27;
        let gt = f64x2::new(2034.5464280380315, 7.28690383001782e-14);
        let output = t.rs_theta(eps);
        assert_close(output, gt, eps, 0.0);

        let t = f64x2::new(10.0, 0.0);
        let gt = f64x2::new(-3.0670743962898954, 1.2987254263207872e-16);
        let eps = 1e-27;
        let output = t.rs_theta(eps);
        assert_close(output, gt, eps, 0.0);
    }
}
