mod bernoulli;
mod complex_functions;
mod erfc;
mod exp_poly_approx;
mod factorial;
mod gabcke_series;
mod gamma;
mod lazy_static;
mod loggamma;
mod riemann_siegel_theta;
mod riemann_siegel_theta_integral;
mod sinc;
mod utils;

pub use bernoulli::Bernoulli;
pub use complex_functions::ComplexFunctions;
pub use erfc::Erfc;
pub use exp_poly_approx::ExpPolyApprox;
pub use factorial::Factorial;
pub use gabcke_series::GabckeSeries;
pub use riemann_siegel_theta::{RiemannSiegelTheta, RiemannSiegelThetaCoeffs};
pub use riemann_siegel_theta_integral::{
    RiemannSiegelThetaIntegral, RiemannSiegelThetaIntegralCoeffs,
};
pub use sinc::Sinc;

pub trait Contexts = Bernoulli + Factorial + Sinc + Erfc + ExpPolyApprox;

pub fn init() {
    use std::sync::Once;
    static ONCE: Once = Once::new();
    ONCE.call_once(|| {
        factorial::init();
        bernoulli::init();
        riemann_siegel_theta::init();
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use F64x2::f64x2;
    use F64x2::test_utils::*;
    type T = f64x2;

    fn is_close(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() / a.abs().max(b.abs()).max(1.0) <= eps
    }

    #[test]
    fn test() {
        init();
        assert!(is_close(f64::factorial(4), 24.0, 1e-9));
        assert!(is_close(f64::bernoulli(4), -0.033333333, 1e-9));
    }

    #[test]
    fn test_rs_theta() {
        init();

        let t = f64x2::new(1000.0, 0.0);

        let eps = 1e-30;
        let gt = f64x2::new(2034.5464280380315, 7.28690383001782e-14);
        let output = t.rs_theta(eps);

        assert_close(output, gt, 0.0, eps);
    }
}
