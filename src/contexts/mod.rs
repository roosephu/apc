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
// mod riemann_siegel_theta_integral;
mod sinc;
// mod utils;

pub use bernoulli::Bernoulli;
pub use complex_functions::ComplexFunctions;
pub use erfc::Erfc;
pub use exp_poly_approx::ExpPolyApprox;
pub use factorial::Factorial;
pub use gabcke_series::GabckeSeries;
pub use riemann_siegel_theta::{RiemannSiegelTheta, RiemannSiegelThetaCoeffs};
// pub use riemann_siegel_theta_integral::{
//     RiemannSiegelThetaIntegral, RiemannSiegelThetaIntegralCoeffs,
// };
pub use sinc::Sinc;

pub trait Contexts = Bernoulli + Factorial + Sinc + Erfc + ExpPolyApprox;

#[cfg(test)]
mod tests {
    use super::*;
    use F64x2::f64x2;
    use F64x2::test_utils::*;
    type T = f64x2;

    #[test]
    fn test() {
        crate::init();
        assert_close(f64::factorial(4), 24.0, 1e-9, 0.0);
        assert_close(f64::bernoulli(4), -0.033333333, 1e-9, 0.0);
    }
}
