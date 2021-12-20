mod bernoulli;
mod complex_functions;
mod erfc;
mod exp_poly_approx;
mod factorial;
mod gabcke_est;
mod gamma;
mod lazy_static;
mod loggamma;
mod riemann_siegel_theta;
mod sinc;
mod utils;

pub use bernoulli::Bernoulli;
pub use complex_functions::ComplexFunctions;
pub use erfc::Erfc;
pub use exp_poly_approx::ExpPolyApprox;
pub use factorial::Factorial;
pub use gabcke_est::GabckeExpansion;
pub use riemann_siegel_theta::{RiemannSiegelTheta, RiemannSiegelThetaCoeffs};
pub use sinc::Sinc;

pub trait Contexts = Bernoulli + Factorial + Sinc + Erfc + ExpPolyApprox;

pub fn init() {
    factorial::init();
    bernoulli::init();
}

#[cfg(test)]
mod tests {
    use super::*;

    fn is_close(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() / a.abs().max(b.abs()).max(1.0) <= eps
    }

    #[test]
    fn test() {
        init();
        assert!(is_close(f64::factorial(4), 24.0, 1e-9));
        assert!(is_close(f64::bernoulli(4), -0.033333333, 1e-9));
    }
}
