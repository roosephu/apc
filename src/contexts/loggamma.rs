use crate::traits::MyReal;
use num::Complex;

use super::{Bernoulli, ComplexFunctions};

pub trait Loggamma {
    fn loggamma(&self, eps: f64) -> Self;
}

/// error estimation
/// See [here](https://www.wikiwand.com/en/Stirling%27s_approximation#Error_bounds) for more details.
fn loggamma_err(ln_z: Complex<f64>, n: usize) -> f64 {
    let arg = ln_z.im;
    let norm = ln_z.re.exp();
    let err_coef = if arg < std::f64::consts::FRAC_PI_4 {
        1.0
    } else if arg < std::f64::consts::FRAC_PI_2 {
        1.0 / arg.sin().abs()
    } else {
        panic!("you should normalize z first!, z = {:?}", ln_z);
        // 1.0 / (arg / 2.0).cos().pow(2 * n as i32)
    };
    f64::bernoulli(n * 2).abs() / ((2 * n) * (2 * n - 1)) as f64 / norm.powi((2 * n - 1) as i32)
        * err_coef
}

impl<T: MyReal + Bernoulli> Loggamma for Complex<T>
where
    Complex<T>: ComplexFunctions,
{
    /// log Gamma function by Stirling series
    ///
    fn loggamma(&self, eps: f64) -> Complex<T> {
        const N: usize = 20;
        let mut z = *self;

        assert!(z.re > T::from_f64(-20.0f64).unwrap(), "beyond impl {:?}", z);
        let mut result = Complex::zero();
        while z.re < T::from_usize(N).unwrap() {
            result -= z.ln();
            z += T::one();
        }

        let ln_z = z.ln();

        result += (z - T::from_f64(0.5).unwrap()) * ln_z - z + (T::PI() * 2.0).ln() / 2.0;
        let z2 = z * z;
        let mut zpow = z;
        for i in 1..N {
            let contrib = T::bernoulli(i * 2) / ((2 * i) * (2 * i - 1)) as f64 / zpow;
            result += contrib;

            zpow *= z2;
            if contrib.l1_norm().to_f64().unwrap() * 10.0 < eps {
                break;
            }
        }
        let err = loggamma_err(ln_z.approx(), N);
        assert!(err <= eps);

        result
    }
}
