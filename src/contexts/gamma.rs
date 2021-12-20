use std::f64::consts::PI;

use num::Complex;

pub trait Gamma {
    fn gamma(&self, eps: f64) -> Self;
}

impl Gamma for f64 {
    #[inline]
    fn gamma(&self, _eps: f64) -> Self { rgsl::gamma_beta::gamma::gamma(*self) }
}

impl Gamma for Complex<f64> {
    /// https://epubs.siam.org/doi/pdf/10.1137/0731050
    /// Spouge's algorithm
    #[inline]
    fn gamma(&self, eps: f64) -> Self {
        let z = *self;
        let z = z - 1.0;
        let mut a = (-z.re).max(2.0).floor() + 0.5;
        loop {
            let err = a.sqrt() / (2.0 * PI).powf(a + 0.5) / (z.re + a);
            if err < eps {
                break;
            }
            a += 1.0;
        }
        let mut coef = Complex::new(1.0, 0.0);
        let mut k = 1.0;
        let mut fac = 1.0;
        while k < a {
            let c_k = (a - k).powf(k - 0.5)
                * (a - k).exp()
                * (if k as i32 % 2 == 1 { 0.5 } else { -0.5 })
                * std::f64::consts::FRAC_2_SQRT_PI
                * std::f64::consts::FRAC_1_SQRT_2
                / fac;
            fac *= k;
            coef += c_k / (z + k);
            k += 1.0;
        }
        // dbg!(coef, (z + a).powc(z + 0.5) / (z + a).exp() * (2.0 * PI).sqrt());
        coef * (z + a).powc(z + 0.5) / (z + a).exp() * (2.0 * PI).sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::Gamma;
    use num::Complex;

    fn _test_gamma(z: Complex<f64>, eps: f64, gt: Complex<f64>) {
        let result = z.gamma(eps);
        let diff = (result - gt).norm();
        println!("{:?} {:?}", result, gt);
        assert!(diff <= eps);
    }

    #[test]
    fn test_gamma() {
        _test_gamma(
            Complex::new(4.0, 10.0),
            1e-10,
            Complex::new(0.0007715342942399662, -0.0010190827990417),
        );
        _test_gamma(Complex::new(4.0, 0.0), 1e-10, Complex::new(6.0, 0.0));
    }
}
