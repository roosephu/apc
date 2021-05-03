#![allow(non_snake_case)]
#![allow(dead_code)]
#![allow(unused_variables)]
#![allow(unused_imports)]
#![allow(clippy::float_cmp)]
#![allow(clippy::needless_range_loop)]
#![allow(clippy::let_and_return)]
#![allow(clippy::approx_constant)]
#![feature(trait_alias)]
#![feature(destructuring_assignment)]
#![feature(min_type_alias_impl_trait)]
#![feature(non_ascii_idents)]
#![feature(try_blocks)]

mod adaptive_interp;
pub mod bandwidth_interp;
pub mod brentq;
pub mod constants;
pub mod context;
pub mod galway;
pub mod gamma;
mod lmfdb;
pub mod platt;
mod platt_integral;
pub mod power_series;
pub mod profiler;
pub mod riemann_siegel;
mod sieve;
pub mod sum_trunc_dirichlet;
pub(crate) mod test_utils;
pub mod traits;
pub mod zeta;

pub use context::Context;
pub use galway::{Galway, GalwayHints};
pub use riemann_siegel::{RiemannSiegelZ, RiemannSiegelZeta};
pub use zeta::ZetaGalway;

#[cfg(test)]
mod tests {
    use crate::test_utils::*;
    use crate::traits::MyReal;
    use num::Complex;
    use F64x2::f64x2;

    #[test]
    fn test_log() {
        let s = f64x2 { hi: 123.0, lo: 0.0 };
        let x = s.ln();
        let gt = f64x2 { hi: 4.812184355372417, lo: 4.291407929980309e-16 };
        assert_close(x, gt, 1e-30);
    }

    #[test]
    fn test_sqrt() {
        let s = f64x2 { hi: 123.0, lo: 0.0 };
        assert_close(s.sqrt(), f64x2 { hi: 11.090536506409418, lo: -5.209651269937913e-16 }, 1e-30);
    }

    #[test]
    fn test_exp() {
        let s = f64x2 { hi: 3.4538776394910684, lo: 0.0000000000000001184757763427252 };
        assert_close(s.exp(), f64x2 { hi: 31.622776601683793, lo: 7.566535620287156e-16 }, 1e-30);
    }

    #[test]
    fn test_complex_log() {
        let s1 = Complex::new(f64x2 { hi: 1.5, lo: 0.0 }, f64x2 { hi: 10.0, lo: 0.0 });
        let gt = Complex::new(
            f64x2 { hi: 2.3137103974614557, lo: -1.1772777930167866e-16 },
            f64x2 { hi: 1.4219063791853994, lo: -4.7442629531916207e-17 },
        );
        assert_complex_close(s1.ln(), gt, 1e-30);
    }

    #[test]
    fn test_cos_sin() {
        let t = f64x2 { hi: 1.034702354811904, lo: 8.06154861650416e-17 };
        assert_close(t.cos(), f64x2 { hi: 0.5107818439368557, lo: 1.0134834604058354e-17 }, 3e-32);

        let s = f64x2 { hi: 23.025850929940457, lo: -0.00000000000000039439938398199903 };
        assert_close(
            s.cos(),
            f64x2 { hi: -0.5107818439368557, lo: -1.0134834604058354e-17 },
            3e-32,
        );
        assert_close(s.sin(), f64x2 { hi: -0.8597103627992777, lo: 4.575936712504712e-17 }, 3e-32);
    }

    #[test]
    fn test_complex_exp() {
        let s1 = Complex::new(f64x2 { hi: 1.5, lo: 0.0 }, f64x2 { hi: 10.0, lo: 0.0 });
        let gt = Complex::new(
            f64x2 { hi: -3.7604577010937845, lo: -1.9567327806486033e-16 },
            f64x2 { hi: -2.438133466706061, lo: -5.786232568383162e-17 },
        );
        assert_complex_close(s1.exp(), gt, 2e-31);

        let s2 = Complex {
            re: f64x2 { hi: 3.4538776394910684, lo: 0.0000000000000001184757763427252 },
            im: f64x2 { hi: 23.025850929940457, lo: -0.00000000000000039439938398199903 },
        };
        let gt = Complex::new(
            f64x2 { hi: -16.1523401430113, lo: -1.3105798758488752e-15 },
            f64x2 { hi: -27.186428744954082, lo: -8.416517329489051e-16 },
        );
        assert_complex_close(s2.exp(), gt, 3e-32);
    }

    #[test]
    fn test_powc() {
        let s = Complex::new(f64x2 { hi: 1.5, lo: 0.0 }, f64x2 { hi: 1.0, lo: 0.0 });
        let b = Complex::new(f64x2 { hi: 10.0, lo: 0.0 }, f64x2::zero());
        let x = b.powc(s);
        let gt = Complex::new(
            f64x2 { hi: -21.130387081656004, lo: 5.226753154184954e-16 },
            f64x2 { hi: 23.52672399165224, lo: -4.0310962858867735e-17 },
        );
        assert_complex_close(x, gt, 1e-31);
    }

    /// can test: add, power for complex numbers.
    #[test]
    fn test_truncated_dirichlet_series() {
        let n = 1000;
        let s = Complex::new(f64x2 { hi: 1.5, lo: 0.0 }, f64x2 { hi: 1.0, lo: 0.0 });
        let mut ans: Complex<f64x2> = Complex::zero();
        for i in 1..=n {
            let x = Complex::new(f64x2 { hi: i as f64, lo: 0.0 }, f64x2::zero()).powc(s);
            ans += x;
        }
        println!("ans = {}", ans);
        let gt = Complex::new(
            f64x2 { hi: 1.1409175704943191e7, lo: 8.194110941294255e-10 },
            f64x2 { hi: 2.8472593984260815e6, lo: 1.4500054775776998e-10 },
        );
        assert_complex_close(ans, gt, 6e-32);
    }
}
