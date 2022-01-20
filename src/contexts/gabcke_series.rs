use F64x2::f64x2;

use crate::traits::MyReal;

pub trait GabckeSeries {
    fn gabcke_series(&self, order: usize, eps: f64) -> Self;
}

pub trait GabckeSeriesCoeffs: Sized {
    fn gabcke_series_coeffs(k: usize) -> &'static [Self];
}

const GABCKE_SERIES_F64X2: &[[f64x2; 26]] = &include!("../../const_table/gabcke_series_f64x2.rs");
const GABCKE_SERIES_F64: &[[f64; 26]] = &include!("../../const_table/gabcke_series_f64.rs");

impl GabckeSeriesCoeffs for f64 {
    fn gabcke_series_coeffs(k: usize) -> &'static [f64] { &GABCKE_SERIES_F64[k] }
}

impl GabckeSeriesCoeffs for f64x2 {
    fn gabcke_series_coeffs(k: usize) -> &'static [f64x2] { &GABCKE_SERIES_F64X2[k] }
}

impl<T: MyReal + GabckeSeriesCoeffs> GabckeSeries for T {
    /// See [Gabcke] Section 2.2. See also Theorem 2.1.6.  For numerical
    /// stability, we use the Chebyshev's version (Eq 2.38) instead of
    /// polynomial version (Eq 2.34, 2.35).
    #[inline(never)]
    fn gabcke_series(&self, order: usize, atol: f64) -> Self {
        let t = *self;
        let a = (t / T::PI() / 2.0).sqrt();
        let n = a.floor().fp();
        let z = T::one() - (a - n) * 2.0;
        let eps = atol * 0.1;

        let mut expansion = T::zero();

        // We evaluate T_{2k}(z) = T_k(T_2(z)).
        let t2_z = z * z * 2.0 - 1.0;
        let mut inv_pow_a = Self::one();
        for k in 0..=order {
            let series = T::gabcke_series_coeffs(k);
            let m = series.len();
            let mut s = series[0];
            let eps = eps / inv_pow_a.fp().abs();

            let mut pre = T::one();
            let mut cur = t2_z;
            for j in 1..m {
                s += series[j] * cur;
                if series[j].fp().abs() < eps {
                    break;
                }
                (pre, cur) = (cur, cur * t2_z * 2.0 - pre);
            }
            if k % 2 == 1 {
                s *= z;
            }
            expansion += s * inv_pow_a;
            inv_pow_a /= a;
            // println!("k = {}, correction = {:?}", k, s);
        }
        expansion / a.sqrt() * (if n as usize % 2 == 0 { -1.0 } else { 1.0 })
    }
}
