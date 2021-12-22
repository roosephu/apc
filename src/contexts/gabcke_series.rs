use F64x2::f64x2;

use super::{lazy_static::UninitCell, utils::*};
use crate::traits::MyReal;

pub trait GabckeSeries {
    fn gabcke_series(&self, eps: f64) -> Self;
}

pub trait GabckeSeriesCoeffs: Sized {
    fn gabcke_series_coeffs() -> &'static Vec<Vec<Self>>;
}

// TODO: determine K wisely
const fn calc_gabcke_n_terms(_t: f64, _eps: f64) -> usize { 7 }

impl<T: MyReal + GabckeSeriesCoeffs> GabckeSeries for T {
    /// See [Gabcke] Section 2.2. See also Theorem 2.1.6.
    /// For numerical stability, we use the Chebyshev's version (Eq 2.38)
    /// instead of polynomial version (Eq 2.34, 2.35).
    fn gabcke_series(&self, eps: f64) -> Self {
        let coeffs = T::gabcke_series_coeffs();

        let t = *self;
        let a = (t / Self::PI() / 2.0).sqrt();
        let n = a.floor();
        let z = Self::one() - (a - n) * 2.0;

        let mut expansion = Self::zero();
        let K = calc_gabcke_n_terms(a.fp(), eps);

        // We evaluate T_{2k}(z) = T_k(T_2(z)).
        let t2_z = z * z * 2.0 - 1.0;
        let mut inv_pow_a = Self::one();
        for k in 0..=K {
            let series = &coeffs[k];
            let m = series.len();
            let mut s = series[0];

            let mut pre = Self::one();
            let mut cur = t2_z;
            for j in 1..m {
                s += series[j] * cur;
                (pre, cur) = (cur, cur * t2_z * 2.0 - pre);
            }
            if k % 2 == 1 {
                s *= z;
            }
            expansion += s * inv_pow_a;
            inv_pow_a /= a;
            // println!("k = {}, correction = {:?}", k, s);
        }
        expansion / a.sqrt() * (if n.to_usize().unwrap() % 2 == 0 { -1.0 } else { 1.0 })
    }
}

#[allow(non_upper_case_globals)]
static TABLE_f64: UninitCell<Vec<Vec<f64>>> = UninitCell::uninit();
#[allow(non_upper_case_globals)]
static TABLE_f64x2: UninitCell<Vec<Vec<f64x2>>> = UninitCell::uninit();

const N_COLS: usize = 26;
const N_ROWS: usize = 11;

pub fn init() {
    let data = read_data("data/gabcke_series.txt", 1000)
        .expect("can't load Bernoulli numbers from `data/gabcke_series.txt`");
    assert!(data.len() == N_COLS * N_ROWS);

    fn group<T: Clone>(data: Vec<T>) -> Vec<Vec<T>> {
        (0..N_COLS).map(|x| data[x * N_ROWS..(x + 1) * N_ROWS].to_vec()).collect()
    }

    TABLE_f64.set(group(data.iter().map(|x| x.to_f64()).collect()));
    TABLE_f64x2.set(group(data.iter().map(mpf_to_f64x2).collect()));
}

impl GabckeSeriesCoeffs for f64 {
    fn gabcke_series_coeffs() -> &'static Vec<Vec<Self>> { &*TABLE_f64 }
}

impl GabckeSeriesCoeffs for f64x2 {
    fn gabcke_series_coeffs() -> &'static Vec<Vec<Self>> { &*TABLE_f64x2 }
}
