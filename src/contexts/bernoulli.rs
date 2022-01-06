use super::utils::{impl_from_uninit_cell, mpf_to_f64x2, read_data};
use F64x2::f64x2;

pub trait Bernoulli {
    fn bernoulli(n: usize) -> Self;
}

const BERNOULLI_F64: [f64; 101] = include!("../../const_table/bernoulli_f64.rs");
const BERNOULLI_F64X2: [f64x2; 101] = include!("../../const_table/bernoulli_f64x2.rs");

impl Bernoulli for f64 {
    fn bernoulli(n: usize) -> Self { BERNOULLI_F64[n] }
}

impl Bernoulli for f64x2 {
    fn bernoulli(n: usize) -> Self { BERNOULLI_F64X2[n] }
}
