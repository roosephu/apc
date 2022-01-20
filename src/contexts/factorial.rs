// use super::utils::{impl_from_uninit_cell, mpf_to_f64x2, read_data};
use F64x2::f64x2;

pub trait Factorial {
    fn factorial(n: usize) -> Self;
}

const FACTORIAL_F64: &[f64] = &include!("../../const_table/factorial_f64.rs");

impl Factorial for f64 {
    fn factorial(n: usize) -> Self { FACTORIAL_F64[n] }
}

const FACTORIAL_F64X2: &[f64x2] = &include!("../../const_table/factorial_f64x2.rs");
impl Factorial for f64x2 {
    fn factorial(n: usize) -> Self { FACTORIAL_F64X2[n] }
}
