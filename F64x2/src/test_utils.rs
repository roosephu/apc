use std::fmt::Display;

use crate::traits::FpOps;
use num::{Complex, Float};
use num_traits::NumOps;

pub fn assert_close<T: FpOps + NumOps + Display + Copy>(a: T, b: T, eps: f64) {
    let rel_abs = (a - b).fp().abs() / b.fp().abs();
    assert!(rel_abs < eps, "rel abs = {:.E}, a = {}, b = {}, diff = {}", rel_abs, a, b, a - b);
}

pub fn assert_complex_close<T: FpOps + Float + Display + Copy>(
    a: Complex<T>,
    b: Complex<T>,
    eps: f64,
) {
    let rel_abs = (a - b).norm().fp() / b.norm().fp();
    assert!(rel_abs < eps, "rel abs = {:.E}, a = {}, b = {}, diff = {}", rel_abs, a, b, a - b);
}
