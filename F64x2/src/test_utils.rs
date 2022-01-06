use std::fmt::Display;

use crate::traits::FpOps;
use num::{Complex, Float};
use num_traits::NumOps;

/// See also `isapprox` in Julia.
pub fn is_approx<T: FpOps + NumOps + Display + Copy>(a: T, b: T, atol: f64, rtol: f64) -> bool {
    let abs_err = (a - b).fp().abs();
    let rel_err = abs_err / a.fp().max(b.fp());
    abs_err < atol || rel_err < rtol
}

pub fn assert_close<T: FpOps + NumOps + Display + Copy>(a: T, b: T, atol: f64, rtol: f64) {
    let abs_err = (a - b).fp().abs();
    let rel_err = abs_err / a.fp().max(b.fp());

    assert!(abs_err < atol || rel_err < rtol, "a = {a}, b = {b}, abs err = {abs_err:.3e} > {atol:.3e}, rel err = {rel_err:.3e} > {rtol:.3e}");
}

pub fn assert_complex_close<T: FpOps + Float + Display + Copy>(
    a: Complex<T>,
    b: Complex<T>,
    eps: f64,
) {
    let rel_abs = (a - b).norm().fp() / b.norm().fp();
    assert!(rel_abs < eps, "rel abs = {:.E}, a = {}, b = {}, diff = {}", rel_abs, a, b, a - b);
}
