use crate::f64xn::f64x2;
use crate::traits::MyReal;
use num::Complex;

pub fn assert_close(a: f64x2, b: f64x2, eps: f64) {
    let rel_abs = (a - b).abs() / b.abs();
    let rel_abs = rel_abs.unchecked_cast::<f64>();
    assert!(rel_abs < eps, "rel abs = {:.E}, a = {}, b = {}, diff = {}", rel_abs, a, b, a - b);
}

pub fn assert_complex_close(a: Complex<f64x2>, b: Complex<f64x2>, eps: f64) {
    let rel_abs = (a - b).norm() / b.norm();
    let rel_abs = rel_abs.unchecked_cast::<f64>();
    assert!(rel_abs < eps, "rel abs = {:.E}, a = {}, b = {}, diff = {}", rel_abs, a, b, a - b);
}
