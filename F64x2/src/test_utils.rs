use num::Complex;
use crate::f64x2;

pub fn assert_close(a: f64x2, b: f64x2, eps: f64) {
    let rel_abs = (a - b).abs().approx() / b.abs().approx();
    assert!(rel_abs < eps, "rel abs = {:.E}, a = {}, b = {}, diff = {}", rel_abs, a, b, a - b);
}

pub fn assert_complex_close(a: Complex<f64x2>, b: Complex<f64x2>, eps: f64) {
    let rel_abs = (a - b).norm().approx() / b.norm().approx();
    assert!(rel_abs < eps, "rel abs = {:.E}, a = {}, b = {}, diff = {}", rel_abs, a, b, a - b);
}
