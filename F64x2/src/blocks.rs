#![allow(clippy::suspicious_arithmetic_impl)]
use std::ops::{Add, Div, Mul, Sub};

use crate::f64x2;

/// Algorithm 1
#[inline]
pub fn two_add_fast(a: f64, b: f64) -> (f64, f64) {
    let hi = a + b;
    let z = hi - a;
    let lo = b - z;
    (hi, lo)
}

#[inline]
pub fn two_add_if(a: f64, b: f64) -> (f64, f64) {
    if a.abs() > b.abs() {
        two_add_fast(a, b)
    } else {
        two_add_fast(b, a)
    }
}

/// Algorithm 2
#[inline]
pub fn two_add(a: f64, b: f64) -> (f64, f64) {
    let hi = a + b;
    let a1 = hi - b;
    let b1 = hi - a1;
    let lo = (a - a1) + (b - b1);
    (hi, lo)
}

#[inline]
pub fn two_sub(a: f64, b: f64) -> (f64, f64) {
    let hi = a - b;
    let a1 = hi + b;
    let b1 = hi - a1;
    let lo = (a - a1) - (b + b1);
    (hi, lo)
}

/// Algorithm 3
#[inline]
pub(crate) fn two_mul(a: f64, b: f64) -> (f64, f64) {
    let s = a * b;
    let t = fma(a, b, -s);
    (s, t)
}

#[inline]
pub fn fma(a: f64, b: f64, c: f64) -> f64 { a.mul_add(b, c) }

impl Add<f64> for f64x2 {
    type Output = f64x2;

    /// Algorithm 4
    #[inline]
    fn add(self, y: f64) -> f64x2 {
        let x = self;
        let s = two_add(x.hi, y);
        let v = x.lo + s.1;
        f64x2::from(two_add_fast(s.0, v))
    }
}

impl Add<f64x2> for f64x2 {
    type Output = f64x2;

    /// Algorithm 6 in [5, 6]
    #[inline]
    fn add(self, y: f64x2) -> f64x2 {
        let x = self;
        let s = two_add(x.hi, y.hi);
        let t = two_add(x.lo, y.lo);
        let c = s.1 + t.0;
        let v = two_add_fast(s.0, c);
        let w = t.1 + v.1;
        Self::from(two_add_fast(v.0, w))
    }
}

impl Sub<f64> for f64x2 {
    type Output = f64x2;

    /// Algorithm 4
    #[inline]
    fn sub(self, y: f64) -> f64x2 {
        let x = self;
        let s = two_sub(x.hi, y);
        let v = x.lo + s.1;
        Self::from(two_add_fast(s.0, v))
    }
}

impl Sub<f64x2> for f64x2 {
    type Output = f64x2;

    /// Algorithm 6 in [5, 6]
    #[inline]
    fn sub(self, y: f64x2) -> f64x2 {
        let x = self;
        let s = two_sub(x.hi, y.hi);
        let t = two_sub(x.lo, y.lo);
        let c = s.1 + t.0;
        let v = two_add_fast(s.0, c);
        let w = t.1 + v.1;
        Self::from(two_add_fast(v.0, w))
    }
}

impl Mul<f64> for f64x2 {
    type Output = f64x2;

    /// Algorithm 9 in [7, 8, 9]
    #[inline]
    fn mul(self, y: f64) -> f64x2 {
        let x = self;
        let c = two_mul(x.hi, y);
        let c_lo2 = x.lo.mul_add(y, c.1);
        Self::from(two_add_fast(c.0, c_lo2))
    }
}

impl Mul<f64x2> for f64x2 {
    type Output = f64x2;

    /// Algorithm 12 in [10, 11, 12]
    #[inline]
    fn mul(self, y: f64x2) -> f64x2 {
        let x = self;
        let c = two_mul(x.hi, y.hi);
        let t_lo0 = x.lo * y.lo;
        let t_lo1 = x.hi.mul_add(y.lo, t_lo0);
        let c_lo2 = x.lo.mul_add(y.hi, t_lo1);
        let c_lo3 = c.1 + c_lo2;
        Self::from(two_add_fast(c.0, c_lo3))
    }
}

impl Div<f64> for f64x2 {
    type Output = f64x2;

    /// Algorithm 15 in [13, 14, 15]
    #[inline]
    fn div(self, y: f64) -> f64x2 {
        let x = self;
        let t_hi = x.hi / y;
        let p = two_mul(t_hi, y);
        let d_hi = x.hi - p.0;
        let d_lo = x.lo - p.1;
        let d = d_hi + d_lo;
        let t_lo = d / y;
        Self::from(two_add_fast(t_hi, t_lo))
    }
}

impl Div<f64x2> for f64x2 {
    type Output = f64x2;

    /// Algorithm 17 in [16, 17, 18]
    #[inline]
    fn div(self, y: f64x2) -> f64x2 {
        let x = self;
        let t_hi = x.hi / y.hi;
        let r = y * t_hi;
        let p_hi = x.hi - r.hi;
        let d_lo = x.lo - r.lo;
        let d = p_hi + d_lo;
        let t_lo = d / y.hi;
        Self::from(two_add_fast(t_hi, t_lo))
    }
}

impl Add<f64x2> for f64 {
    type Output = f64x2;

    #[inline]
    fn add(self, y: f64x2) -> f64x2 {
        let x = self;
        y + x
    }
}

impl Sub<f64x2> for f64 {
    type Output = f64x2;

    #[inline]
    fn sub(self, y: f64x2) -> f64x2 {
        let x = self;
        -y + x
    }
}

impl Mul<f64x2> for f64 {
    type Output = f64x2;

    #[inline]
    fn mul(self, y: f64x2) -> f64x2 {
        let x = self;
        y * x
    }
}

impl Div<f64x2> for f64 {
    type Output = f64x2;

    #[inline]
    fn div(self, y: f64x2) -> f64x2 {
        let x = self;
        f64x2 { hi: x, lo: 0.0 } / y
    }
}
