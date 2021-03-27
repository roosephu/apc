#![allow(clippy::suspicious_arithmetic_impl)]
use std::ops::{Add, Div, Mul, Sub};

use crate::f64xn::f64x2;

/// Algorithm 1
#[inline]
pub(crate) fn two_add_fast(a: f64, b: f64) -> f64x2 {
    let hi = a + b;
    let z = hi - a;
    let lo = b - z;
    f64x2 { hi, lo }
}

/// Algorithm 2
#[inline]
pub(crate) fn two_add(a: f64, b: f64) -> f64x2 {
    let hi = a + b;
    let a1 = hi - b;
    let b1 = hi - a1;
    let lo = (a - a1) + (b - b1);
    f64x2 { hi, lo }
}

#[inline]
pub(crate) fn two_sub(a: f64, b: f64) -> f64x2 {
    let hi = a - b;
    let a1 = hi + b;
    let b1 = hi - a1;
    let lo = (a - a1) - (b + b1);
    f64x2 { hi, lo }
}

/// Algorithm 3
#[inline]
pub(crate) fn two_mul(a: f64, b: f64) -> f64x2 {
    let s = a * b;
    let t = fma(a, b, -s);
    f64x2 { hi: s, lo: t }
}

#[inline]
pub(crate) fn fma(a: f64, b: f64, c: f64) -> f64 { a.mul_add(b, c) }

impl f64x2 {
    #[inline]
    pub fn new(hi: f64, lo: f64) -> Self { Self { hi, lo } }
}

impl Add<f64> for f64x2 {
    type Output = f64x2;

    /// Algorithm 4
    #[inline]
    fn add(self, y: f64) -> f64x2 {
        let x = self;
        let s = two_add(x.hi, y);
        let v = x.lo + s.lo;
        two_add_fast(s.hi, v)
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
        let c = s.lo + t.hi;
        let v = two_add_fast(s.hi, c);
        let w = t.lo + v.lo;
        two_add_fast(v.hi, w)
    }
}

impl Sub<f64> for f64x2 {
    type Output = f64x2;

    /// Algorithm 4
    #[inline]
    fn sub(self, y: f64) -> f64x2 {
        let x = self;
        let s = two_sub(x.hi, y);
        let v = x.lo + s.lo;
        two_add_fast(s.hi, v)
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
        let c = s.lo + t.hi;
        let v = two_add_fast(s.hi, c);
        let w = t.lo + v.lo;
        two_add_fast(v.hi, w)
    }
}

impl Mul<f64> for f64x2 {
    type Output = f64x2;

    /// Algorithm 9 in [7, 8, 9]
    #[inline]
    fn mul(self, y: f64) -> f64x2 {
        let x = self;
        let c = two_mul(x.hi, y);
        let c_lo2 = x.lo.mul_add(y, c.lo);
        two_add_fast(c.hi, c_lo2)
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
        let c_lo3 = c.lo + c_lo2;
        two_add_fast(c.hi, c_lo3)
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
        let d_hi = x.hi - p.hi;
        let d_lo = x.lo - p.lo;
        let d = d_hi + d_lo;
        let t_lo = d / y;
        two_add_fast(t_hi, t_lo)
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
        two_add_fast(t_hi, t_lo)
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
