use crate::f64x2;
use std::ops::{Add, Div, Mul, Sub};

pub trait FpOps:
    From<f64>
    + Add<f64, Output = Self>
    + Sub<f64, Output = Self>
    + Mul<f64, Output = Self>
    + Div<f64, Output = Self>
{
    fn fp(&self) -> f64;
    fn mp(x: f64) -> Self;
}

impl FpOps for f64 {
    #[inline]
    fn fp(&self) -> f64 { *self }

    #[inline]
    fn mp(x: f64) -> Self { x }
}

impl FpOps for f64x2 {
    #[inline]
    fn fp(&self) -> f64 { Self::approx(self) }

    #[inline]
    fn mp(x: f64) -> Self { Self::new(x, 0.0) }
}
