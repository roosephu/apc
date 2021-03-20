use num::{Complex, Signed};
use num_traits::{AsPrimitive, Float, FloatConst, FromPrimitive, NumAssignOps};
use std::{
    fmt::{Debug, Display},
    ops::{Add, Div, Mul, Sub},
};

pub trait GenericFloat = Float
    + FloatConst
    + Signed
    + NumAssignOps
    + FromPrimitive
    + AsPrimitive<f64>
    + Display
    + Debug
    + Default
    + Add<Complex<Self>, Output = Complex<Self>>
    + Sub<Complex<Self>, Output = Complex<Self>>
    + Mul<Complex<Self>, Output = Complex<Self>>
    + Div<Complex<Self>, Output = Complex<Self>>
    + 'static;

pub trait Rotate90 {
    fn rotate90(self) -> Self;
}

impl<T: GenericFloat> Rotate90 for Complex<T> {
    #[inline]
    fn rotate90(self) -> Self {
        Self {
            re: -self.im,
            im: self.re,
        }
    }
}
