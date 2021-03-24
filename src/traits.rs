use num::{Complex, Signed};
use num_traits::{AsPrimitive, Float, FloatConst, FromPrimitive, NumAssignOps, Pow};
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
    + Pow<i32, Output = Self>
    + 'static;

pub trait ComplexFunctions {
    fn mul_i(self) -> Self;
}

impl<T: GenericFloat> ComplexFunctions for Complex<T> {
    #[inline]
    fn mul_i(self) -> Self {
        Self {
            re: -self.im,
            im: self.re,
        }
    }
}
