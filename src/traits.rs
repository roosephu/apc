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
    + AsPrimitive<i64>
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
    fn approx(self) -> Complex<f64>;
}

impl<T: GenericFloat> ComplexFunctions for Complex<T> {
    #[inline]
    fn mul_i(self) -> Self {
        Self {
            re: -self.im,
            im: self.re,
        }
    }

    #[inline]
    fn approx(self) -> Complex<f64> {
        Complex::<f64> {
            re: AsPrimitive::<f64>::as_(self.re),
            im: AsPrimitive::<f64>::as_(self.im),
        }
    }
}
