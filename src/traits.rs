use num::{Complex, Signed};
use num_traits::{AsPrimitive, Float, FloatConst, NumAssignOps, Pow};
use rustfft::FftNum;
use std::{
    fmt::{Debug, Display},
    ops::{Add, Div, Mul, Sub},
};

use crate::sum_trunc_dirichlet::ExpPolyApprox;

pub trait MyFloat = Float
    + FloatConst
    + Signed
    + NumAssignOps
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
    + ExpPolyApprox
    + FftNum
    + 'static;

pub trait ComplexFunctions {
    fn mul_i(self) -> Self;
    fn approx(self) -> Complex<f64>;
}

impl<T: MyFloat> ComplexFunctions for Complex<T> {
    #[inline]
    fn mul_i(self) -> Self { Self { re: -self.im, im: self.re } }

    #[inline]
    fn approx(self) -> Complex<f64> {
        Complex::<f64> {
            re: AsPrimitive::<f64>::as_(self.re),
            im: AsPrimitive::<f64>::as_(self.im),
        }
    }
}
