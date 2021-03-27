use num::{Complex, Signed};
use num_traits::{Float, FloatConst, NumAssignOps, Pow};
use rustfft::FftNum;
use std::{
    fmt::{Debug, Display, LowerExp},
    ops::{Add, Div, Mul, Sub},
};

use crate::{
    sum_trunc_dirichlet::ExpPolyApprox,
    unchecked_cast::{UncheckedCast, UncheckedFrom, UncheckedInto},
};

pub trait Erfc {
    fn erfc(self, eps: f64) -> Self;
}

pub trait MyFloat = Float
    + FloatConst
    + Signed
    + NumAssignOps
    + Display // might not need
    + Debug // might not need
    + LowerExp // might not need
    + Default
    + UncheckedInto<f64>
    + UncheckedFrom<f64>
    + UncheckedFrom<i64>
    + UncheckedInto<i64>
    + UncheckedFrom<i32>
    + UncheckedInto<i32>
    + UncheckedCast
    + Add<Complex<Self>, Output = Complex<Self>>
    + Sub<Complex<Self>, Output = Complex<Self>>
    + Mul<Complex<Self>, Output = Complex<Self>>
    + Div<Complex<Self>, Output = Complex<Self>>
    + Pow<i32, Output = Self>
    + ExpPolyApprox
    + FftNum
    + Erfc
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
        Complex::<f64> { re: self.re.unchecked_cast(), im: self.im.unchecked_cast() }
    }
}

impl Erfc for f64 {
    fn erfc(self, eps: f64) -> Self { rgsl::error::erfc(self) }
}
