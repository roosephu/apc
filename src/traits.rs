use num::{Complex, Signed};
use num_traits::{Float, FloatConst, FromPrimitive, NumAssignOps, Pow};
use std::{
    fmt::{Debug, Display, LowerExp},
    ops::{Add, Div, Mul, Sub},
};
use F64x2::f64x2;

pub trait Approximate {
    fn approx(&self) -> f64;
}

impl Approximate for f64 {
    #[inline]
    fn approx(&self) -> f64 { *self }
}

impl Approximate for f64x2 {
    #[inline]
    fn approx(&self) -> f64 { self.hi }
}

pub trait MyReal = Float
    + FloatConst
    + FromPrimitive
    + Signed
    + NumAssignOps
    + Display // might not need
    + Debug // might not need
    + LowerExp // might not need
    + Default
    + Add<f64, Output = Self>
    + Sub<f64, Output = Self>
    + Mul<f64, Output = Self>
    + Div<f64, Output = Self>
    + Add<Complex<Self>, Output = Complex<Self>>
    + Sub<Complex<Self>, Output = Complex<Self>>
    + Mul<Complex<Self>, Output = Complex<Self>>
    + Div<Complex<Self>, Output = Complex<Self>>
    + Pow<i32, Output = Self>
    + Sync
    + Send
    + Approximate
    + 'static;
