use num::{Complex, Signed};
use num_traits::{Float, FloatConst, FromPrimitive, Num, NumAssignOps, NumOps, Pow};
use std::{
    fmt::{Debug, Display, LowerExp},
    ops::{Add, Div, Mul, Sub},
};
use F64x2::traits::FpOps;

pub trait ComplexOps = Float
    + Add<Complex<Self>, Output = Complex<Self>>
    + Sub<Complex<Self>, Output = Complex<Self>>
    + Mul<Complex<Self>, Output = Complex<Self>>
    + Div<Complex<Self>, Output = Complex<Self>>;

pub trait BaseReal = Copy
    + Clone
    + FpOps
    + NumAssignOps
    + Signed
    + Num
    + Sync
    + Send
    + 'static
    + Display
    + Debug
    + Default;

pub trait MyReal =
    BaseReal + Float + FloatConst + FromPrimitive + ComplexOps + Pow<i32, Output = Self>;

pub trait ZetaZerosDatabase<T> = Iterator<Item = (T, f64)>;
