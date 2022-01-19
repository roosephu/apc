use num::{Complex, Signed};
use num_traits::{Float, FloatConst, FromPrimitive, NumAssignOps, Pow};
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

pub trait ErgonomicOps = NumAssignOps
    + Display // might not need
    + Debug // might not need
    + LowerExp // might not need
    + Default;

pub trait MyReal = Float
    + FloatConst
    + Signed
    + FromPrimitive
    + FpOps
    + ComplexOps
    + ErgonomicOps
    + Pow<i32, Output = Self>
    + Sync
    + Send
    + 'static;

pub trait ZetaZerosDatabase<T> = Iterator<Item = (T, f64)>;
