use crate::constants::*;
use num::{Complex, Signed};
use num_traits::{Float, FloatConst, NumAssignOps, Pow};
use rustfft::FftNum;
use std::{
    fmt::{Debug, Display, LowerExp},
    ops::{Add, Div, Mul, Sub},
};
use F64x2::f64x2;

pub trait Erfc {
    fn erfc(self, eps: f64) -> Self;
}

pub trait Sinc {
    fn sinc(self) -> Self;
}

pub trait ExpPolyApprox: Sized {
    type Output: Iterator<Item = (usize, Complex<Self>)>;

    fn get_poly_approx() -> Self::Output;
}

pub trait GabckeExpansion {
    fn expand(a: Self, z: Self, k: usize, eps: f64) -> Self;
}

pub trait HighPrecMod2PI {
    fn mod_2pi(self) -> Self;
}

pub trait MyReal = Float
    + FloatConst
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
    // + UncheckedInto<f64>
    // + UncheckedFrom<f64>
    // + UncheckedFrom<i64>
    // + UncheckedInto<i64>
    // + UncheckedFrom<i32>
    // + UncheckedInto<i32>
    // + UncheckedCast
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
    fn mul_pow_i(self, k: usize) -> Self;
    fn approx(self) -> Complex<f64>;
    fn exp_simul(self) -> Self;
}

impl<T: MyReal> ComplexFunctions for Complex<T> {
    #[inline]
    fn mul_i(self) -> Self { Self { re: -self.im, im: self.re } }

    fn mul_pow_i(self, k: usize) -> Self {
        match k % 4 {
            0 => self,
            1 => Self { re: -self.im, im: self.re },
            2 => Self { re: -self.re, im: -self.im },
            3 => Self { re: self.im, im: -self.re },
            _ => unreachable!(),
        }
    }

    #[inline]
    fn approx(self) -> Complex<f64> {
        Complex::<f64> { re: self.re.to_f64().unwrap(), im: self.im.to_f64().unwrap() }
    }

    #[inline(never)]
    fn exp_simul(self) -> Self {
        let r = self.re.exp();
        let (sin_im, cos_im) = self.im.sin_cos();
        Self { re: r * cos_im, im: r * sin_im }
    }
}

impl Erfc for f64 {
    fn erfc(self, eps: f64) -> Self { rgsl::error::erfc(self) }
}

impl ExpPolyApprox for f64 {
    type Output = impl Iterator<Item = (usize, Complex<Self>)>;

    fn get_poly_approx() -> Self::Output {
        EXP_POLY_EXP_F64.iter().enumerate().map(|(idx, &x)| (idx, x))
    }
}

impl Sinc for f64 {
    fn sinc(self) -> Self {
        if self == 0.0 {
            1.0
        } else {
            self.sin() / self
        }
    }
}

impl ExpPolyApprox for f64x2 {
    type Output = impl Iterator<Item = (usize, Complex<Self>)>;

    #[inline]
    fn get_poly_approx() -> Self::Output {
        F64X2_EXP_POLY_COEFFS.iter().enumerate().map(|(idx, &x)| (idx, x))
    }
}

impl Erfc for f64x2 {
    fn erfc(self, eps: f64) -> Self { self.erfc_eps(eps) }
}

impl Sinc for f64x2 {
    fn sinc(self) -> Self {
        if self.is_zero() {
            Self::one()
        } else {
            self.sin() / self
        }
    }
}

impl GabckeExpansion for f64x2 {
    fn expand(a: Self, z: Self, K: usize, eps: f64) -> Self {
        let mut expansion = Self::zero();
        let t2_z = z * z * 2.0 - 1.0;
        for k in 0..=K {
            let m = RS_GABCKE_GAMMA[k].len();
            let mut s = RS_GABCKE_GAMMA[k][0];

            let mut pre = Self::one();
            let mut cur = t2_z;
            for j in 1..m {
                s += RS_GABCKE_GAMMA[k][j] * cur;
                (pre, cur) = (cur, cur * t2_z * 2.0 - pre);
            }
            if k % 2 == 1 {
                s *= z;
            }
            expansion += s / a.powi(k as i64);
            println!("k = {}, correction = {:?}", k, s);
        }
        expansion
    }
}
