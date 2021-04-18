use num::{Complex, Signed};
use num_traits::{Float, FloatConst, NumAssignOps, Pow};
use rustfft::FftNum;
use std::{
    fmt::{Debug, Display, LowerExp},
    ops::{Add, Div, Mul, Sub},
};

use crate::unchecked_cast::{UncheckedCast, UncheckedFrom, UncheckedInto};

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

const COEFFS: [Complex<f64>; 18] = [
    Complex { re: 1.0, im: 1.0812144107799266e-23 },
    Complex { re: -4.479610786345225e-21, im: 1.0 },
    Complex { re: -0.5, im: 3.027415987199093e-19 },
    Complex { re: -8.413828845781633e-18, im: -0.16666666666666669 },
    Complex { re: 0.04166666666666679, im: 1.1656663448809618e-16 },
    Complex { re: -1.0603511152022404e-15, im: 0.008333333333332324 },
    Complex { re: -0.0013888888888827402, im: 5.789264486508273e-15 },
    Complex { re: -2.491995923872859e-14, im: -0.00019841269843586228 },
    Complex { re: 2.4801587374768556e-5, im: 6.704175576866034e-14 },
    Complex { re: -1.594987515102099e-13, im: 2.7557317787217356e-6 },
    Complex { re: -2.755729303110001e-7, im: 2.3127460502103687e-13 },
    Complex { re: -3.2663668749921504e-13, im: -2.5052389834713885e-8 },
    Complex { re: 2.087985316554709e-9, im: 2.5867211760028217e-13 },
    Complex { re: -2.2167241850689593e-13, im: 1.6041263496425594e-10 },
    Complex { re: -1.1352710114429515e-11, im: 8.943908448871146e-14 },
    Complex { re: -4.542339711641447e-14, im: -7.962911435347713e-13 },
    Complex { re: 5.979573239083729e-14, im: 7.185782517642856e-15 },
    Complex { re: -1.970149077208406e-15, im: 1.9701490772084063e-15 },
];

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
    fn mul_pow_i(self, k: usize) -> Self;
    fn approx(self) -> Complex<f64>;
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
        Complex::<f64> { re: self.re.unchecked_cast(), im: self.im.unchecked_cast() }
    }
}

impl Erfc for f64 {
    fn erfc(self, eps: f64) -> Self { rgsl::error::erfc(self) }
}

impl ExpPolyApprox for f64 {
    type Output = impl Iterator<Item = (usize, Complex<Self>)>;

    fn get_poly_approx() -> Self::Output { COEFFS.iter().enumerate().map(|(idx, &x)| (idx, x)) }
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
