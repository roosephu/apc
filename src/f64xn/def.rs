use std::{
    fmt::{Display, Error},
    num::FpCategory,
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Rem, RemAssign, Sub, SubAssign},
};

use num::{Complex, Float, FromPrimitive, Num, One, Signed, ToPrimitive, Zero};
use num_traits::{AsPrimitive, FloatConst, Pow};

use crate::sum_trunc_dirichlet::ExpPolyApprox;

#[allow(non_camel_case_types)]
#[derive(Clone, Copy, PartialEq, PartialOrd, Debug, Default)]
pub struct f64x2 {
    pub hi: f64,
    pub lo: f64,
}

impl From<f64x2> for String {
    fn from(mut a: f64x2) -> Self {
        // TODO
        if a.is_nan() {
            String::from("NaN")
        } else if a.is_infinite() {
            if a.is_sign_negative() {
                String::from("-Inf")
            } else {
                String::from("Inf")
            }
        } else if a.is_zero() {
            String::from("0")
        } else {
            let mut ret = String::from("");
            if a.is_sign_negative() {
                a = -a;
                ret.push('-');
            }
            let mut e = 0;
            while a / 10.0 >= f64x2::one() {
                a = a / 10.0;  // TODO: avoid division
                e += 1;
            }
            while a * 10.0 < f64x2::one() {
                a = a * 10.0;
                e -= 1;
            }

            let mut dec_point = false;

            for _ in 0..30 {
                let d = a.floor().hi as u8;  // a is in [0, 8] so it's fine to use floor
                ret.push((b'0' + d) as char);
                a = (a - f64x2::from_f64(d as f64).unwrap()) * 10.;
                if a.is_zero() {
                    break;
                }

                if !dec_point {
                    ret.push('.');
                    dec_point = true;
                }
            }
            if e != 0 {
                ret.push('E');
                ret.push_str(e.to_string().as_str());
            }

            ret
        }
    }
}

impl Display for f64x2 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", String::from(*self))
    }
}

impl Zero for f64x2 {
    fn zero() -> Self {
        Self { hi: 0.0, lo: 0.0 }
    }

    fn is_zero(&self) -> bool {
        self.hi == 0.0 && self.lo == 0.0
    }

    fn set_zero(&mut self) {
        self.hi = 0.0;
        self.lo = 0.0;
    }
}

impl One for f64x2 {
    fn one() -> Self {
        Self { hi: 1.0, lo: 0.0 }
    }

    fn is_one(&self) -> bool {
        self.hi == 1.0 && self.lo == 0.0
    }

    fn set_one(&mut self) {
        self.hi = 1.0;
        self.lo = 0.0;
    }
}

impl Neg for f64x2 {
    type Output = f64x2;

    fn neg(self) -> Self::Output {
        f64x2 {
            hi: -self.hi,
            lo: -self.lo,
        }
    }
}

impl Rem for f64x2 {
    type Output = f64x2;

    fn rem(self, rhs: f64x2) -> Self::Output {
        unimplemented!()
    }
}

impl Num for f64x2 {
    type FromStrRadixErr = Error;
    fn from_str_radix(str: &str, radix: u32) -> Result<Self, Self::FromStrRadixErr> {
        todo!()
    }
}

impl ToPrimitive for f64x2 {
    fn to_f32(&self) -> Option<f32> {
        todo!()
    }
    fn to_f64(&self) -> Option<f64> {
        Some(self.as_())
    }
    fn to_i128(&self) -> Option<i128> {
        todo!()
    }
    fn to_i16(&self) -> Option<i16> {
        todo!()
    }
    fn to_i32(&self) -> Option<i32> {
        todo!()
    }
    fn to_i64(&self) -> Option<i64> {
        Some(self.as_())
    }
    fn to_i8(&self) -> Option<i8> {
        todo!()
    }
    fn to_isize(&self) -> Option<isize> {
        todo!()
    }
    fn to_u128(&self) -> Option<u128> {
        todo!()
    }
    fn to_u16(&self) -> Option<u16> {
        todo!()
    }
    fn to_u32(&self) -> Option<u32> {
        todo!()
    }
    fn to_u64(&self) -> Option<u64> {
        todo!()
    }
    fn to_u8(&self) -> Option<u8> {
        todo!()
    }
    fn to_usize(&self) -> Option<usize> {
        todo!()
    }
}

impl num::traits::NumCast for f64x2 {
    fn from<T: ToPrimitive>(n: T) -> Option<Self> {
        todo!()
    }
}

impl Float for f64x2 {
    fn classify(self) -> FpCategory {
        todo!()
    }
    fn epsilon() -> Self {
        todo!()
    }
    fn hypot(self, other: Self) -> Self {
        todo!()
    }
    fn nan() -> Self {
        todo!()
    }
    fn infinity() -> Self {
        todo!()
    }
    fn neg_infinity() -> Self {
        todo!()
    }
    fn neg_zero() -> Self {
        todo!()
    }
    fn min_value() -> Self {
        todo!()
    }
    fn min_positive_value() -> Self {
        todo!()
    }
    fn max_value() -> Self {
        todo!()
    }
    fn is_nan(self) -> bool {
        todo!()
    }
    fn is_finite(self) -> bool {
        todo!()
    }
    fn is_infinite(self) -> bool {
        todo!()
    }
    fn is_normal(self) -> bool {
        todo!()
    }
    fn floor(self) -> Self {
        todo!()
    }
    fn ceil(self) -> Self {
        todo!()
    }
    fn round(self) -> Self {
        todo!()
    }
    fn trunc(self) -> Self {
        todo!()
    }
    fn fract(self) -> Self {
        todo!()
    }
    fn abs(self) -> Self {
        todo!()
    }
    fn signum(self) -> Self {
        todo!()
    }
    fn is_sign_negative(self) -> bool {
        todo!()
    }
    fn is_sign_positive(self) -> bool {
        todo!()
    }
    fn mul_add(self, a: Self, b: Self) -> Self {
        todo!()
    }
    fn recip(self) -> Self {
        todo!()
    }
    fn powi(self, n: i32) -> Self {
        todo!()
    }
    fn powf(self, n: Self) -> Self {
        todo!()
    }
    fn sqrt(self) -> Self {
        todo!()
    }
    fn exp(self) -> Self {
        todo!()
    }
    fn exp2(self) -> Self {
        todo!()
    }
    fn ln(self) -> Self {
        todo!()
    }
    fn log10(self) -> Self {
        todo!()
    }
    fn ln_1p(self) -> Self {
        todo!()
    }
    fn log(self, base: Self) -> Self {
        todo!()
    }
    fn log2(self) -> Self {
        todo!()
    }
    fn max(self, other: Self) -> Self {
        todo!()
    }
    fn min(self, other: Self) -> Self {
        todo!()
    }
    fn exp_m1(self) -> Self {
        todo!()
    }
    fn abs_sub(self, other: Self) -> Self {
        todo!()
    }
    fn acos(self) -> Self {
        todo!()
    }
    fn acosh(self) -> Self {
        todo!()
    }
    fn asin(self) -> Self {
        todo!()
    }
    fn asinh(self) -> Self {
        todo!()
    }
    fn atan(self) -> Self {
        todo!()
    }
    fn atan2(self, other: Self) -> Self {
        todo!()
    }
    fn atanh(self) -> Self {
        todo!()
    }
    fn tan(self) -> Self {
        todo!()
    }
    fn tanh(self) -> Self {
        todo!()
    }
    fn to_degrees(self) -> Self {
        todo!()
    }
    fn to_radians(self) -> Self {
        todo!()
    }
    fn sin(self) -> Self {
        todo!()
    }
    fn sin_cos(self) -> (Self, Self) {
        todo!()
    }
    fn sinh(self) -> Self {
        todo!()
    }
    fn cbrt(self) -> Self {
        todo!()
    }
    fn cos(self) -> Self {
        todo!()
    }
    fn cosh(self) -> Self {
        todo!()
    }
    fn integer_decode(self) -> (u64, i16, i8) {
        todo!()
    }
}

impl Signed for f64x2 {
    fn abs(&self) -> Self {
        todo!()
    }
    fn abs_sub(&self, other: &Self) -> Self {
        todo!()
    }
    fn signum(&self) -> Self {
        todo!()
    }
    fn is_negative(&self) -> bool {
        todo!()
    }
    fn is_positive(&self) -> bool {
        todo!()
    }
}

impl Add<Complex<f64x2>> for f64x2 {
    type Output = Complex<f64x2>;

    #[inline]
    fn add(self, rhs: Complex<f64x2>) -> Complex<f64x2> {
        Complex::<f64x2>::new(rhs.re + self, rhs.im)
    }
}

impl Sub<Complex<f64x2>> for f64x2 {
    type Output = Complex<f64x2>;

    #[inline]
    fn sub(self, rhs: Complex<f64x2>) -> Self::Output {
        Complex::<f64x2>::new(self - rhs.re, rhs.im)
    }
}

impl Mul<Complex<f64x2>> for f64x2 {
    type Output = Complex<f64x2>;

    #[inline]
    fn mul(self, rhs: Complex<f64x2>) -> Complex<f64x2> {
        Complex::<f64x2>::new(self * rhs.re, self * rhs.im)
    }
}

impl Div<Complex<f64x2>> for f64x2 {
    type Output = Complex<f64x2>;

    #[inline]
    fn div(self, rhs: Complex<f64x2>) -> Self::Output {
        let norm_sqr = rhs.norm_sqr();
        Complex::<f64x2>::new(self * rhs.re / norm_sqr, -self * rhs.im / norm_sqr)
    }
}

impl AddAssign<f64x2> for f64x2 {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs
    }
}

impl SubAssign<f64x2> for f64x2 {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs
    }
}

impl MulAssign<f64x2> for f64x2 {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs
    }
}

impl DivAssign<f64x2> for f64x2 {
    fn div_assign(&mut self, rhs: Self) {
        *self = *self / rhs
    }
}

impl RemAssign<f64x2> for f64x2 {
    fn rem_assign(&mut self, rhs: Self) {
        *self = *self % rhs
    }
}

impl Pow<i32> for f64x2 {
    type Output = Self;

    fn pow(self, rhs: i32) -> Self::Output {
        todo!()
    }
}

impl FloatConst for f64x2 {
    fn PI() -> Self {
        todo!()
    }
    fn LOG2_E() -> Self {
        todo!()
    }
    fn LN_10() -> Self {
        todo!()
    }
    fn LN_2() -> Self {
        todo!()
    }
    fn LOG10_2() -> Self {
        todo!()
    }
    fn LOG10_E() -> Self {
        todo!()
    }
    fn LOG2_10() -> Self {
        todo!()
    }
    fn FRAC_1_PI() -> Self {
        todo!()
    }
    fn FRAC_1_SQRT_2() -> Self {
        todo!()
    }
    fn FRAC_2_PI() -> Self {
        todo!()
    }
    fn FRAC_2_SQRT_PI() -> Self {
        todo!()
    }
    fn FRAC_PI_2() -> Self {
        todo!()
    }
    fn FRAC_PI_3() -> Self {
        todo!()
    }
    fn FRAC_PI_4() -> Self {
        todo!()
    }
    fn FRAC_PI_6() -> Self {
        todo!()
    }
    fn FRAC_PI_8() -> Self {
        todo!()
    }
    fn SQRT_2() -> Self {
        todo!()
    }
    fn TAU() -> Self {
        todo!()
    }
    fn E() -> Self {
        todo!()
    }
}

impl AsPrimitive<f64> for f64x2 {
    fn as_(self) -> f64 {
        self.hi
    }
}

impl AsPrimitive<i64> for f64x2 {
    fn as_(self) -> i64 {
        self.hi as i64 + self.lo as i64
    }
}

impl FromPrimitive for f64x2 {
    fn from_i64(n: i64) -> Option<Self> {
        // TODO: might over flow
        Some(Self {  hi: n as f64, lo: 0.0 })
    }
    fn from_u64(n: u64) -> Option<Self> {
        // TODO: might over flow
        Some(Self {  hi: n as f64, lo: 0.0 })
    }
    fn from_f64(n: f64) -> Option<Self> {
        Some(Self {  hi: n, lo: 0.0 })
    }
}

const COEFFS: [Complex<f64x2>; 0] = [];

impl ExpPolyApprox for f64x2 {
    type Output = impl Iterator<Item = (usize, Complex<Self>)>;

    fn get_poly_approx() -> Self::Output {
        COEFFS.iter().enumerate().map(|(idx, &x)| (idx, x))
    }
}
