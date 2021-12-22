#![allow(clippy::approx_constant)]

use std::{
    fmt::{Display, LowerExp},
    num::FpCategory,
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Rem, RemAssign, Sub, SubAssign},
};

use crate::f64x2;
use num::{Complex, Float, FromPrimitive, Num, One, Signed, ToPrimitive, Zero};
use num_traits::{AsPrimitive, FloatConst, Pow};

impl From<f64x2> for String {
    fn from(mut a: f64x2) -> Self {
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
                a = a / 10.0; // TODO: avoid division
                e += 1;
            }
            while a < f64x2::one() {
                a = a * 10.0;
                e -= 1;
            }

            let mut dec_point = false;

            for _ in 0..30 {
                let d = a.floor().hi as u8; // a is in [0, 8] so it's fine to use floor
                ret.push((b'0' + d) as char);
                a = (a - f64x2::from(d as i32)) * 10.;
                if a.is_zero() {
                    break;
                }

                if !dec_point {
                    ret.push('.');
                    dec_point = true;
                }
            }
            if e != 0 {
                ret.push('e');
                ret.push_str(e.to_string().as_str());
            }

            ret
        }
    }
}

impl Display for f64x2 {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", String::from(*self))
    }
}

impl Zero for f64x2 {
    #[inline]
    fn zero() -> Self { Self { hi: 0.0, lo: 0.0 } }

    #[inline]
    fn is_zero(&self) -> bool { self.hi == 0.0 && self.lo == 0.0 }
}

impl One for f64x2 {
    #[inline]
    fn one() -> Self { Self { hi: 1.0, lo: 0.0 } }

    #[inline]
    fn is_one(&self) -> bool { self.hi == 1.0 && self.lo == 0.0 }
}

impl Neg for f64x2 {
    type Output = f64x2;

    #[inline]
    fn neg(self) -> Self::Output { f64x2 { hi: -self.hi, lo: -self.lo } }
}

impl Rem for f64x2 {
    type Output = f64x2;

    #[inline]
    fn rem(self, rhs: f64x2) -> Self::Output { todo!() }
}

impl Num for f64x2 {
    type FromStrRadixErr = ();

    #[inline]
    fn from_str_radix(str: &str, radix: u32) -> Result<Self, Self::FromStrRadixErr> {
        let mut ret = f64x2::zero();
        let mut after_point = false;
        let mut pow = f64x2::one();
        for &c in str.as_bytes() {
            if c == b'.' {
                after_point = true
            } else {
                if after_point {
                    pow = pow / 10.0;
                }
                let d = (c - b'0') as f64;
                if after_point {
                    ret += pow * d;
                } else {
                    ret = ret * 10.0 + d;
                }
            }
        }
        Ok(ret)
    }
}

impl ToPrimitive for f64x2 {
    #[inline]
    fn to_f32(&self) -> Option<f32> { todo!() }

    #[inline]
    fn to_f64(&self) -> Option<f64> { Some(self.as_()) }

    #[inline]
    fn to_i128(&self) -> Option<i128> { todo!() }

    #[inline]
    fn to_i16(&self) -> Option<i16> { todo!() }

    #[inline]
    fn to_i32(&self) -> Option<i32> { Some(self.hi as i32 + self.lo as i32) }

    #[inline]
    fn to_i64(&self) -> Option<i64> { Some(self.as_()) }

    #[inline]
    fn to_i8(&self) -> Option<i8> { todo!() }

    #[inline]
    fn to_isize(&self) -> Option<isize> { todo!() }

    #[inline]
    fn to_u128(&self) -> Option<u128> { todo!() }

    #[inline]
    fn to_u16(&self) -> Option<u16> { todo!() }

    #[inline]
    fn to_u32(&self) -> Option<u32> { todo!() }

    #[inline]
    fn to_u64(&self) -> Option<u64> { Some((self.hi as i64 + self.lo as i64) as u64) }

    #[inline]
    fn to_u8(&self) -> Option<u8> { todo!() }

    #[inline]
    fn to_usize(&self) -> Option<usize> { Some((self.hi as i64 + self.lo as i64) as usize) }
}

impl num::traits::NumCast for f64x2 {
    #[inline]
    fn from<T: ToPrimitive>(n: T) -> Option<Self> { todo!() }
}

impl Float for f64x2 {
    #[inline]
    fn classify(self) -> FpCategory { todo!() }

    #[inline]
    fn epsilon() -> Self { Self::mp(1.232595164407831e-32) }

    #[inline]
    fn hypot(self, other: Self) -> Self { f64x2::hypot(self, other) }

    #[inline]
    fn nan() -> Self { todo!() }

    #[inline]
    fn infinity() -> Self { todo!() }

    #[inline]
    fn neg_infinity() -> Self { todo!() }

    #[inline]
    fn neg_zero() -> Self { todo!() }

    #[inline]
    fn min_value() -> Self { todo!() }

    #[inline]
    fn min_positive_value() -> Self { todo!() }

    #[inline]
    fn max_value() -> Self { todo!() }

    #[inline]
    fn is_nan(self) -> bool { self.hi.is_nan() }

    #[inline]
    fn is_finite(self) -> bool { self.hi.is_finite() }

    #[inline]
    fn is_infinite(self) -> bool { self.hi.is_infinite() }

    #[inline]
    fn is_normal(self) -> bool { todo!() }

    #[inline]
    fn floor(self) -> Self { self.floor() }

    #[inline]
    fn ceil(self) -> Self { self.ceil() }

    #[inline]
    fn round(self) -> Self { self.round() }

    #[inline]
    fn trunc(self) -> Self { todo!() }

    #[inline]
    fn fract(self) -> Self { todo!() }

    #[inline]
    fn abs(self) -> Self { self.abs() }

    #[inline]
    fn signum(self) -> Self { todo!() }

    #[inline]
    fn is_sign_negative(self) -> bool { self.hi < 0.0 }

    #[inline]
    fn is_sign_positive(self) -> bool { self.hi > 0.0 }

    #[inline]
    fn mul_add(self, a: Self, b: Self) -> Self { todo!() }

    #[inline]
    fn recip(self) -> Self { f64x2::recip(self) }

    #[inline]
    fn powi(self, n: i32) -> Self { self.powi(n as i64) }

    #[inline]
    fn powf(self, n: Self) -> Self { self.powf(n) }

    #[inline]
    fn sqrt(self) -> Self { self.sqrt() }

    #[inline]
    fn exp(self) -> Self { self.exp() }

    #[inline]
    fn exp2(self) -> Self { todo!() }

    #[inline]
    fn ln(self) -> Self { self.ln() }

    #[inline]
    fn log10(self) -> Self { todo!() }

    #[inline]
    fn ln_1p(self) -> Self { todo!() }

    #[inline]
    fn log(self, base: Self) -> Self { todo!() }

    #[inline]
    fn log2(self) -> Self { todo!() }

    #[inline]
    fn max(self, other: Self) -> Self { todo!() }

    #[inline]
    fn min(self, other: Self) -> Self {
        if self.hi > other.hi {
            other
        } else if self.hi < other.hi {
            self
        } else if self.lo > other.lo {
            other
        } else {
            self
        }
    }

    #[inline]
    fn exp_m1(self) -> Self { todo!() }

    #[inline]
    fn abs_sub(self, other: Self) -> Self { todo!() }

    #[inline]
    fn acos(self) -> Self { todo!() }

    #[inline]
    fn acosh(self) -> Self { todo!() }

    #[inline]
    fn asin(self) -> Self { todo!() }

    #[inline]
    fn asinh(self) -> Self { todo!() }

    #[inline]
    fn atan(self) -> Self { self.atan() }

    #[inline]
    fn atan2(self, other: Self) -> Self { self.atan2(other) }

    #[inline]
    fn atanh(self) -> Self { todo!() }

    #[inline]
    fn tan(self) -> Self { todo!() }

    #[inline]
    fn tanh(self) -> Self { todo!() }

    #[inline]
    fn to_degrees(self) -> Self { todo!() }

    #[inline]
    fn to_radians(self) -> Self { todo!() }

    #[inline]
    fn sin(self) -> Self { self.sin() }

    #[inline]
    fn sin_cos(self) -> (Self, Self) { self.sin_cos() }

    #[inline]
    fn sinh(self) -> Self { self.sinh() }

    #[inline]
    fn cbrt(self) -> Self { todo!() }

    #[inline]
    fn cos(self) -> Self { self.cos() }

    #[inline]
    fn cosh(self) -> Self { self.cosh() }

    #[inline]
    fn integer_decode(self) -> (u64, i16, i8) { todo!() }
}

impl Signed for f64x2 {
    #[inline]
    fn abs(&self) -> Self {
        if self.hi < 0.0 {
            -*self
        } else {
            *self
        }
    }

    #[inline]
    fn abs_sub(&self, other: &Self) -> Self { todo!() }

    #[inline]
    fn signum(&self) -> Self { todo!() }

    #[inline]
    fn is_negative(&self) -> bool { self.hi < 0.0 }

    #[inline]
    fn is_positive(&self) -> bool { self.hi > 0.0 }
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
        Complex::<f64x2>::new(self - rhs.re, -rhs.im)
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
    #[inline]
    fn add_assign(&mut self, rhs: Self) { *self = *self + rhs }
}

impl SubAssign<f64x2> for f64x2 {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) { *self = *self - rhs }
}

impl MulAssign<f64x2> for f64x2 {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) { *self = *self * rhs }
}

impl DivAssign<f64x2> for f64x2 {
    #[inline]
    fn div_assign(&mut self, rhs: Self) { *self = *self / rhs }
}

impl RemAssign<f64x2> for f64x2 {
    #[inline]
    fn rem_assign(&mut self, rhs: Self) { *self = *self % rhs }
}

impl Pow<i32> for f64x2 {
    type Output = Self;

    #[inline]
    fn pow(self, rhs: i32) -> Self::Output { self.powi(rhs as i64) }
}

impl FloatConst for f64x2 {
    #[inline]
    fn PI() -> Self { Self { hi: 3.141592653589793, lo: 1.2246467991473532e-16 } }

    #[inline]
    fn E() -> Self { Self { hi: 2.718281828459045, lo: 1.4456468917292502e-16 } }

    #[inline]
    fn LOG2_E() -> Self { Self { hi: 1.4426950408889634, lo: 2.0355273740931033e-17 } }

    #[inline]
    fn LN_10() -> Self { Self { hi: 2.302585092994046, lo: -2.1707562233822494e-16 } }

    #[inline]
    fn LN_2() -> Self { Self { hi: 0.6931471805599453, lo: 2.3190468138462996e-17 } }

    #[inline]
    fn LOG10_2() -> Self { Self { hi: 0.3010299956639812, lo: -2.8037281277851704e-18 } }

    #[inline]
    fn LOG10_E() -> Self { Self { hi: 0.4342944819032518, lo: 1.098319650216765e-17 } }

    #[inline]
    fn LOG2_10() -> Self { Self { hi: 3.321928094887362, lo: 1.661617516973592e-16 } }

    #[inline]
    fn FRAC_1_PI() -> Self { Self { hi: 0.3183098861837907, lo: -1.9678676675182486e-17 } }

    #[inline]
    fn FRAC_1_SQRT_2() -> Self { Self { hi: 0.7071067811865476, lo: -4.833646656726457e-17 } }

    #[inline]
    fn FRAC_2_PI() -> Self { Self { hi: 0.6366197723675814, lo: -3.935735335036497e-17 } }

    #[inline]
    fn FRAC_2_SQRT_PI() -> Self { Self { hi: 1.1283791670955126, lo: 1.533545961316588e-17 } }

    #[inline]
    fn FRAC_PI_2() -> Self { Self { hi: 1.5707963267948966, lo: 6.123233995736766e-17 } }

    #[inline]
    fn FRAC_PI_3() -> Self { Self { hi: 1.0471975511965979, lo: -1.072081766451091e-16 } }

    #[inline]
    fn FRAC_PI_4() -> Self { Self { hi: 0.7853981633974483, lo: 3.061616997868383e-17 } }

    #[inline]
    fn FRAC_PI_6() -> Self { Self { hi: 0.5235987755982989, lo: -5.360408832255455e-17 } }

    #[inline]
    fn FRAC_PI_8() -> Self { Self { hi: 0.39269908169872414, lo: 1.5308084989341915e-17 } }

    #[inline]
    fn SQRT_2() -> Self { Self { hi: 1.4142135623730951, lo: -9.667293313452913e-17 } }

    #[inline]
    fn TAU() -> Self { Self { hi: 6.283185307179586, lo: 2.4492935982947064e-16 } }
}

impl AsPrimitive<f64> for f64x2 {
    #[inline]
    fn as_(self) -> f64 { self.hi }
}

impl AsPrimitive<i64> for f64x2 {
    #[inline]
    fn as_(self) -> i64 { self.hi as i64 + self.lo as i64 }
}

impl FromPrimitive for f64x2 {
    #[inline]
    fn from_i64(n: i64) -> Option<Self> {
        let hi = n as f64;
        let lo = (n - hi as i64) as f64;
        Some(Self { hi, lo })
    }

    #[inline]
    fn from_u64(n: u64) -> Option<Self> {
        let hi = n as f64;
        let lo = n.wrapping_sub(hi as u64) as i64 as f64;
        Some(Self { hi, lo })
    }

    #[inline]
    fn from_f64(n: f64) -> Option<Self> { Some(Self { hi: n, lo: 0.0 }) }

    #[inline]
    fn from_u128(n: u128) -> Option<Self> {
        let hi = n as f64;
        let lo = n.wrapping_sub(hi as u128) as i128 as f64;
        Some(Self { hi, lo })
    }
}

impl From<i32> for f64x2 {
    #[inline]
    fn from(x: i32) -> Self { Self { hi: x as f64, lo: 0.0 } }
}

impl From<u64> for f64x2 {
    #[inline]
    fn from(n: u64) -> Self {
        let hi = n as f64;
        let lo = n.wrapping_sub(hi as u64) as i64 as f64;
        Self { hi, lo }
    }
}

// impl From<f64> for f64x2 {
//     #[inline]
//     fn from(x: f64) -> Self { f64x2 { hi: x, lo: 0.0 } }
// }

impl From<(f64, f64)> for f64x2 {
    #[inline]
    fn from((hi, lo): (f64, f64)) -> Self { f64x2 { hi, lo } }
}

impl Into<i64> for f64x2 {
    #[inline]
    fn into(self) -> i64 { self.hi as i64 + self.lo as i64 }
}

impl LowerExp for f64x2 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", String::from(*self))
    }
}
