#![allow(clippy::approx_constant)]

use std::{
    fmt::{Display, LowerExp},
    num::FpCategory,
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Rem, RemAssign, Sub, SubAssign},
};

use num::{Complex, Float, FromPrimitive, Num, One, Signed, ToPrimitive, Zero};
use num_traits::{AsPrimitive, FloatConst, Pow};

use crate::{
    constants::RS_GABCKE_GAMMA,
    f64x2,
    traits::{Erfc, ExpPolyApprox, GabckeExpansion, Sinc},
    unchecked_cast::{UncheckedCast, UncheckedFrom, UncheckedInto},
};

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
                a = (a - (d as i32).unchecked_cast::<f64x2>()) * 10.;
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

    #[inline]
    fn set_zero(&mut self) {
        self.hi = 0.0;
        self.lo = 0.0;
    }
}

impl One for f64x2 {
    #[inline]
    fn one() -> Self { Self { hi: 1.0, lo: 0.0 } }

    #[inline]
    fn is_one(&self) -> bool { self.hi == 1.0 && self.lo == 0.0 }

    #[inline]
    fn set_one(&mut self) {
        self.hi = 1.0;
        self.lo = 0.0;
    }
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
    fn to_i32(&self) -> Option<i32> { todo!() }

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
    fn to_u64(&self) -> Option<u64> { todo!() }

    #[inline]
    fn to_u8(&self) -> Option<u8> { todo!() }

    #[inline]
    fn to_usize(&self) -> Option<usize> { todo!() }
}

impl num::traits::NumCast for f64x2 {
    #[inline]
    fn from<T: ToPrimitive>(n: T) -> Option<Self> { todo!() }
}

impl Float for f64x2 {
    #[inline]
    fn classify(self) -> FpCategory { todo!() }

    #[inline]
    fn epsilon() -> Self { Self::from(1.232595164407831e-32) }

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
    fn sin_cos(self) -> (Self, Self) { (self.sin(), self.cos()) }

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

const COEFFS: [Complex<f64x2>; 28] = [
    Complex::<f64x2> {
        re: f64x2 { hi: 1.0, lo: -7.41592015742533e-33 },
        im: f64x2 { hi: 0.0, lo: 0.0 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: -4.1405847908859944e-150, lo: -2.7044188792406287e-166 },
        im: f64x2 { hi: 1.0, lo: -7.278167722215507e-33 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: -0.5, lo: 1.1784183606192051e-30 },
        im: f64x2 { hi: -2.8844479592092057e-150, lo: 1.1454158016738223e-166 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: 1.380525771201648e-148, lo: 8.507802891506658e-165 },
        im: f64x2 { hi: -0.16666666666666666, lo: -9.251858538542578e-18 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: 0.041666666666666664, lo: 2.312964634604693e-18 },
        im: f64x2 { hi: 1.1906452414789667e-148, lo: -3.1085875021191995e-166 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: -1.664101265258138e-147, lo: -1.2349498604220513e-163 },
        im: f64x2 { hi: 0.008333333333333333, lo: 1.1564823172544488e-19 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: -0.001388888888888889, lo: 5.300543986595874e-20 },
        im: f64x2 { hi: -1.427683349621088e-147, lo: 2.445946111198975e-164 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: 9.937823711263713e-147, lo: -1.4453039209474945e-163 },
        im: f64x2 { hi: -0.0001984126984126984, lo: -1.7209553527463326e-22 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: 2.48015873015873e-5, lo: 2.151020312414919e-23 },
        im: f64x2 { hi: 7.523947356746099e-147, lo: -3.358929767724308e-163 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: -3.432799909335047e-146, lo: -2.0424307324901113e-162 },
        im: f64x2 { hi: 2.7557319223985893e-6, lo: -1.858395301572312e-22 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: -2.755731922398589e-7, lo: -2.3762056217062672e-23 },
        im: f64x2 { hi: -2.1723874828075946e-146, lo: -1.3627541464984915e-162 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: 7.515140858682638e-146, lo: -1.0757049173505409e-163 },
        im: f64x2 { hi: -2.505210838544172e-8, lo: 1.4493557661413441e-24 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: 2.08767569878681e-9, lo: -1.3262132406487786e-25 },
        im: f64x2 { hi: 3.8530912735021028e-146, lo: -1.9511596812789716e-162 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: -1.101480375319775e-145, lo: -2.0805076937678002e-162 },
        im: f64x2 { hi: 1.6059043836821613e-10, lo: 1.1618661939412562e-26 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: -1.1470745597729708e-11, lo: 5.8438510559461995e-28 },
        im: f64x2 { hi: -4.468201307595334e-146, lo: -8.975208601010553e-163 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: 1.1153698317092626e-145, lo: 7.829517776374605e-162 },
        im: f64x2 { hi: -7.647163731819804e-13, lo: -2.3009856010063025e-29 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: 4.779477332385702e-14, lo: -3.0710780807421415e-30 },
        im: f64x2 { hi: 3.4960514869622984e-146, lo: -2.1553825540367132e-162 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: -7.913231597925261e-146, lo: -1.676809158510425e-162 },
        im: f64x2 { hi: 2.811457254344474e-15, lo: 1.6032115815813728e-32 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: -1.561920696740829e-16, lo: 8.75341825032837e-33 },
        im: f64x2 { hi: -1.8616951964141138e-146, lo: -5.81157555159622e-163 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: 3.9264631351239826e-146, lo: 2.159639405941261e-162 },
        im: f64x2 { hi: -8.22063524597162e-18, lo: -6.246424667017666e-34 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: 4.1103175654770674e-19, lo: -6.632635734786154e-36 },
        im: f64x2 { hi: 6.656055583590879e-147, lo: -1.8468467153724547e-164 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: -1.3353889097920745e-146, lo: -6.684296744870731e-163 },
        im: f64x2 { hi: 1.957294077520665e-20, lo: -1.1296847467718588e-37 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: -8.89677188795462e-22, lo: -1.2759138958027608e-38 },
        im: f64x2 { hi: -1.5293832734904152e-147, lo: 1.155341504112017e-163 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: 2.9671574961042376e-147, lo: 2.7323001122765724e-163 },
        im: f64x2 { hi: -3.8681613701503016e-23, lo: 1.0579945255033147e-39 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: 1.6113071539190217e-24, lo: -2.952627456105686e-41 },
        im: f64x2 { hi: 2.040650482520839e-148, lo: -1.6781413118635922e-164 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: -3.8808736845917945e-148, lo: -2.1606761273305543e-164 },
        im: f64x2 { hi: 6.445182048382404e-26, lo: -1.7354019709548277e-43 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: -2.4235676574866644e-27, lo: 2.6354374780345496e-44 },
        im: f64x2 { hi: -1.201611139025832e-149, lo: -8.561291650615996e-166 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: 2.266317284126442e-149, lo: -1.6661220503998478e-165 },
        im: f64x2 { hi: -8.973190255604533e-29, lo: -2.0194856621200598e-45 },
    },
];

impl ExpPolyApprox for f64x2 {
    type Output = impl Iterator<Item = (usize, Complex<Self>)>;

    #[inline]
    fn get_poly_approx() -> Self::Output { COEFFS.iter().enumerate().map(|(idx, &x)| (idx, x)) }
}

impl From<i32> for f64x2 {
    fn from(x: i32) -> Self { Self { hi: x as f64, lo: 0.0 } }
}

impl From<f64> for f64x2 {
    fn from(x: f64) -> Self { f64x2 { hi: x, lo: 0.0 } }
}

impl LowerExp for f64x2 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", String::from(*self))
    }
}

impl Erfc for f64x2 {
    fn erfc(self, eps: f64) -> Self { self.erfc_eps(eps) }
}

impl UncheckedFrom<f64> for f64x2 {
    fn unchecked_from(x: f64) -> Self { Self { hi: x, lo: 0.0 } }
}

impl UncheckedFrom<i64> for f64x2 {
    fn unchecked_from(x: i64) -> Self {
        let hi = x as f64;
        let lo = (x - x as i64) as f64;
        Self { hi, lo }
    }
}

impl UncheckedFrom<i32> for f64x2 {
    fn unchecked_from(x: i32) -> Self { Self { hi: x as f64, lo: 0.0 } }
}

impl UncheckedInto<f64> for f64x2 {
    fn unchecked_into(self) -> f64 { self.hi }
}

impl UncheckedInto<i64> for f64x2 {
    fn unchecked_into(self) -> i64 { self.hi as i64 + self.lo as i64 }
}

impl UncheckedInto<i32> for f64x2 {
    fn unchecked_into(self) -> i32 { self.hi as i32 }
}

impl UncheckedCast for f64x2 {}

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
