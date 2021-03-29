#![allow(clippy::approx_constant)]

use std::{
    fmt::{Display, LowerExp},
    num::FpCategory,
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Rem, RemAssign, Sub, SubAssign},
};

use num::{Complex, Float, FromPrimitive, Num, One, Signed, ToPrimitive, Zero};
use num_traits::{AsPrimitive, FloatConst, Pow};

use crate::{
    f64x2,
    sum_trunc_dirichlet::ExpPolyApprox,
    traits::Erfc,
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
    type FromStrRadixErr = &'static str;

    #[inline]
    fn from_str_radix(str: &str, radix: u32) -> Result<Self, Self::FromStrRadixErr> { todo!() }
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
    fn epsilon() -> Self { Self::from(f64::epsilon()) }

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
    fn floor(self) -> Self { todo!() }

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
    fn powi(self, n: i32) -> Self { todo!() }

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
    fn min(self, other: Self) -> Self { todo!() }

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
    fn PI() -> Self { Self { hi: 3.141592653589793, lo: 1.2246467991473515e-16 } }

    #[inline]
    fn E() -> Self { Self { hi: 2.718281828459045, lo: 1.4456468917292359e-16 } }

    #[inline]
    fn LOG2_E() -> Self { Self { hi: 1.4426950408889634, lo: 2.0355273740930133e-17 } }

    #[inline]
    fn LN_10() -> Self { Self { hi: 2.302585092994046, lo: -2.170756223382255e-16 } }

    #[inline]
    fn LN_2() -> Self { Self { hi: 0.6931471805599453, lo: 2.31904681384628e-17 } }

    #[inline]
    fn LOG10_2() -> Self { Self { hi: 0.3010299956639812, lo: -2.803728127785059e-18 } }

    #[inline]
    fn LOG10_E() -> Self { Self { hi: 0.4342944819032518, lo: 1.0983196502167326e-17 } }

    #[inline]
    fn LOG2_10() -> Self { Self { hi: 3.321928094887362, lo: 1.6616175169735967e-16 } }

    #[inline]
    fn FRAC_1_PI() -> Self { Self { hi: 0.3183098861837907, lo: -1.9678676675182615e-17 } }

    #[inline]
    fn FRAC_1_SQRT_2() -> Self { Self { hi: 0.7071067811865476, lo: -4.833646656726404e-17 } }

    #[inline]
    fn FRAC_2_PI() -> Self { Self { hi: 0.6366197723675814, lo: -3.935735335036523e-17 } }

    #[inline]
    fn FRAC_2_SQRT_PI() -> Self { Self { hi: 1.1283791670955126, lo: 1.5335459613166276e-17 } }

    #[inline]
    fn FRAC_PI_2() -> Self { Self { hi: 1.5707963267948966, lo: 6.123233995736757e-17 } }

    #[inline]
    fn FRAC_PI_3() -> Self { Self { hi: 1.0471975511965979, lo: -1.0720817664510863e-16 } }

    #[inline]
    fn FRAC_PI_4() -> Self { Self { hi: 0.7853981633974483, lo: 3.061616997868379e-17 } }

    #[inline]
    fn FRAC_PI_6() -> Self { Self { hi: 0.5235987755982989, lo: -5.3604088322554315e-17 } }

    #[inline]
    fn FRAC_PI_8() -> Self { Self { hi: 0.39269908169872414, lo: 1.5308084989341894e-17 } }

    #[inline]
    fn SQRT_2() -> Self { Self { hi: 1.4142135623730951, lo: -9.667293313452965e-17 } }

    #[inline]
    fn TAU() -> Self { Self { hi: 6.283185307179586, lo: 2.449293598294703e-16 } }
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
    // TODO: might over flow
    #[inline]
    fn from_i64(n: i64) -> Option<Self> { Some(Self { hi: n as f64, lo: 0.0 }) }

    // TODO: might overflow
    #[inline]
    fn from_u64(n: u64) -> Option<Self> { Some(Self { hi: n as f64, lo: 0.0 }) }

    #[inline]
    fn from_f64(n: f64) -> Option<Self> { Some(Self { hi: n, lo: 0.0 }) }
}

const COEFFS: [Complex<f64x2>; 24] = [
    Complex::<f64x2> {
        re: f64x2 { hi: 1.0, lo: -4.098967624221445e-34 },
        im: f64x2 { hi: -4.098967624221445e-34, lo: -1.5555341444609406e-50 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: 3.013639102622782e-31, lo: 1.955165282442337e-48 },
        im: f64x2 { hi: 1.0, lo: 2.9978991830669402e-31 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: -0.5, lo: -3.686440395906948e-29 },
        im: f64x2 { hi: -3.647933289584011e-29, lo: -1.7688387812116778e-45 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: 1.7944312320610667e-27, lo: 1.204247011048231e-43 },
        im: f64x2 { hi: -0.16666666666666666, lo: -9.251858536776753e-18 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: 0.041666666666666664, lo: 2.3129645882491033e-18 },
        im: f64x2 { hi: -4.540745872903123e-26, lo: 2.3501738525509453e-43 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: 7.369657246850729e-25, lo: 1.3905305042833041e-41 },
        im: f64x2 { hi: 0.008333333333333333, lo: 1.1564894904668587e-19 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: -0.001388888888888889, lo: 5.2997583852953265e-20 },
        im: f64x2 { hi: -7.600749353559733e-24, lo: -1.4295452929877914e-40 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: 5.952935691552932e-23, lo: -4.8473575958472206e-39 },
        im: f64x2 { hi: -0.0001984126984126984, lo: -1.148609310509542e-22 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: 2.48015873015873e-5, lo: -3.124037069753179e-22 },
        im: f64x2 { hi: -3.1887890759822763e-22, lo: -5.0724289498974544e-39 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: 1.4275560892801856e-21, lo: 1.7937983076802776e-38 },
        im: f64x2 { hi: 2.7557319223985905e-6, lo: -1.03072795026984e-22 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: -2.7557319223986365e-7, lo: -1.1994468685385469e-23 },
        im: f64x2 { hi: -4.469574115662231e-21, lo: -1.4269807213263557e-37 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: 1.2519942771424863e-20, lo: 6.493343756430319e-37 },
        im: f64x2 { hi: -2.505210838543005e-8, lo: 1.0110223535016303e-24 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: 2.0876756987604126e-9, lo: 1.596410632745976e-25 },
        im: f64x2 { hi: -2.4358528608849723e-20, lo: -4.7764500560843e-37 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: 4.4900371785054185e-20, lo: -1.7502963760208312e-36 },
        im: f64x2 { hi: 1.6059043840917783e-10, lo: -2.0035561424541567e-27 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: -1.1470745659638183e-11, lo: -1.6244361659394653e-28 },
        im: f64x2 { hi: -5.573373427988895e-20, lo: -1.4956134151510146e-36 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: 6.93156930534892e-20, lo: -2.4015923876407123e-36 },
        im: f64x2 { hi: -7.647163117497653e-13, lo: -4.5278084264034546e-29 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: 4.779471038020406e-14, lo: -1.1886747323088784e-30 },
        im: f64x2 { hi: -5.474019274931799e-20, lo: -2.473921669251002e-36 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: 4.614348173584815e-20, lo: 2.0776561890169216e-37 },
        im: f64x2 { hi: 2.811496455585735e-15, lo: -7.581326552442112e-32 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: -1.562191434065831e-16, lo: 1.8968334145071065e-34 },
        im: f64x2 { hi: -2.2322810711198874e-20, lo: -1.17247052214804e-36 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: 1.2537156366849241e-20, lo: -3.4871226451665094e-37 },
        im: f64x2 { hi: -8.210701448480838e-18, lo: -6.658269703479169e-34 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: 4.0654755277215486e-19, lo: -4.773989908786719e-36 },
        im: f64x2 { hi: -3.359690826492369e-21, lo: -1.40728387339678e-37 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: 1.19777004853017e-21, lo: -6.132839002942767e-38 },
        im: f64x2 { hi: 2.0397006812988565e-20, lo: 8.748308795144894e-37 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: -1.1156950768713551e-21, lo: -2.2840707511532131e-38 },
        im: f64x2 { hi: -1.340923713021694e-22, lo: -3.6414216793947587e-39 },
    },
    Complex::<f64x2> {
        re: f64x2 { hi: 2.7169899607602477e-23, lo: 8.768090448822874e-40 },
        im: f64x2 { hi: -2.7169899607602477e-23, lo: -8.768090448822874e-40 },
    },
];

impl ExpPolyApprox for f64x2 {
    type Output = impl Iterator<Item = (usize, Complex<Self>)>;

    #[inline]
    fn get_poly_approx() -> Self::Output { COEFFS.iter().enumerate().map(|(idx, &x)| (idx, x)) }
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
    fn erfc(self, eps: f64) -> Self {
        if self.hi.abs() > 2.0 {
            let z = self;
            let h = f64::PI() / (6.0 / eps).ln().sqrt();
            let K = ((1.0 / eps).ln().sqrt() / h).ceil() as i32;

            let z_sq = z.square();
            let mut ret = Self::one() / z_sq;

            let h = h.unchecked_cast::<Self>();
            let h_sq = h.square();
            for k in 1..=K {
                let w = h_sq * (k * k).unchecked_cast::<Self>();
                ret += 2.0.unchecked_cast::<Self>() * (-w).exp() / (z_sq + w);
            }
            ret * (-z_sq).exp() * h * z / Self::PI()
                + 2.0.unchecked_cast::<Self>() / (Self::one() - (Self::TAU() * z / h).exp())
        } else {
            let s;
            let z;
            if self.hi >= 0.0 {
                s = 1;
                z = self;
            } else {
                s = -1;
                z = -self;
            }

            let eps0 = eps / 2.0;

            let mut t = z;
            let mut k = 0i32;
            let mut S = Self::zero();
            let z_sq = z.square();
            loop {
                let ds = t / (2 * k + 1).unchecked_cast::<Self>();
                S += ds;
                if ds.abs().hi < eps0 {
                    break;
                }
                k += 1;
                t *= -z_sq / k.unchecked_cast::<Self>();
            }
            if s == 1 {
                Self::one() - S * Self::FRAC_2_SQRT_PI()
            } else {
                Self::one() + S * Self::FRAC_2_SQRT_PI()
            }
        }
    }
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
