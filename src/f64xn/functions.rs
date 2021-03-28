use crate::unchecked_cast::UncheckedCast;

use super::{blocks::*, f64x2};
use num::traits::FloatConst;
use num::{Float, One, Zero};
use num_traits::AsPrimitive;

impl f64x2 {
    /// see https://github.com/JuliaMath/DoubleFloats.jl/blob/master/src/math/ops/op_dd_dd.jl#L45-L52
    #[inline]
    pub fn square(self) -> Self {
        let a = self;
        let Self { hi, lo } = two_mul(a.hi, a.hi);
        let lo = lo + (a.hi * a.lo) * 2.;
        let lo = lo + a.lo * a.lo;
        Self { hi, lo }
    }

    fn modf(self) -> (Self, Self) {
        let f_hi = self.hi.floor();
        let f_lo = self.lo.floor();
        let i_hi = self.hi.fract();
        let i_lo = self.lo.fract();
        let i = two_add(i_hi, i_lo);
        let f = two_add(f_hi, f_lo);
        (f, i)
    }

    pub fn exp(self) -> Self {
        let factor2 = (self / Self::LN_2()).round();
        let rem = self - Self::LN_2() * factor2;
        let factor2: i64 = factor2.as_();
        let exp_factor2 = 2.0.powi(factor2 as i32);
        let exp_rem = rem.exp_remez();
        Self { hi: exp_rem.hi * exp_factor2, lo: exp_rem.lo * exp_factor2 }
    }

    pub fn powi(self, n: i64) -> Self {
        if n == 0 {
            Self::one()
        } else {
            let mut s = Self::one();
            let mut n_abs = n.abs();
            let mut r = self;
            if n_abs > 1 {
                while n_abs > 0 {
                    if n_abs % 2 == 1 {
                        s *= r;
                    }
                    n_abs >>= 1;
                    if n_abs > 0 {
                        r = r.square();
                    }
                }
            } else {
                s = r;
            }
            if n < 0 {
                s.recip()
            } else {
                s
            }
        }
    }

    pub fn ln(self) -> Self {
        let y = Self::from(self.hi.ln());
        let y = y + self * (-y).exp() - 1.0;
        y + self * (-y).exp() - 1.0
    }

    pub fn cos(self) -> Self {
        let x = (self / Self::FRAC_PI_2()).floor();
        let y = self - x * Self::FRAC_PI_2();
        let x: i64 = x.unchecked_cast();
        match x.rem_euclid(4) {
            0 => y.cos_remez(),
            1 => -y.sin_remez(),
            2 => -y.cos_remez(),
            _ => y.sin_remez(),
        }
    }

    #[inline]
    pub fn sin(self) -> Self {
        let x = (self / Self::FRAC_PI_2()).floor();
        let y = self - x * Self::FRAC_PI_2();
        let x: i64 = x.unchecked_cast();
        match x.rem_euclid(4) {
            0 => y.sin_remez(),
            1 => y.cos_remez(),
            2 => -y.sin_remez(),
            _ => -y.cos_remez(),
        }
    }

    /// See https://github.com/ajtribick/twofloat/blob/master/src/functions/power.rs#L30-L40
    #[inline]
    pub fn sqrt(self) -> Self {
        assert!(self.hi >= 0.0, "self = {:?}", self);
        let x = self.hi.sqrt().recip();
        let y = self.hi * x;
        two_add(y, (self.hi - two_mul(y, y).hi) * (x * 0.5))
    }

    #[inline]
    pub fn recip(self) -> Self { 1.0 / self }

    #[inline]
    pub fn hypot(self, other: Self) -> Self { (self.square() + other.square()).sqrt() }

    #[inline]
    pub fn cosh(self) -> Self { (self.exp() + (-self).exp()) / 2.0 }

    #[inline]
    pub fn sinh(self) -> Self { (self.exp() - (-self).exp()) / 2.0 }

    pub fn floor(self) -> Self {
        let ret = if self.hi.fract() == 0.0 {
            two_add_fast(self.hi, self.lo.floor())
        } else {
            Self { hi: self.hi.floor(), lo: 0.0 }
        };
        // let diff = self - ret;
        // assert!(0.0 <= diff.hi && diff.hi <= 1.0, "self = {:?}", self);
        ret
    }

    pub fn ceil(self) -> Self {
        -(-self).floor()
        // let ret = two_add_fast(self.hi.ceil(), self.lo.ceil());
        // let ret;
        // if self.lo >= 0.5 || self.lo <= -0.5 {
        //     // see analysis for floor
        //     ret = two_add_fast(self.hi, self.lo.ceil())
        // } else {
        //     // hi.ceil doesn't overflow.
        //     ret = Self { hi: self.hi.ceil(), lo: 0.0 }
        // }
        // let diff = ret - self;
        // assert!(0.0 <= diff.hi && diff.hi < 1.0, "self = {:?}", self);
        // ret
    }

    pub fn round(self) -> Self {
        let ret = two_add_fast(self.hi.round(), self.lo.round());
        // let ret;
        // if self.lo >= 0.5 || self.lo <= -0.5 {
        //     ret = two_add_fast(self.hi, self.lo.round())
        // } else {
        //     ret = Self { hi: self.hi.round(), lo: 0.0 }
        // }
        // let diff = ret - self;
        // assert!(-0.5 <= diff.hi && diff.hi <= 0.5, "self = {:?}", self);
        ret
    }

    /// The idea is similar to libm: We split [0, \pi/2) into 4 parts, and approximate atan(x) in each
    /// libm uses 7/16 for thresholding. I guess its due to the precision limit of variable `ix`
    /// here we just partition them equally.
    pub fn atan(self) -> Self {
        let approx = self.hi;
        if approx < 0.0 {
            -(-self).atan()
        } else {
            let base;
            let tan_base;

            if approx <= 1.0 {
                // atan(x) < PI/4
                if approx < 0.41421356237309503 {
                    // tan(PI / 8)
                    // atan(x) < PI/8 => base = PI / 16
                    (base, tan_base) = (0.19634954084936207, 0.198912367379658);
                } else {
                    // PI/8 <= atan(x) < PI/4 => base = 3 PI / 16
                    (base, tan_base) = (0.5890486225480862, 0.6681786379192989);
                }
            } else {
                // another division
                if approx < 2.414213562373095 {
                    // tan(3 PI / 8)
                    // PI/4 <= atan(x) < 3PI/8 => base = 5PI / 16
                    (base, tan_base) = (0.9817477042468103, 1.496605762665489);
                } else {
                    // 3 PI / 8 <= atan(x) < PI / 2 => base = 7 PI / 16
                    (base, tan_base) = (1.3744467859455345, 5.027339492125846);
                }
            }

            base + ((self - tan_base) / (1.0 + self * tan_base)).atan_remez()
        }
    }

    pub fn atan2(self, other: Self) -> Self {
        // y = self, x = rhs
        if self.hi < 0.0 {
            // y < 0
            -(-self).atan2(other)
        } else if self.is_zero() {
            // y == 0
            if other.is_sign_positive() {
                // x > 0, y = 0
                Self::zero()
            } else {
                // x <= 0, y = 0
                Self::PI()
            }
        } else if other.is_sign_negative() {
            // y > 0, x < 0
            Self::PI() - (-self / other).atan()
        } else if other.is_sign_positive() {
            // y > 0, x > 0
            (self / other).atan()
        } else {
            // y > 0, x == 0
            Self::FRAC_PI_2()
        }
    }

    #[inline]
    pub fn powf(self, other: Self) -> Self { (self.ln() * other).exp() }

    pub fn abs(self) -> Self {
        if self.hi < 0.0 {
            -self
        } else {
            self
        }
    }
}
