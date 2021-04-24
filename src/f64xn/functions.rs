use crate::unchecked_cast::UncheckedCast;

use super::{blocks::*, f64x2};
use num::traits::FloatConst;
use num::{Float, One, Zero};
use num_traits::AsPrimitive;

#[inline]
fn pow2(n: i32) -> f64 { f64::from_bits(((0x3ff + n) as u64) << 52) }

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
        let exp_factor2 = pow2(factor2 as i32);
        let exp_rem = rem.expm1_remez() + 1.0;
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

    pub fn ln_newton(self) -> Self {
        let y = Self::from(self.hi.ln());
        let y = y + self * (-y).exp() - 1.0;
        let ret = y + self * (-y).exp() - 1.0;

        let diff = self.ln() - ret;
        assert!(-1e-20 <= diff.hi && diff.hi <= 1e-20);

        ret
    }

    pub fn ln(self) -> Self {
        let bits = self.hi.to_bits();
        let mut e = ((bits >> 52) & 0x7ff) as i32 - 0x3ff;
        let mantissa = bits & 0xfffffffffffff;
        if mantissa > 1865452045155277u64 {
            // 2.0f64.sqrt().to_bits() & 0xfffffffffffff
            // in this case: self.hi.log2().round() = e + 1;
            e += 1;
        }

        let p2 = pow2(e);
        let normalized = Self { hi: self.hi / p2, lo: self.lo / p2 };
        normalized.ln_remez() + e.unchecked_cast::<Self>() * Self::LN_2()
    }

    pub fn cos(self) -> Self {
        let x = (self / Self::FRAC_PI_2()).floor();
        let y = self - x * Self::FRAC_PI_2();
        let x: i64 = x.unchecked_cast();
        match x.rem_euclid(4) {
            0 => y.cos_remez(),
            1 => -y.sin_remez(),
            2 => -y.cos_remez(),
            3 => y.sin_remez(),
            _ => unreachable!(),
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
            3 => -y.cos_remez(),
            _ => unreachable!(),
        }
    }

    /// See https://github.com/ajtribick/twofloat/blob/master/src/functions/power.rs#L30-L40
    #[inline]
    pub fn sqrt(self) -> Self {
        if self.is_zero() {
            return self;
        }
        assert!(self.hi >= 0.0, "self = {:?}", self);
        let x = self.hi.sqrt().recip();
        let y = self.hi * x;
        two_add(y, (self - two_mul(y, y)).hi * (x * 0.5))
    }

    #[inline]
    pub fn recip(self) -> Self { 1.0 / self }

    #[inline]
    pub fn hypot(self, other: Self) -> Self { (self.square() + other.square()).sqrt() }

    pub fn exp_m1(self) -> Self {
        if -0.34657359028 <= self.hi && self.hi <= 0.34657359028 {
            self.expm1_remez()
        } else {
            self.exp() - 1.0
        }
    }

    #[inline(never)]
    pub fn cosh(self) -> Self {
        let t = self.exp();
        (t + t.recip()) * 0.5
    }

    #[inline(never)]
    pub fn sinh(self) -> Self {
        let threshold = 0.34657359028; // ln(2) / 2
        if -threshold <= self.hi && self.hi <= threshold {
            let t = self.expm1_remez();
            t - t * t / (t + 1.0) * 0.5
            // (t + t / (t + 1.0)) * 0.5
        } else {
            let t = self.exp();
            (t - t.recip()) * 0.5
        }
    }

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

    pub fn ceil(self) -> Self { -(-self).floor() }

    pub fn round(self) -> Self {
        let ret = two_add_fast(self.hi.round(), self.lo.round());
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
                    (base, tan_base) = (
                        0.19634954084936207,
                        f64x2 { hi: 0.198912367379658, lo: -7.117703886485398e-18 },
                    );
                } else {
                    // PI/8 <= atan(x) < PI/4 => base = 3 PI / 16
                    (base, tan_base) = (
                        0.5890486225480862,
                        f64x2 { hi: 0.6681786379192989, lo: 7.828409495095202e-18 },
                    );
                }
            } else {
                // another division
                if approx < 2.414213562373095 {
                    // tan(3 PI / 8)
                    // PI/4 <= atan(x) < 3PI/8 => base = 5PI / 16
                    (base, tan_base) = (
                        0.9817477042468103,
                        f64x2 { hi: 1.496605762665489, lo: -5.424792800215555e-17 },
                    );
                } else {
                    // 3 PI / 8 <= atan(x) < PI / 2 => base = 7 PI / 16
                    (base, tan_base) = (
                        1.3744467859455345,
                        f64x2 { hi: 5.027339492125846, lo: 3.9817094207573838e-16 },
                    );
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

impl f64x2 {
    fn erfc_1(self) -> Self { self.erfc_1_remez() + 1.0 }

    fn erfc_1_5(self) -> Self {
        debug_assert!(self.hi >= 0.0);
        self.erfc_1_5_remez() / self.square().exp()
    }

    fn erfc_5_10(self) -> Self {
        debug_assert!(self.hi >= 0.0);
        self.erfc_5_10_remez() / self.square().exp() / self
    }

    fn erfc_eps_large(self, eps: f64) -> Self {
        let z = self;
        let h = f64::PI() / (6.0 / eps).ln().sqrt();
        let K = ((1.0 / eps).ln().sqrt() / h).ceil() as i32;

        let z_sq = z.square();
        let mut ret = Self::one() / z_sq;

        unsafe {
            crate::profiler::ERFC_EPS_LARGE_K += K as usize;
            crate::profiler::ERFC_EPS_LARGE_CNT += 1;
        }

        let h = h.unchecked_cast::<Self>();
        let h_sq = h.square();
        for k in 1..=K {
            let w = h_sq * (k * k).unchecked_cast::<Self>();
            ret += (-w).exp() * 2.0 / (z_sq + w);
        }
        ret * (-z_sq).exp() * h * z / Self::PI() + 2.0 / (Self::one() - (Self::TAU() * z / h).exp())
    }

    fn erfc_eps_small(self, eps: f64) -> Self {
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
            let ds = t / (2 * k + 1) as f64;
            S += ds;
            if ds.abs().hi < eps0 {
                break;
            }
            k += 1;
            t *= -z_sq / k as f64;
        }

        unsafe {
            crate::profiler::ERFC_EPS_SMALL_CNT += 1;
            crate::profiler::ERFC_EPS_SMALL_K += k as usize;
        }

        if s == 1 {
            Self::one() - S * Self::FRAC_2_SQRT_PI()
        } else {
            Self::one() + S * Self::FRAC_2_SQRT_PI()
        }
    }

    fn erfc_eps_old(self, eps: f64) -> Self {
        if self.hi.abs() < 2.0 {
            self.erfc_eps_small(eps)
        } else {
            self.erfc_eps_large(eps)
        }
    }

    pub fn erfc_eps(self, eps: f64) -> Self {
        let x = self;
        let x_abs = x.abs();

        let ret = if x_abs.hi < 1.0 {
            self.erfc_1()
        } else {
            let y;
            if x_abs.hi < 2.0 {
                y = x_abs.erfc_eps_small(eps)
            } else if x_abs.hi < 5.0 {
                y = x_abs.erfc_1_5()
            } else if x_abs.hi < 10.0 {
                y = x_abs.erfc_5_10()
            } else {
                y = Self::zero()
            }

            if x.hi > 0.0 {
                y
            } else {
                Self { hi: 2.0, lo: 0.0 } - y
            }
        };

        // let check = x.erfc_eps_old(1e-33);
        // assert!((ret - check).abs().hi < 1e-30, "x = {}, ret = {}, check = {}", x, ret, check);

        ret
    }

    fn erfc(self) -> Self { todo!() }
}

#[cfg(test)]
mod tests {}
