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

    // |x| <= 0.34657359028 = ln(2) / 2
    fn exp_remez(self) -> Self {
        const C2: f64x2 = f64x2 { hi: 0.16666666666666666, lo: 9.251858538542447e-18 };
        const C4: f64x2 = f64x2 { hi: -0.002777777777777778, lo: 1.0601087929995308e-19 };
        const C6: f64x2 = f64x2 { hi: 6.613756613756614e-5, lo: -4.460173646997389e-21 };
        const C8: f64x2 = f64x2 { hi: -1.6534391534391535e-6, lo: 7.121962972677988e-23 };
        const C10: f64x2 = f64x2 { hi: 4.1753513975736114e-8, lo: 1.158547249215353e-24 };
        const C12: f64x2 = f64x2 { hi: -1.0568380277354652e-9, lo: 5.58404280005523e-26 };
        const C14: f64x2 = f64x2 { hi: 2.6765073029312422e-11, lo: 9.12829168899536e-28 };
        const C16: f64x2 = f64x2 { hi: -6.779357355296518e-13, lo: -1.936460166148485e-29 };
        const C18: f64x2 = f64x2 { hi: 1.71700970829093e-14, lo: -1.033795194722353e-30 };
        const C20: f64x2 = f64x2 { hi: -4.278007864569552e-16, lo: 1.5907971763038017e-32 };

        let x = self;
        let x2 = x * x;
        let x4 = x2 * x2;
        let x6 = x2 * x4;
        let x8 = x4 * x4;
        let x14 = x8 * x6;

        let r1 = x2 * C2 + x4 * C4 + x6 * C6;
        let r2 = x8 * (C8 + x2 * C10 + x4 * C12);
        let r3 = x14 * (C16 + x2 * C18 + x4 * C20);
        // r = 2.0 + r1 + r2 + r3;
        // r = x (exp(x) + 1) / (exp(x) - 1) => exp(x) = 1 + 2r / (r - x)
        let c = x - (r1 + r2 + r3);
        1.0 + x + x * c / (2.0 - c)
    }

    pub fn exp(self) -> Self {
        let factor2 = (self / Self::LN_2()).round();
        let rem = self - Self::LN_2() * factor2;
        let factor2: i64 = factor2.as_();
        let exp_factor2 = 2.0.powi(factor2 as i32);
        let exp_rem = rem.exp_remez();
        Self { hi: exp_rem.hi * exp_factor2, lo: exp_rem.lo * exp_factor2 }
    }

    fn powi(self, n: i64) -> Self {
        if n == 0 {
            Self::zero()
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
        let y = y + self * y.exp() - 1.0;
         y + self * y.exp() - 1.0
    }

    pub fn cos(self) -> Self {
        let x = (self / Self::FRAC_PI_2()).floor();
        let y = self - x * Self::FRAC_PI_2();
        let x: i64 = x.as_();
        match x.rem_euclid(4) {
            0 => y.cos_remez(),
            1 => -y.sin_remez(),
            2 => -y.cos_remez(),
            _ => y.sin_remez(),
        }
    }

    pub fn sin(self) -> Self {
        let x = (self / Self::FRAC_PI_2()).floor();
        let y = self - x * Self::FRAC_PI_2();
        let x: i64 = x.as_();
        match x.rem_euclid(4) {
            0 => y.sin_remez(),
            1 => y.cos_remez(),
            2 => -y.sin_remez(),
            _ => -y.cos_remez(),
        }
    }

    #[allow(clippy::float_cmp)]
    pub fn powf(self, n: f64) -> Self {
        if n.floor() == n {
            self.powi(n as i64)
        } else {
            (self.ln() * n).exp()
        }
    }

    // 1 / sqrt(2) < self < sqrt(2)
    fn log_remez(self) -> Self {
        const C1: f64x2 = f64x2 { hi: 2.0, lo: 2.531693403050348e-32 };
        const C3: f64x2 = f64x2 { hi: 0.6666666666666666, lo: 3.700743415403453e-17 };
        const C5: f64x2 = f64x2 { hi: 0.4, lo: -2.220446027080773e-17 };
        const C7: f64x2 = f64x2 { hi: 0.2857142857142857, lo: 1.586016141254861e-17 };
        const C9: f64x2 = f64x2 { hi: 0.2222222222222222, lo: 1.2407739663861083e-17 };
        const C11: f64x2 = f64x2 { hi: 0.1818181818181818, lo: 3.2065997154550064e-18 };
        const C13: f64x2 = f64x2 { hi: 0.1538461538461574, lo: -3.1965060154741107e-18 };
        const C15: f64x2 = f64x2 { hi: 0.13333333333287886, lo: 7.44820907006051e-18 };
        const C17: f64x2 = f64x2 { hi: 0.11764705886515148, lo: -1.17519451979169e-18 };
        const C19: f64x2 = f64x2 { hi: 0.1052631551291071, lo: 4.4806067994093946e-18 };
        const C21: f64x2 = f64x2 { hi: 0.095238228662164, lo: -3.577683177336467e-18 };
        const C23: f64x2 = f64x2 { hi: 0.08695190137972356, lo: -5.409228774321587e-18 };
        const C25: f64x2 = f64x2 { hi: 0.08011166789102804, lo: -2.092932986143882e-18 };
        const C27: f64x2 = f64x2 { hi: 0.07229374987314603, lo: -4.021248728080285e-18 };
        const C29: f64x2 = f64x2 { hi: 0.0855816562700506, lo: 6.675796100148133e-19 };

        let s = self;
        let x = (s - 1.0) / (s + 1.0);
        let x2 = x * x;
        let x4 = x2 * x2;
        let x6 = x2 * x4;
        let x8 = x4 * x4;
        let x10 = x4 * x6;
        let x20 = x10 * x10;

        let r1 = C1 + x2 * C3 + x4 * C5 + x6 * C7 + x8 * C9;
        let r2 = x10 * (C11 + x2 * C11 + x4 * C11 + x6 * C15 + x8 * C19);
        let r3 = x20 * (C21 + x2 * C23 + x4 * C25 + x6 * C27 + x8 * C29);

        x * (r1 + r2 + r3)
    }

    // |x| < pi / 4
    fn cos_remez(self) -> Self {
        const C0: f64x2 = f64x2 { hi: 1.0, lo: -5.795385665153811e-34 };
        const C2: f64x2 = f64x2 { hi: -0.5, lo: 2.7060102269624583e-31 };
        const C4: f64x2 = f64x2 { hi: 0.041666666666666664, lo: 2.3129646346148305e-18 };
        const C6: f64x2 = f64x2 { hi: -0.001388888888888889, lo: 5.3005440176619425e-20 };
        const C8: f64x2 = f64x2 { hi: 2.48015873015873e-5, lo: 2.1502053426601327e-23 };
        const C10: f64x2 = f64x2 { hi: -2.755731922398589e-7, lo: -2.367645329977107e-23 };
        const C12: f64x2 = f64x2 { hi: 2.087675698786809e-9, lo: 1.7286815346236624e-25 };
        const C14: f64x2 = f64x2 { hi: -1.1470745597727671e-11, lo: -7.352837032889324e-29 };
        const C16: f64x2 = f64x2 { hi: 4.77947733186016e-14, lo: 3.0400747811634378e-30 };
        const C18: f64x2 = f64x2 { hi: -1.5619206074477628e-16, lo: -7.726721637678806e-33 };
        const C20: f64x2 = f64x2 { hi: 4.110221447608981e-19, lo: -1.9180683347243325e-35 };
        const C22: f64x2 = f64x2 { hi: -8.837331187274062e-22, lo: -6.189968272642389e-38 };

        let x = self;
        let x2 = x * x;
        let x4 = x2 * x2;
        let x6 = x2 * x4;
        let x8 = x4 * x4;
        let x16 = x8 * x8;

        let r1 = C0 + x2 * C2 + x4 * C4 + x6 * C6;
        let r2 = x8 * (C8 + x2 * C10 + x4 * C12 + x6 * C14);
        let r3 = x16 * (C16 + x2 * C18 + x4 * C20 + x6 * C22);

        r1 + r2 + r3
    }

    // |x| < pi / 4
    fn sin_remez(self) -> Self {
        const C1: f64x2 = f64x2 { hi: 1.0, lo: -2.8984825337707473e-34 };
        const C3: f64x2 = f64x2 { hi: -0.16666666666666666, lo: -9.251858538542923e-18 };
        const C5: f64x2 = f64x2 { hi: 0.008333333333333333, lo: 1.1564823172934677e-19 };
        const C7: f64x2 = f64x2 { hi: -0.0001984126984126984, lo: -1.7209552641261292e-22 };
        const C9: f64x2 = f64x2 { hi: 2.7557319223985893e-6, lo: -1.8584006050460325e-22 };
        const C11: f64x2 = f64x2 { hi: -2.505210838544172e-8, lo: 1.4546921392186555e-24 };
        const C13: f64x2 = f64x2 { hi: 1.605904383682161e-10, lo: 7.644272326832948e-27 };
        const C15: f64x2 = f64x2 { hi: -7.647163731818732e-13, lo: -4.7818697387989105e-29 };
        const C17: f64x2 = f64x2 { hi: 2.8114572540870227e-15, lo: -1.8417162789441038e-31 };
        const C19: f64x2 = f64x2 { hi: -8.22063483478287e-18, lo: 1.7103320315576008e-34 };
        const C21: f64x2 = f64x2 { hi: 1.9572521175911013e-20, lo: -5.661527325082782e-38 };
        const C23: f64x2 = f64x2 { hi: -3.843391876512855e-23, lo: -2.614405314562012e-40 };

        let x = self;
        let x2 = x * x;
        let x4 = x2 * x2;
        let x6 = x2 * x4;
        let x8 = x4 * x4;
        let x16 = x8 * x8;

        let r1 = C1 + x2 * C3 + x4 * C5 + x6 * C7;
        let r2 = x8 * (C9 + x2 * C11 + x4 * C11 + x6 * C15);
        let r3 = x16 * (C17 + x2 * C19 + x4 * C21 + x6 * C23);

        x * (r1 + r2 + r3)
    }

    /// See https://github.com/ajtribick/twofloat/blob/master/src/functions/power.rs#L30-L40
    fn sqrt(self) -> Self {
        assert!(self.hi >= 0.0);
        let x = self.hi.sqrt().recip();
        let y = self.hi * x;
        two_add(y, (self.hi - two_mul(y, y).hi) * (x * 0.5))
    }

    fn recip(self) -> Self { 1.0 / self }

    fn hypot(self, other: Self) -> Self { (self.square() + other.square()).sqrt() }

    fn cosh(self) -> Self { (self.exp() + (-self).exp()) / 2.0 }

    fn sinh(self) -> Self { (self.exp() - (-self).exp()) / 2.0 }
}
