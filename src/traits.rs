use num::{Complex, Signed};
use num_traits::{Float, FloatConst, NumAssignOps, Pow};
use rustfft::FftNum;
use std::{
    fmt::{Debug, Display, LowerExp},
    ops::{Add, Div, Mul, Sub},
};

use crate::{
    constants::EXP_POLY_EXP_F64,
    unchecked_cast::{UncheckedCast, UncheckedFrom, UncheckedInto},
};

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

use F64x2::f64x2;

use crate::constants::RS_GABCKE_GAMMA;

const F64X2_COEFFS: [Complex<f64x2>; 28] = [
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
    fn get_poly_approx() -> Self::Output {
        F64X2_COEFFS.iter().enumerate().map(|(idx, &x)| (idx, x))
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
