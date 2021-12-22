use crate::{f64xn::f64x, traits::FpOps};
use num::{Float, One, Zero};
use std::ops::{Add, Mul};

impl<const N: usize> Zero for f64x<N>
where
    f64x<N>: Add<Output = Self>,
{
    #[inline]
    fn zero() -> Self { Self { data: [0.0; N] } }

    #[inline]
    fn is_zero(&self) -> bool { *self == Self::ZERO }
}

impl<const N: usize> One for f64x<N>
where
    f64x<N>: Mul<Output = Self>,
{
    #[inline]
    fn one() -> Self { Self::ONE }

    #[inline]
    fn is_one(&self) -> bool { *self == Self::ONE }
}

impl<const N: usize> From<f64x<N>> for String
where
    f64x<N>: Float + FpOps,
{
    fn from(a: f64x<N>) -> Self {
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
            let mut a = a;
            let mut ret = String::from("");
            if a.is_sign_negative() {
                a = -a;
                ret.push('-');
            }
            let mut e = 0;
            while a / 10.0 >= f64x::<N>::one() {
                a = a / 10.0; // TODO: avoid division
                e += 1;
            }
            while a < f64x::<N>::one() {
                a = a * 10.0;
                e -= 1;
            }

            let mut dec_point = false;

            for _ in 0..30 {
                let d = a.floor().fp(); // a is in [0, 8] so it's fine to use floor
                ret.push((b'0' + d as u8) as char);
                a = (a - d) * 10.;
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

impl<const N: usize> std::fmt::Display for f64x<N>
where
    f64x<N>: Into<String>,
{
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", (*self).into())
    }
}
