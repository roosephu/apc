use crate::{f64xn::f64x, traits::FpOps};
use num::{Float, One, Zero};
use num_traits::NumOps;
use std::ops::{Add, Div, Mul, Sub};

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

impl<const N: usize> FpOps for f64x<N>
where
    f64x<N>: Add<f64, Output = Self>
        + Sub<f64, Output = Self>
        + Mul<f64, Output = Self>
        + Div<f64, Output = Self>,
{
    fn fp(&self) -> f64 { Self::fp(self) }

    fn mp(x: f64) -> Self { Self::mp(x) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand;
    use rug::{ops::CompleteRound, Float};

    const N: usize = 3;
    const PREC: u32 = N as u32 * 53;
    type f64xn = f64x<N>;

    fn assert_rug_close(a: &f64xn, b: &Float, rtol: f64, atol: f64) {
        if b.to_f64() == 0.0 || b.is_nan() || a.data[0].is_nan() {
            return;
        }
        let diff = f64xn_to_float(*a) - b;
        let abs = diff.to_f64().abs();
        let rel = abs / b.to_f64().abs();
        assert!(rel <= rtol || abs <= atol, "rel = {:.3e}, abs = {:.3e}", rel, abs);
    }

    fn rand_f64xn(k: usize, p: f64) -> (f64x<N>, Float) {
        let mut x = f64x::<N>::zero();
        let mut rug_float = Float::with_val(PREC, 0.0f64);
        for _ in 0..k {
            let mut bits = 0x7ff0000000000000u64;
            while (bits & 0x7ff0000000000000u64) == 0x7ff0000000000000u64 {
                bits = 0;
                for i in 0..62 {
                    bits |= ((rand::random::<f64>() < p) as u64) << i;
                }
            }
            let f = f64::from_bits(bits) * (if rand::random::<u8>() % 2 == 0 { 1.0 } else { -1.0 });
            assert_rug_close(&(x + f), &(rug_float.clone() + f), 1e-40, 1e-300);
            assert_rug_close(&(x - f), &(rug_float.clone() - f), 1e-40, 1e-300);
            assert_rug_close(&(x * f), &(rug_float.clone() * f), 1e-40, 1e-300);
            assert_rug_close(&(x / f), &(rug_float.clone() / f), 1e-40, 1e-300 / f.abs());

            x = x + f;

            rug_float += f;
        }

        (x, rug_float)
    }

    fn f64xn_to_float(x: f64x<N>) -> Float {
        let mut float = Float::with_val(PREC, 0.0);
        for i in 0..N {
            float += x.data[i];
        }
        float
    }

    #[test]
    fn basic_ops() {
        const M: usize = 10000;
        // creation, f64xn + f64
        for _ in 0..M {
            rand_f64xn(100, 0.7);
        }

        // f64xn + f64xn
        for _ in 0..M {
            let a = rand_f64xn(10, 0.8);
            let b = rand_f64xn(10, 0.8);
            assert_rug_close(&(a.0 + b.0), &(a.1.clone() + &b.1), 1e-40, 1e-300);
            assert_rug_close(&(a.0 - b.0), &(a.1.clone() - &b.1), 1e-40, 1e-300);
            assert_rug_close(&(a.0 * b.0), &(a.1.clone() * &b.1), 1e-40, 1e-300);
            assert_rug_close(
                &(a.0 / b.0),
                &(a.1.clone() / &b.1),
                1e-40,
                1e-300 / &b.1.to_f64().abs(),
            );
        }
    }
}
