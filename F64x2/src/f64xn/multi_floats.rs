use crate::f64xn::f64x;
use num::{One, Zero};
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

// impl<const N: usize> From<f64x<N>> for String {
//     fn from(mut a: f64x::<N>) -> Self {
//         // TODO
//         if a.is_nan() {
//             String::from("NaN")
//         } else if a.is_infinite() {
//             if a.is_sign_negative() {
//                 String::from("-Inf")
//             } else {
//                 String::from("Inf")
//             }
//         } else if a.is_zero() {
//             String::from("0")
//         } else {
//             let mut ret = String::from("");
//             if a.is_sign_negative() {
//                 a = -a;
//                 ret.push('-');
//             }
//             let mut e = 0;
//             while a >= f64x::<N>::from(10.0) {
//                 a = a / 10.;  // TODO: avoid division
//                 e += 1;
//             }
//             while a * 10. < f64x::<N>::from(1.) {
//                 a = a * 10.;
//                 e -= 1;
//             }

//             let mut dec_point = false;

//             for _ in 0..30 {
//                 let d = a.floor() as u8;
//                 ret.push((b'0' + d) as char);
//                 a = (a - f64x::<N>::from(d as f64)) * 10.;
//                 if a.is_zero() {
//                     break;
//                 }

//                 if !dec_point {
//                     ret.push('.');
//                     dec_point = true;
//                 }
//             }
//             if e != 0 {
//                 ret.push('E');
//                 ret.push_str(e.to_string().as_str());
//             }

//             ret
//         }
//     }
// }

// extern crate proc_macro;
