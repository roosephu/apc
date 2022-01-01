use std::ops::{Div, Neg};

use f64xn_impl_arith::f64xn_impl_arith;
mod multi_floats;

use crate::blocks::*;

#[allow(non_camel_case_types)]
#[derive(Copy, Clone, PartialEq, PartialOrd, Debug)]
struct f64x<const N: usize> {
    pub(crate) data: [f64; N],
}

impl<const N: usize> f64x<N> {
    const ONE: Self = Self::mp(1.0);
    const ZERO: Self = Self::mp(0.0);

    #[inline]
    const fn mp(x: f64) -> Self {
        let mut data = [0.0; N];
        data[0] = x;
        Self { data }
    }

    #[inline]
    const fn fp(&self) -> f64 { self.data[0] }
}

impl<const N: usize> Neg for f64x<N> {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self::Output { Self { data: self.data.map(|x| -x) } }
}

impl<const N: usize> Div<f64x<N>> for f64
where
    f64x<N>: Div<f64x<N>, Output = f64x<N>>,
{
    type Output = f64x<N>;

    #[inline]
    fn div(self, rhs: f64x<N>) -> Self::Output { f64x::<N>::mp(self) / rhs }
}

macro_rules! f64xn_assoc_std_ops {
    ($ty: ty, $trait: ident, $method: ident, $inner_method: ident) => {
        impl std::ops::$trait for $ty {
            type Output = Self;
            #[inline]
            fn $method(self, rhs: Self) -> Self {
                Self { data: Self::$method(self.data, rhs.data) }
            }
        }

        impl std::ops::$trait<f64> for $ty {
            type Output = Self;
            #[inline]
            fn $method(self, rhs: f64) -> Self {
                Self { data: Self::$inner_method(self.data, rhs) }
            }
        }
    };
}

macro_rules! f64xn_define {
    ($ty: ident, $n: expr) => {
        f64xn_impl_arith!($n);

        #[allow(non_camel_case_types)]
        type $ty = f64x<$n>;
        f64xn_assoc_std_ops!($ty, Add, add, add_f64);
        f64xn_assoc_std_ops!($ty, Sub, sub, sub_f64);
        f64xn_assoc_std_ops!($ty, Mul, mul, mul_f64);
        f64xn_assoc_std_ops!($ty, Div, div, div_f64);
    };
}

f64xn_define!(f64x4, 4);

// f64xn_impl_arith!(3);

// include!(concat!(env!("OUT_DIR"), "/impl_f64xn.rs"));
