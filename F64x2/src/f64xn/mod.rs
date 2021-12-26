use std::ops::{Add, Div, Mul, Neg, Sub};

// use multi_floats_impl::impl_f64x;
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

    const fn mp(x: f64) -> Self {
        let mut data = [0.0; N];
        data[0] = x;
        Self { data }
    }

    const fn fp(&self) -> f64 { self.data[0] }
}

impl<const N: usize> Neg for f64x<N> {
    type Output = Self;
    fn neg(self) -> Self::Output { Self { data: self.data.map(|x| -x) } }
}

impl<const N: usize> Div<f64x<N>> for f64
where
    f64x<N>: Div<f64x<N>, Output = f64x<N>>,
{
    type Output = f64x<N>;

    fn div(self, rhs: f64x<N>) -> Self::Output { f64x::<N>::mp(self) / rhs }
}

// impl_f64x!(1);
// impl_f64x!(2);
// impl_f64x!(3);
// impl_f64x!(4);
// impl_f64x!(5);
// impl_f64x!(6);
// impl_f64x!(7);
// impl_f64x!(8);

include!(concat!(env!("OUT_DIR"), "/impl_f64xn.rs"));
