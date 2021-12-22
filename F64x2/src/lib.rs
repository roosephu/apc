#![allow(non_snake_case)]
#![allow(dead_code)]
#![allow(unused_variables)]
#![allow(clippy::float_cmp)]
#![allow(clippy::needless_range_loop)]
#![allow(clippy::let_and_return)]
#![allow(clippy::approx_constant)]
#![feature(type_alias_impl_trait)]
#![feature(destructuring_assignment)]
#![feature(const_float_bits_conv)]

#[allow(non_camel_case_types)]
#[derive(Clone, Copy, PartialEq, PartialOrd, Debug, Default)]
pub struct f64x2 {
    pub(crate) hi: f64,
    pub(crate) lo: f64,
}

mod blocks;
mod f64xn;
mod functions;
mod impl_traits;
mod remez;
pub mod test_utils;
