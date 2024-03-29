#![allow(non_snake_case)]
#![allow(dead_code)]
#![allow(unused_variables)]
#![allow(clippy::float_cmp)]
#![allow(clippy::needless_range_loop)]
#![allow(clippy::let_and_return)]
#![allow(clippy::approx_constant)]
#![feature(type_alias_impl_trait)]
#![feature(const_float_bits_conv)]
#![feature(trait_alias)]

#[allow(non_camel_case_types)]
#[derive(Clone, Copy, PartialEq, PartialOrd, Debug, Default)]
pub struct f64x2 {
    pub hi: f64,
    pub lo: f64,
}

pub mod blocks;
mod f64xn;
mod functions;
mod impl_traits;
mod remez;
pub mod traits;

pub mod test_utils;
