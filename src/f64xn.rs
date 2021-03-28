#[allow(non_camel_case_types)]
#[derive(Clone, Copy, PartialEq, PartialOrd, Debug, Default)]
pub struct f64x2 {
    pub hi: f64,
    pub lo: f64,
}

mod blocks;
mod functions;
mod impl_traits;
mod remez;
