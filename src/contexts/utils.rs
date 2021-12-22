use log::info;
pub use paste;
use rug::{ops::CompleteRound, Float};
use std::{
    io::{BufRead, BufReader},
    path::Path,
};
use F64x2::f64x2;

#[macro_export]
macro_rules! impl_from_uninit_cell {
    ($trait: ident, $method: ident, $ty: ty) => {
        $crate::contexts::utils::paste::paste! {
            #[allow(non_upper_case_globals)]
            static [<TABLE_ $ty>]: $crate::contexts::lazy_static::UninitCell<Vec<$ty>> = $crate::contexts::lazy_static::UninitCell::uninit();

            impl $trait for $ty {
                #[inline]
                fn $method(n: usize) -> Self {
                    [<TABLE_ $ty>][n]
                }
            }
        }
    };
}

pub use impl_from_uninit_cell;

pub fn read_data(path: impl AsRef<Path>, prec: u32) -> Result<Vec<Float>, std::io::Error> {
    let file = std::fs::File::open(&path)?;
    let reader = BufReader::new(file);
    let ret: Vec<Float> =
        reader.lines().map(|x| Float::parse(x.unwrap()).unwrap().complete(prec)).collect();
    info!("read {} numbers from {}", ret.len(), path.as_ref().display());
    Ok(ret)
}

pub fn mpf_to_f64x2(x: &Float) -> f64x2 {
    let hi = x.to_f64();
    let lo = (x - hi).complete(x.prec()).to_f64();
    f64x2::new(hi, lo)
}
