#![allow(non_snake_case)]
#![allow(dead_code)]
#![allow(clippy::float_cmp)]
#![allow(clippy::many_single_char_names)]
#![allow(clippy::needless_range_loop)]
#![allow(unused_variables)]
#![feature(trait_alias)]
#![feature(destructuring_assignment)]
#![feature(min_type_alias_impl_trait)]

pub mod brentq;
pub mod context;
pub mod f64xn;
pub mod galway;
pub mod gamma;
pub mod power_series;
pub mod riemann_siegel;
pub mod sum_trunc_dirichlet;
pub mod traits;
pub mod unchecked_cast;
pub mod zeta;

pub use context::Context;
pub use f64xn::f64x2;
pub use galway::{Galway, GalwayHints};
pub use zeta::ZetaGalway;
