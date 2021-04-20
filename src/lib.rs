#![allow(non_snake_case)]
#![allow(dead_code)]
#![allow(unused_variables)]
#![allow(unused_imports)]
#![allow(clippy::float_cmp)]
#![allow(clippy::needless_range_loop)]
#![allow(clippy::let_and_return)]
#![feature(trait_alias)]
#![feature(destructuring_assignment)]
#![feature(min_type_alias_impl_trait)]

pub mod bandwidth_interp;
pub mod brentq;
pub mod constants;
pub mod context;
pub mod f64xn;
pub mod galway;
pub mod gamma;
pub mod power_series;
pub mod riemann_siegel;
pub mod sum_trunc_dirichlet;
pub(crate) mod test_utils;
pub mod traits;
pub mod unchecked_cast;
pub mod zeta;

pub use context::Context;
pub use f64xn::f64x2;
pub use galway::{Galway, GalwayHints};
pub use riemann_siegel::{RiemannSiegelZ, RiemannSiegelZeta};
pub use zeta::ZetaGalway;
