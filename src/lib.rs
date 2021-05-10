#![allow(non_snake_case)]
#![allow(dead_code)]
#![allow(confusable_idents)]
#![allow(uncommon_codepoints)]
#![allow(clippy::float_cmp)]
#![allow(clippy::needless_range_loop)]
#![allow(clippy::let_and_return)]
#![allow(clippy::approx_constant)]
#![feature(trait_alias)]
#![feature(destructuring_assignment)]
#![feature(min_type_alias_impl_trait)]

mod bandwidth_interp;
mod brentq;
mod cache_stat;
mod constants;
mod context;
mod fast_phi;
mod galway;
mod gamma;
mod lmfdb;
mod platt;
mod platt_integral;
mod power_series;
mod riemann_siegel;
mod sieve;
mod sum_trunc_dirichlet;
mod traits;
mod zeta;

pub use context::Context;
pub use fast_phi::LittlePhiSum;
pub use galway::{Galway, GalwayHints};
pub use platt::PlattBuilder;
pub use riemann_siegel::{RiemannSiegelZ, RiemannSiegelZeta};
pub use zeta::ZetaGalway;
