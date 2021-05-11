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

cfg_if::cfg_if! {
    if #[cfg(feature = "zeta")] {
        mod bandwidth_interp;
        mod gamma;
        mod riemann_siegel;
        mod sum_trunc_dirichlet;
        mod zeta;
        pub use riemann_siegel::{RiemannSiegelZ, RiemannSiegelZeta};
        pub use zeta::ZetaGalway;
    }
}

cfg_if::cfg_if! {
    if #[cfg(feature = "galway")] {
        mod galway;
        pub use galway::{Galway, GalwayHints};
    }
}

mod brentq;
mod cache_stat;
mod constants;
mod context;
mod fast_phi;
mod lmfdb;
mod platt;
mod platt_integral;
mod power_series;
mod sieve;
mod traits;

pub use context::Context;
pub use fast_phi::LittlePhiSum;
pub use platt::PlattBuilder;

