#![feature(trait_alias)]
#![allow(non_snake_case)]
#![feature(non_ascii_idents)]
#![feature(try_blocks)]

use crate::platt::PlattHints;
use apc::*;
use clap::{crate_authors, crate_version, Clap};
use F64x2::f64x2;

type T = f64x2;

#[derive(Clap, Debug)]
#[clap(version = crate_version!(), author = crate_authors!())]
struct Opts {
    n: u64,
    #[clap(long)]
    lambda: Option<f64>,
}

fn main() {
    env_logger::init();
    let opts: Opts = Opts::parse();

    let n = opts.n;
    println!("======= computing pi({}) ======", n);

    // let ctx = Context::<T>::new(100);
    // let mut zeta_galway = ZetaGalway::new(&ctx);
    // let mut galway = Galway::new(&ctx, &mut zeta_galway);
    // let hints = GalwayHints { lambda: opts.lambda };
    // let ans = galway.compute(n, hints);
    // println!("[Galway] ans = {}", ans);
    // println!("[ZetaGalway] complexity = {}", zeta_galway.complexity);

    let mut platt = crate::platt::Platt::<T>::new();
    let hints = PlattHints { λ: opts.lambda };
    let ans = platt.compute(n, hints);
    println!("[Platt] ans = {}", ans);

    unsafe {
        println!(
            "[erfc] large = [count = {}, K = {}], small = [count = {}, K = {}]",
            crate::profiler::ERFC_EPS_LARGE_CNT,
            crate::profiler::ERFC_EPS_LARGE_K,
            crate::profiler::ERFC_EPS_SMALL_CNT,
            crate::profiler::ERFC_EPS_SMALL_K,
        );
    }
}

#[cfg(test)]
mod tests {
    use crate::unchecked_cast::*;
    use crate::zeta::FnZeta;
    use crate::*;
    use num::Complex;
    use F64x2::f64x2;

    #[test]
    fn test_tmp() {
        let ctx1 = Context::<f64>::new(100);
        let mut zeta_galway1 = ZetaGalway::new(&ctx1);
        zeta_galway1.prepare_multi_eval(0.23939822958279525, 1e-10);

        let ctx2 = Context::<f64x2>::new(100);
        let mut zeta_galway2 = ZetaGalway::new(&ctx2);
        zeta_galway2.prepare_multi_eval(0.23939822958279525.unchecked_into(), 1e-10);

        let s1 = Complex::<f64>::new(1.5, 0.23939822958279525);
        let s2 = Complex::<f64x2>::new(1.5.unchecked_into(), 0.23939822958279525.unchecked_into());
        let result1 = zeta_galway1.zeta(s1, 1e-10);
        let result2 = zeta_galway2.zeta(s2, 1e-10);
        println!("{:?}", result1);
        println!("{:?}", result2);
    }
}
