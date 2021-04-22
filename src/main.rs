#![feature(trait_alias)]
#![allow(non_snake_case)]

use analytic::*;
use clap::{crate_authors, crate_version, Clap};
use platt::PlattHints;

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
    let hints = PlattHints { lambda: opts.lambda };
    let ans = platt.compute(n, hints);
    println!("[Platt] ans = {}", ans);
}

#[cfg(test)]
mod tests {
    use analytic::unchecked_cast::*;
    use analytic::zeta::FnZeta;
    use analytic::*;
    use num::Complex;

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
