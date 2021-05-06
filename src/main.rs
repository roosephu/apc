#![feature(trait_alias)]
#![allow(non_snake_case)]
#![feature(non_ascii_idents)]
#![feature(try_blocks)]

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
    #[clap(long)]
    poly_order: Option<usize>,
}

fn main() {
    env_logger::init();
    let opts: Opts = Opts::parse();

    let n = opts.n;
    println!("======= computing pi({}) ======", n);
    assert!(n >= 100000);

    // let ctx = Context::<T>::new(100);
    // let mut zeta_galway = ZetaGalway::new(&ctx);
    // let mut galway = Galway::new(&ctx, &mut zeta_galway);
    // let hints = GalwayHints { lambda: opts.lambda };
    // let ans = galway.compute(n, hints);
    // println!("[Galway] ans = {}", ans);
    // println!("[ZetaGalway] complexity = {}", zeta_galway.complexity);

    let mut platt = Platt::<T>::new();
    let hints = PlattHints { λ: opts.lambda, poly_order: opts.poly_order };
    let ans = platt.compute(n, hints);
    println!("[Platt] ans = {}", ans);
}
