use apc::*;
use clap::{crate_authors, crate_version, Clap};
use F64x2::f64x2;

type T = f64x2;

#[derive(Clap, Debug)]
#[clap(version = crate_version!(), author = crate_authors!())]
struct Opts {
    n: u64,
    #[clap(long)]
    lambda_hint: Option<f64>,
    #[clap(long)]
    poly_order: Option<usize>,
}

fn main() {
    env_logger::builder().format_timestamp_micros().init();
    let opts: Opts = Opts::parse();

    let n = opts.n;
    println!("======= computing pi({}) ======", n);
    // assert!(n >= 100000);

    // let ctx = Context::<T>::new(100);
    // let mut zeta_galway = ZetaGalway::new(&ctx);
    // let mut galway = Galway::new(&ctx, &mut zeta_galway);
    // let hints = GalwayHints { lambda: opts.lambda };
    // let ans = galway.compute(n, hints);
    // println!("[Galway] ans = {}", ans);
    // println!("[ZetaGalway] complexity = {}", zeta_galway.complexity);

    let hints = PlattHints { λ: opts.lambda_hint, poly_order: opts.poly_order };
    let mut platt = Platt::new(n, hints);
    let ans = platt.compute::<T>(n);
    println!("[Platt] ans = {}", ans);
}
