use std::path::PathBuf;

use apc::*;
use clap::{clap_app, crate_authors, crate_description, crate_version};
use log::info;

type T = F64x2::f64x2;

fn main() {
    let opts = clap_app!(apc =>
        (@arg x: * "Compute π(x)")
        (@arg poly_order: --("poly-order") default_value("15") "Degree of polynomial approximation used in Φ̂(s)")
        (@arg lambda_hint: --("lambda-hint") default_value("10") "Hint of λ")
        (@arg zeta_zeros_path: --("zeta-zeros-path")  default_value("./data/zeros") "Path to LMFDB ζ zeros")
        // (@arg log_level: --"log-level" +takes_value "Show debugging information")
    )
    .version(crate_version!())
    .author(crate_authors!())
    .about(crate_description!())
    .get_matches();

    env_logger::builder().format_timestamp_micros().init();

    let x = parse_int::parse::<u64>(opts.value_of("x").unwrap()).unwrap();
    info!("======= computing pi({}) ======", x);
    // assert!(n >= 100000);

    // let ctx = Context::<T>::new(100);
    // let mut zeta_galway = ZetaGalway::new(&ctx);
    // let mut galway = Galway::new(&ctx, &mut zeta_galway);
    // let hints = GalwayHints { lambda: opts.lambda };
    // let ans = galway.compute(n, hints);
    // println!("[Galway] ans = {}", ans);

    let mut builder = PlattBuilder::default();
    opts.value_of("lambda_hint").map(|v| builder.hint_λ(v.parse::<f64>().unwrap()));
    opts.value_of("poly_order").map(|v| builder.poly_order(v.parse::<usize>().unwrap()));
    opts.value_of("zeta_zero_path").map(|v| builder.ζ_zeros_path(PathBuf::from(v)));
    let mut platt = builder.build(x);

    let ans = platt.compute::<T>();
    println!("{}", ans);
}
