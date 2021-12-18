use std::{num::ParseIntError, path::PathBuf};

use apc::*;
use clap::{clap_app, crate_authors, crate_description, crate_version};
use log::info;

type T = F64x2::f64x2;

#[allow(clippy::single_char_pattern)]
fn parse_int(s: &str) -> Result<u64, ParseIntError> {
    let s = String::from(s);
    let s = s.replace("_", "");
    if s.contains("^") {
        let v: Vec<_> = s.split("^").collect();
        let a = v[0].parse::<u64>()?;
        let b = v[1].parse::<u32>()?;
        Ok(a.pow(b))
    } else {
        s.parse::<u64>()
    }
}

fn main() {
    let opts = clap_app!(apc =>
        (@arg x: * "Compute π(x)")
        (@arg poly_order: --("poly-order") default_value("15") "Degree of polynomial approximation used in Φ̂(s)")
        (@arg lambda_hint: --("lambda-hint") default_value("10") "Hint of λ")
        (@arg sieve_segment: --("sieve-segment") default_value("2^32") "Segment size in sieve")
        (@arg zeta_zeros_path: --("zeta-zeros-path")  default_value("./data/zeros") "Path to LMFDB ζ zeros")
        (@arg verbose: --verbose "Shows debugging information")
    )
    .version(crate_version!())
    .author(crate_authors!())
    .about(crate_description!())
    .get_matches();

    let mut builder = env_logger::builder();
    builder.format_timestamp_micros();
    if opts.is_present("verbose") {
        builder.filter_level(log::LevelFilter::Info);
    }
    builder.init();

    let x = parse_int(opts.value_of("x").unwrap()).unwrap();
    info!("======= computing π({}) ======", x);
    // assert!(n >= 100000);

    let mut builder = PlattBuilder::default();
    opts.value_of("lambda_hint").map(|v| builder.hint_λ(v.parse::<f64>().unwrap()));
    opts.value_of("sieve_segment").map(|v| builder.sieve_segment(parse_int(v).unwrap()));
    opts.value_of("poly_order").map(|v| builder.poly_order(v.parse::<usize>().unwrap()));
    opts.value_of("zeta_zero_path").map(|v| builder.ζ_zeros_path(PathBuf::from(v)));
    let mut platt = builder.build(x);

    let ans = platt.compute::<T>();
    println!("{}", ans);
}
