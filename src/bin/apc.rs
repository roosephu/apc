use std::{num::ParseIntError, path::PathBuf};

use apc::*;
use clap::Parser;
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

#[derive(Parser, Debug)]
#[clap(about, version, author)]
struct Args {
    #[clap(short, long, help = "Compute π(x)")]
    x: String,

    #[clap(
        short,
        long,
        default_value_t = 15,
        help = "Degree of polynomial approximation used in Φ̂(s)"
    )]
    poly_order: usize,

    #[clap(short, long, default_value_t = 10.0, help = "Hint of λ")]
    lambda_hint: f64,

    #[clap(short, long, default_value = "2^32", help = "Segment size in sieve")]
    sieve_segment: String,

    #[clap(short, long, default_value = "./data/zeros", help = "Path to LMFDB ζ zeros")]
    zeta_zeros_path: String,

    #[clap(short, long, help = "Shows debugging information")]
    verbose: bool,
}

fn main() {
    let opts = Args::parse();

    let mut builder = env_logger::builder();
    builder.format_timestamp_micros();
    if opts.verbose {
        builder.filter_level(log::LevelFilter::Info);
    }
    builder.init();

    let x = parse_int(&opts.x).unwrap();
    info!("======= computing π({}) ======", x);
    // assert!(n >= 100000);

    let mut platt = PlattBuilder::default()
        .hint_λ(opts.lambda_hint)
        .sieve_segment(parse_int(&opts.sieve_segment).unwrap())
        .poly_order(opts.poly_order)
        .ζ_zeros_path(PathBuf::from(opts.zeta_zeros_path))
        .build(x);

    let ans = platt.compute::<T>();
    println!("{}", ans);
}
