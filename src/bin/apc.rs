#![feature(explicit_generic_args_with_impl_trait)]

use std::{num::ParseIntError, path::PathBuf};

use apc::{LMFDB, APCDB};
use clap::{AppSettings, Parser, Subcommand};
use log::info;

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
#[clap(global_setting(AppSettings::DeriveDisplayOrder))]
#[clap(global_setting(AppSettings::PropagateVersion))]
#[clap(global_setting(AppSettings::UseLongFormatForHelpSubcommand))]
// #[clap(global_setting(AppSettings::SubcommandRequiredElseHelp))]
struct App {
    #[clap(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
#[allow(non_snake_case)]
enum Commands {
    Pi {
        #[clap(name = "x", help = "Compute π(x)")]
        x: String,

        #[clap(
            long,
            name = "k",
            default_value_t = 15,
            help = "Degree of polynomial approximation used in Φ̂(s)"
        )]
        poly_order: usize,

        #[clap(long, name = "lambda", default_value_t = 10.0, help = "Hint of λ")]
        lambda_hint: f64,

        #[clap(long, name = "seg", default_value = "2^32", help = "Segment size in sieve")]
        sieve_segment: String,

        #[clap(long, name = "dir", help = "Directory to LMFDB ζ zeros. Override `APCDB`.")]
        lmfdb: Option<String>,

        #[clap(
            long,
            name = "dir",
            help = "Directory to APCDB ζ zeros. If both `LMFDB` and `APCDB` are unset, generate zeros online."
        )]
        apcdb: Option<String>,
    },
    Zero {
        #[clap(name = "t", help = "To which height we compute the zeros of Riemann ζ function")]
        height: f64,

        #[clap(long, name = "eps", default_value = "1e-18", help = "Absolute tolerance")]
        atol: f64,

        #[clap(long, name = "eps", default_value = "1e-18", help = "Relative tolerance")]
        rtol: f64,

        #[clap(long, name = "k", default_value_t = 10, help = "Max Riemann-Siegel formula order")]
        max_rs_order: usize,

        #[clap(
            long,
            name = "dir",
            default_value = "./data/apcdb",
            help = "Directory to save the zeros"
        )]
        apcdb: String,
    },
}

fn main() {
    apc::init();
    let args = App::parse();

    match args.command {
        Commands::Pi { x, poly_order, lambda_hint, sieve_segment, lmfdb, apcdb } => {
            use apc::PlattBuilder;
            use F64x2::f64x2;

            let x = parse_int(&x).unwrap();
            info!("======= computing π({}) ======", x);
            // assert!(n >= 100000);

            let mut platt = PlattBuilder::default()
                .hint_λ(lambda_hint)
                .sieve_segment(parse_int(&sieve_segment).unwrap())
                .poly_order(poly_order)
                .build(x);

            let ans = match (lmfdb, apcdb) {
                (Some(p), _) => platt.compute::<f64x2>(LMFDB::<f64x2>::directory(p)),
                (None, Some(p)) => platt.compute::<f64x2>(APCDB::<f64x2>::directory(p)),
                (None, None) => unreachable!(),  // TODO: generate it online
            };

            println!("{}", ans);
        }
        Commands::Zero { height: t, atol, rtol, max_rs_order, apcdb } => {
            use apc::write_APCDB;
            use apc::zeta_zeros::{try_isolate, HybridPrecHardyZ};
            use F64x2::f64x2;

            let mut hardy_z = HybridPrecHardyZ::<f64x2>::new(t, max_rs_order, atol);
            let stats = try_isolate(&mut hardy_z, -1, t, atol, rtol, true);
            assert!(stats.gram_end.1 >= t);
            let roots = stats.roots.iter().copied().filter(|&x| x.fp() <= t).collect::<Vec<_>>();

            let n_calls_separate = stats.count_separate;
            let n_calls_locate = stats.count_locate;
            let n_zeros = stats.roots.len();
            info!("Seek range: [0, {t:.3}]");
            info!(
                "Certified : [g({}), g({}))] = [{:.3}, {:.3}]. {} zeros located. ",
                stats.gram_start.0, stats.gram_end.0, stats.gram_start.1, stats.gram_end.1, n_zeros
            );
            info!(
                "{:.3} calls to separate, {:.3} calls to locate, total = {:.3}",
                n_calls_separate as f64 / n_zeros as f64,
                n_calls_locate as f64 / n_zeros as f64,
                (n_calls_separate + n_calls_locate) as f64 / n_zeros as f64,
            );

            let r1 = roots[0];
            let rn = roots[n_zeros - 1];
            info!("the first zero is {}, and the last zero is {}", r1, rn);

            let file_path = format!("{}/height_100000.dat", apcdb);
            let _ = write_APCDB(&roots, file_path, 0, roots.len(), atol, rtol);
        }
    }
}
