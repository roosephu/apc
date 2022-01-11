use std::{num::ParseIntError, path::PathBuf};

use clap::{AppSettings, Parser, Subcommand};
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
#[clap(global_setting(AppSettings::DeriveDisplayOrder))]
#[clap(global_setting(AppSettings::PropagateVersion))]
#[clap(global_setting(AppSettings::SubcommandRequiredElseHelp))]
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
        #[clap(name = "t", help = "to which height we compute the zeros of ζ")]
        height: f64,

        #[clap(long, name = "eps", default_value = "1e-18", help = "absolute tolerance")]
        atol: f64,

        #[clap(long, name = "eps", default_value = "1e-18", help = "relative tolerance")]
        rtol: f64,

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

            // TODO: incorporate apcdb and lmfdb. Need to check precision.
            let _ = apcdb;

            let x = parse_int(&x).unwrap();
            info!("======= computing π({}) ======", x);
            // assert!(n >= 100000);

            let mut platt = PlattBuilder::default()
                .hint_λ(lambda_hint)
                .sieve_segment(parse_int(&sieve_segment).unwrap())
                .poly_order(poly_order)
                .ζ_zeros_path(PathBuf::from(lmfdb.unwrap()))
                .build(x);

            let ans = platt.compute::<T>();
            println!("{}", ans);
        }
        Commands::Zero { height: t, atol, rtol, apcdb } => {
            use apc::apcdb::write_APCDB;
            use apc::zeta_zeros::{try_isolate, HybridPrecHardyZ};
            use F64x2::f64x2;

            let mut roots = vec![f64x2::new(14.134725141734695, -8.407109157053214e-16)];

            const PI: f64 = std::f64::consts::PI;
            let mut hardy_z = HybridPrecHardyZ::<f64x2>::new(t, 10, atol);
            let n0 = 0;
            let n1 = t / 2.0 / PI * (t / 2.0 / PI).ln()
                - (0.112 * t.ln() + 0.278 * t.ln().ln() + 3.385 + 0.2 / t);
            let n1 = n1 as usize + 200;
            let stats = try_isolate(&mut hardy_z, n0, n1, atol, rtol);
            assert!(stats.height >= t);
            roots.extend(stats.roots.iter().copied().filter(|&x| x.fp() <= t));

            let n_calls_separate = stats.count_separate;
            let n_calls_locate = stats.count_locate;
            let n_zeros = stats.roots.len();
            info!(
                "{} zeros located. The next call should starts with n0 = {}",
                n_zeros,
                n0 + n_zeros
            );
            info!(
                "{:.3} calls to separate, {:.3} calls to locate, total = {:.3}",
                n_calls_separate as f64 / n_zeros as f64,
                n_calls_locate as f64 / n_zeros as f64,
                (n_calls_separate + n_calls_locate) as f64 / n_zeros as f64,
            );

            let r1 = stats.roots[0];
            let rn = stats.roots[n_zeros - 1];
            info!("the first zero is {}, and the last zero is {}", r1, rn);

            let _ = write_APCDB(&roots, format!("{}/height_100000.dat", apcdb), 0, roots.len());
        }
    }
}
