#![allow(dead_code)]
use apc::zeta_zeros::{try_isolate, HybridPrecHardyZ};
use log::info;
use F64x2::f64x2;

fn main() {
    apc::init();

    let t = 1.5e5;
    const PI: f64 = std::f64::consts::PI;
    let mut hardy_z = HybridPrecHardyZ::<f64x2>::new(t, 10, 1e-18);
    let n1 = t / 2.0 / PI * (t / 2.0 / PI).ln()
        - (0.112 * t.ln() + 0.278 * t.ln().ln() + 3.385 + 0.2 / t);
    let n1 = n1 as usize;
    let (roots, stats) = try_isolate(&mut hardy_z, 3, n1, 1e-18, 1e-30);
    for i in 0..10 {
        println!("{}", roots[i]);
    }
    let n_calls_separate = stats.count[0];
    let n_calls_locate = stats.count[1];
    let n_zeros = roots.len();
    info!(
        "To separate {} zeros: {:.3} calls to separate, {:.3} calls to locate, total = {:.3}",
        n_zeros,
        n_calls_separate as f64 / n_zeros as f64,
        n_calls_locate as f64 / n_zeros as f64,
        (n_calls_separate + n_calls_locate) as f64 / n_zeros as f64,
    );

    let r1 = roots[0];
    let rn = roots[n_zeros - 1];
    info!("the first zero is {}, and the last zero is {}", r1, rn);
}
