#![allow(dead_code)]
use apc::zeta_zeros::{try_isolate, HybridPrecHardyZ};
use log::info;
use F64x2::f64x2;

fn main() {
    apc::init();

    // let mut hardy_z = HardyZ::<f64x2>::new(1e8, 10, 1e-18);
    let mut hardy_z = HybridPrecHardyZ::<f64x2>::new(1e8, 10, 1e-18);
    let (roots, stats) = try_isolate(&mut hardy_z, 100002, 300000, 1e-18, 1e-30);
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
