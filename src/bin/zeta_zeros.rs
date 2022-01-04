#![allow(dead_code)]
use apc::zeta_zeros::{try_isolate, RiemannSiegelZ};
use log::info;
use F64x2::f64x2;

fn main() {
    apc::init();

    let mut rsz = RiemannSiegelZ::<f64x2>::new(1e8, 1e-12);

    let roots = try_isolate(&mut rsz, 100002, 300000, 1e-5, 1e-30);
    let n_calls_separate = rsz.counts[0];
    let n_calls_locate = rsz.counts[1];
    let n_zeros = roots.len();
    info!(
        "To separate {} zeros: {:.3} calls to separate, {:.3} calls to locate",
        n_zeros,
        n_calls_separate as f64 / n_zeros as f64,
        n_calls_locate as f64 / n_zeros as f64
    );

    let r1 = roots[0];
    let rn = roots[n_zeros - 1];
    info!("the first zero is {}, and the last zero is {}", r1, rn);
}
