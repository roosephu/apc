#![allow(uncommon_codepoints)]
#![allow(non_snake_case)]

use F64x2::f64x2;
use apc::traits::MyReal;
use apc::bandwidth_interp::BandwidthInterp;
use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn bench_bandwidth_interp(c: &mut Criterion) {
    type T = f64x2;

    let k = black_box(100);
    let sigma = black_box(T::mp(0.5));
    let eps = black_box(1e-26);
    let min_t = black_box(T::mp(2.0 * f64::PI() * (k as f64).powi(2)));
    let max_t = black_box(T::mp(2.0 * f64::PI() * (k as f64 + 1.0).powi(2)));
    let ds = BandwidthInterp::<T>::new(k, min_t, max_t, sigma, eps);

    c.bench_function("bandwidth interpolation", |b| {
        b.iter(|| {
            let t = black_box(T::PI() * 2.0 * 10100.0);
            ds.query(t)
        })
    });
}

criterion_group!(benches, bench_bandwidth_interp);
criterion_main!(benches);
