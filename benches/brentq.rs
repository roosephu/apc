#![allow(uncommon_codepoints)]
#![allow(non_snake_case)]

use apc::brentq::*;
use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn bench_brentq(c: &mut Criterion) {
    let xtol = black_box(1e-13);
    let rtol = black_box(0.0);
    let n_iters = black_box(100);

    c.bench_function("find_root-Brentq", |b| {
        b.iter(|| {
            let mut ans = 0.0;
            for i in 1..=100 {
                let f = black_box(|x: f64| x * x.exp() - i as f64);
                let a = EvalPoint::new(black_box(0.0), f);
                let b = EvalPoint::new(black_box((i as f64).ln() * 2.0 + 1.0), f);

                let result1 = Brentq::new(a, b, xtol, rtol).solve(f, n_iters);
                ans += result1.x;
            }
            ans
        })
    });

    c.bench_function("find_root-brentq2", |b| {
        b.iter(|| {
            let mut ans = 0.0;
            for i in 1..=100 {
                let f = black_box(|x: f64| x * x.exp() - i as f64);
                let a = EvalPoint::new(black_box(0.0), f);
                let b = EvalPoint::new(black_box((i as f64).ln() * 2.0 + 1.0), f);

                let result2 = brentq2(f, a, b, xtol, rtol, n_iters);
                ans += result2.x;
            }
            ans
        })
    });
}

criterion_group!(benches, bench_brentq);
criterion_main!(benches);
