use criterion::{black_box, criterion_group, criterion_main, Criterion};

use F64x2::f64x2;

fn bench_exp(c: &mut Criterion) {
    let s = black_box(f64x2::new(3.4538776394910684, 0.0000000000000001184757763427252));

    c.bench_function("exp", |b| {
        b.iter(|| s.exp());
    });
}

criterion_group!(benches, bench_exp);
criterion_main!(benches);
