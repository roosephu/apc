use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};

use rand;
use F64x2::blocks::{two_add, two_add_fast, two_add_if};

fn bench_two_add(c: &mut Criterion) {
    let mut group = c.benchmark_group("two add");

    let setup = || {
        let a = f64::from_bits(rand::random::<u64>());
        let b = f64::from_bits(rand::random::<u64>());
        (a, b)
    };

    group.bench_function("two add", |b| {
        b.iter_batched(setup, |(a, b)| black_box(two_add(a, b)), criterion::BatchSize::SmallInput)
    });
    group.bench_function("two add if", |b| {
        b.iter_batched(
            setup,
            |(a, b)| black_box(two_add_if(a, b)),
            criterion::BatchSize::SmallInput,
        )
    });
    group.bench_function("two_add_fast", |b| {
        b.iter_batched(
            setup,
            |(a, b)| black_box(two_add_fast(a, b)),
            criterion::BatchSize::SmallInput,
        )
    });
    group.finish();
}

criterion_group!(benches, bench_two_add);
criterion_main!(benches);
