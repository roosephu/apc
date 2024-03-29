#![allow(uncommon_codepoints)]
#![allow(non_snake_case)]

use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn bench_fast_phi(c: &mut Criterion) {
    let x: u64 = black_box(1000000000000000);
    let λ = black_box(5.728162e-7);
    let eps = black_box(0.1);
    let x1: u64 = black_box(999997603247000);
    let x2: u64 = black_box(1000002396758745);

    c.bench_function("fast ϕ", |b| {
        b.iter(|| {
            let mut n_primes_less_than_x = 0usize;
            let mut ϕ = apc::LittlePhiSum::new(λ, x as f64, eps / (x2 - x1) as f64);

            let mut calc = |p: u64| {
                let t = (p - x) as i64;
                if t <= 0 {
                    n_primes_less_than_x += 1;
                }
                ϕ.add(t)
            };

            let step = (x as f64).ln() as usize;
            for p in (x1..x2).step_by(step) {
                calc(p);
            }
            n_primes_less_than_x as f64 - ϕ.sum()
        })
    });
}

criterion_group!(benches, bench_fast_phi);
criterion_main!(benches);
