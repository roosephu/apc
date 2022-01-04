//! https://arxiv.org/abs/1410.7176
use rug::Float;

fn repr(x: &Float) -> String {
    let hi = x.to_f64();
    let lo = (x.clone() - hi).to_f64();
    String::from(format!("f64x2 {{ hi: {:?}, lo: {:?} }}", hi, lo))
}

/// Outputs exp(ln(2) / 2 / 2^k * i) for every i in [-2^k, 2^k]
/// To comptue exp(x) for |x| < ln(2) / 2, we can always find a $x'$
/// such that $|x'| < ln(2) / 2 / 2^(k+1)$ and $x - x' = i / 2^k$ for some i. .
fn gen_exp_table(k: usize, prec: u32) {
    assert!(k <= 20);
    let limit = 1i32 << k;
    let delta = 2.0f64.ln() / 2.0f64.powi(k as i32 + 1);
    eprintln!("delta = {:?}", delta);
    println!("[");
    for i in -limit..=limit {
        let s = Float::with_val(prec, delta * i as f64).exp();
        println!("{},", repr(&s));
    }
    println!("]");
}

fn main() {
    gen_exp_table(14, 200);
}
