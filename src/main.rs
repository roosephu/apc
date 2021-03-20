
#![allow(non_snake_case)]

use analytic::{context::Context, galway::Galway, zeta::{FnZeta, ZetaGalway}};
use num::Complex;


fn main() {
    env_logger::init();

    let ctx = Context::<f64>::new(100);
    let zeta_galway = ZetaGalway::new(&ctx);

    // let s = Complex::new(1.5, 21846.405929161672);
    // let lambda = 0.000200;
    // let z = zeta_galway.zeta(s, 0.0000000000000000021714724095162593);
    // let ln_x = (1e10f64).ln();
    // let Psi = ((lambda * s).powi(2) / 2.0 + s * ln_x).exp() * z.ln() / s;
    // println!("Psi = {}", Psi);

    let mut galway = Galway::new(&ctx, &zeta_galway);
    let ans = galway.compute(10_000_000_000);
    println!("[Galway] ans = {}", ans);
    println!("[ZetaGalway] complexity = {}", zeta_galway.complexity);
}
