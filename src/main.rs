
#![allow(non_snake_case)]

use analytic::{context::Context, galway::Galway, zeta::ZetaGalway};


fn main() {
    env_logger::init();

    let ctx = Context::<f64>::new(100);
    let zeta_galway = ZetaGalway::new(&ctx);
    let mut galway = Galway::new(&ctx, &zeta_galway);
    let ans = galway.compute(10_000_000);
    println!("[Galway] ans = {}", ans);
}
