#![feature(trait_alias)]
#![allow(non_snake_case)]

use analytic::*;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    env_logger::init();

    assert!(args.len() == 2);
    let n = args[1].to_string().parse::<u64>().unwrap();
    println!("======= computing pi({}) ======", n);

    let ctx = Context::<f64>::new(100);
    let mut zeta_galway = ZetaGalway::new(&ctx);
    let mut galway = Galway::new(&ctx, &mut zeta_galway);
    let ans = galway.compute(n);
    println!("[Galway] ans = {}", ans);
    println!("[ZetaGalway] complexity = {}", zeta_galway.complexity);
}
