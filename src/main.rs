#![feature(trait_alias)]
#![allow(non_snake_case)]

use analytic::*;

type T = f64x2;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    env_logger::init();

    assert!(args.len() == 2);
    let n = args[1].to_string().parse::<u64>().unwrap();
    println!("======= computing pi({}) ======", n);

    let ctx = Context::<T>::new(100);
    let mut zeta_galway = ZetaGalway::new(&ctx);
    let mut galway = Galway::new(&ctx, &mut zeta_galway);
    let ans = galway.compute(n);
    println!("[Galway] ans = {}", ans);
    println!("[ZetaGalway] complexity = {}", zeta_galway.complexity);
}

// fn main() {
//     let x = f64x2 { hi: -4959337563492667.0, lo: -0.2204808216865093 };
//     let y = x.floor();
//     println!("x = {:?}, y = {:?}", x, y);
// }
