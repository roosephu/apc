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

    let ctx = Context::<f64>::new(100);
    let mut zeta_galway = ZetaGalway::new(&ctx);
    let mut galway = Galway::new(&ctx, &mut zeta_galway);
    let ans = galway.compute(n);
    println!("[Galway] ans = {}", ans);
    println!("[ZetaGalway] complexity = {}", zeta_galway.complexity);
}


#[cfg(test)]
mod tests {
    use analytic::f64x2;
    use num::{Complex, One};

    #[test]
    fn test() {
        // let x = f64x2 { hi: 3.0, lo: -0.2204808216865093e-4 };
        let s = Complex { re: f64x2 { hi: 1.5, lo: 0.0 }, im: f64x2 { hi: 0.1364376340480039, lo: 0.0 } };
        println!("s = {:?}, 1 - s = {:?}", s, f64x2::one() - s);
    }
}
