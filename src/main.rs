#![feature(trait_alias)]
#![allow(non_snake_case)]

use analytic::{context::Context, f64xn::f64x2, galway::Galway, traits::MyFloat, zeta::ZetaGalway};

fn test<T: MyFloat>(x: T) -> T {
    x + T::one()
}

fn the_default() {
    println!("default implementation");
}

trait Foo {
    fn method(&self) {
        the_default()
    }
}

struct Bar;

impl Foo for Bar {
    fn method(&self) {
        the_default();
        println!("Hey, I'm doing something entirely different!");
    }
}

fn main() {
    let b = Bar;
    b.method();
    let x = f64x2::new(1.0, 2.0);
    let y = test(x);

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
