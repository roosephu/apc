
#![allow(non_snake_case)]

use std::f64::consts::PI;
use analytic::brentq::brentq;
use analytic::gamma::gamma;
use analytic::zeta::zeta;

type Float = f64;
type Complex = num::Complex<Float>;

struct Context {
    prec: u32,
    pi: Float,
    pi2: Float,
    pi_i: Complex,
    pi_i2: Complex,
}

impl Context {
    fn new() -> Context {
        Context {
            prec: 30,
            pi: PI,
            pi2: PI * 2.0,
            pi_i: Complex::new(0.0, PI),
            pi_i2: Complex::new(0.0, PI * 2.0),
        }
    }

    fn exp_pi_i(&self, mut s: Complex) -> Complex {
        s *= &self.pi_i;
        s.exp()
    }
}

use analytic::traits::Alloc;

fn test() {
    let prec = 53;
    // let a = rug::Complex::with_val(prec, (1.0, 1.0));
    let a = rug::Complex::with_val(prec, 1.0);
    let b = a.sin_ref().alloc::<rug::Complex>(prec);
}

fn main() {
    let eps = 1e-10;
    let s = Complex::new(1.11, 100.0);
    let result = zeta(s, eps);
    dbg!(result);
    println!("Hello, world!");
}
