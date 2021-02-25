use std::ops::MulAssign;

use rug::{float::SmallFloat, ops::*};
use rug::Assign;

type Float = rug::Float;
type Complex = rug::Complex;

struct Context {
    prec: u32,
    fbuf: Vec<Float>,
    cbuf: Vec<Complex>,
}

impl Context {
    fn f<T>(&mut self, val: T) -> &mut Float where Float: Assign<T> {
        &mut Float::with_val(self.prec, val)
    }

    fn c<T>(&mut self, val: T) -> &mut Complex where Complex: Assign<T> {
        &mut Complex::with_val(self.prec, val)
    }

    fn pi(&self) -> &Float {
        &Float::with_val(self.prec, rug::float::Constant::Pi)
    }

    fn enter(&mut self) {

    }

    fn leave(&mut self) {

    }
}

// https://epubs.siam.org/doi/pdf/10.1137/0731050
pub fn gamma2(ctx: &mut Context, z: &Complex, eps: Float) -> Complex {
    let PI: &Float = ctx.pi();
    let PI2: Float = PI * 2;
    let z: &Complex = ctx.c(z - &1.0);
    let mut a: Float = z.real().clone().max(ctx.f(2.0)).floor() + 0.5;
    loop {
        let err = a.sqrt() / (PI2).pow(a + 0.5) / ctx.f(z.real() + &a);
        if err < eps { break }
        a += 1.0;
    }
    let mut coef = ctx.c(1);
    let mut k = 1;
    let mut fac = ctx.f(1);
    while k < a {
        let diff: &Float = ctx.f(&a - k);
        let t1 = ctx.f(&a);
        let c_k = ctx.f(diff.pow(k as f64 - 0.5)) * diff.exp() / PI2.sqrt() / fac;
        fac.mul_assign(k);
        if k % 2 == 0 {
            coef += c_k / ctx.c(&z + k);
        } else {
            coef -= c_k / ctx.c(&z + k);
        }
        k += 1;
    }
    // dbg!(coef, (z + a).powc(z + 0.5) / (z + a).exp() * (2.0 * PI).sqrt());
    coef * ctx.c(&z + a).pow(&z + 0.5) / ctx.c(&z + &a).exp() * PI2.sqrt()
}

#[cfg(test)]
mod tests {
    use super::{gamma2, Complex, Float};

    fn _test_gamma(z: Complex, eps: Float, gt: Complex) {
        let result = gamma2(z, eps);
        let diff = (result - gt).norm();
        println!("{:?} {:?}", result, gt);
        assert!(diff <= eps);
    }
    #[test]
    fn test_gamma() {
        _test_gamma(Complex::new(4.0, 10.0), 1e-10, Complex::new(0.0007715342942399662, -0.0010190827990417));
        _test_gamma(Complex::new(4.0, 0.0), 1e-10, Complex::new(6.0, 0.0));
    }
}
