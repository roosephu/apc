use num::Complex;

use crate::traits::MyReal;

pub trait ComplexFunctions {
    fn mul_i(self) -> Self;
    fn mul_pow_i(self, k: usize) -> Self;
    fn approx(self) -> Complex<f64>;
    fn exp_simul(self) -> Self;
}

impl<T: MyReal> ComplexFunctions for Complex<T> {
    #[inline]
    fn mul_i(self) -> Self { Self { re: -self.im, im: self.re } }

    fn mul_pow_i(self, k: usize) -> Self {
        match k % 4 {
            0 => self,
            1 => Self { re: -self.im, im: self.re },
            2 => Self { re: -self.re, im: -self.im },
            3 => Self { re: self.im, im: -self.re },
            _ => unreachable!(),
        }
    }

    #[inline]
    fn approx(self) -> Complex<f64> {
        Complex::<f64> { re: self.re.to_f64().unwrap(), im: self.im.to_f64().unwrap() }
    }

    #[inline(never)]
    fn exp_simul(self) -> Self {
        let r = self.re.exp();
        let (sin_im, cos_im) = self.im.sin_cos();
        Self { re: r * cos_im, im: r * sin_im }
    }
}
