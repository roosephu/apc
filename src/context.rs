use crate::traits::GenericFloat;
use num_complex::Complex;
use log::{info, debug};

#[derive(Default)]
pub struct Context<T> {
    binomial: Vec<Vec<T>>,
    bernoulli: Vec<T>,
    euler: Vec<T>,
    factorial: Vec<T>,
}

impl<T: GenericFloat> Context<T> {
    pub fn new(n: usize) -> Self {
        let mut ret = Self::default();
        ret.init_binomial(n * 2);
        ret.init_bernoulli(n);
        ret.init_euler(n);
        ret
    }

    #[inline] pub fn binom(&self, n: usize, m: usize) -> T { self.binomial[n][m] }
    #[inline] pub fn bernoulli(&self, n: usize) -> T { self.bernoulli[n] }
    #[inline] pub fn euler(&self, n: usize) -> T { self.euler[n] }
    #[inline] pub fn factorial(&self, n: usize) -> T { self.factorial[n] }
}

impl<T: GenericFloat> Context<T> {
    pub fn loggamma(&self, mut z: Complex<T>, eps: f64) -> Complex<T> {
        let mut result = Complex::zero();
        while z.re < T::from(20).unwrap() {
            result -= z.ln();
            z += T::one();
        }

        result += (z - T::from(0.5).unwrap()) * z.ln() - z + (T::PI() * T::from(2).unwrap()).ln() / T::from(2).unwrap();
        let z2 = z * z;
        let mut zpow = z;
        for i in 1..20 {
            result += self.bernoulli(i * 2) / T::from((2 * i) * (2 * i - 1)).unwrap() / zpow;
            zpow *= z2;
        }
        result
    }
}

impl<T: GenericFloat> Context<T> {
    fn init_binomial(&mut self, n: usize) {
        let mut binomial = vec![vec![T::zero(); n + 1]; n + 1];
        for i in 0..=n {
            for j in 0..=i {
                if j == 0 {
                    binomial[i][j] = T::one();
                } else {
                    binomial[i][j] = binomial[i - 1][j - 1] + binomial[i - 1][j];
                }
            }
        }
        self.binomial = binomial
    }

    fn init_bernoulli(&mut self, n: usize) {
        info!("initialize Bernoulli numbers up to {}", n);
        let mut bernoulli = vec![T::zero(); n + 1];
        bernoulli[0] = T::from(1.0).unwrap();

        for i in 1..=n {
            let mut b = T::zero();
            for j in 0..i {
                b += self.binom(i, j) * bernoulli[j] / T::from(i - j + 1).unwrap();
            }
            bernoulli[i] = -b;
        }

        debug!("bernoulli numbers = {:?}", &bernoulli[..10]);

        self.bernoulli = bernoulli;
    }

    fn init_euler(&mut self, n: usize) {
        let mut euler = vec![T::zero(); n + 1];
        euler[0] = T::one();
        for i in 1..=n {
            let mut s = T::zero();
            for j in 0..i {
                s += (if (i + j) % 2 == 0 { T::one() } else { -T::one()}) * euler[j] * self.binom(2 * i, 2 * j);
            }
            euler[i] = -s;
        }
        self.euler = euler;
    }
}
