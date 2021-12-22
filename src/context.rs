use crate::traits::MyReal;
use log::{debug, info};

const N: usize = 256;

pub struct Context<T> {
    binomial: Vec<Vec<T>>,
    bernoulli: Vec<T>,
    euler: Vec<T>,
    factorial: Vec<T>,
    PI: T,
    pow_pi: [T; N],
    pow2: [T; N],
}

impl<T: MyReal> Context<T> {
    pub fn new(n: usize) -> Self {
        let mut ret = Self {
            binomial: vec![],
            bernoulli: vec![],
            euler: vec![],
            factorial: vec![],
            PI: T::PI(),
            pow_pi: [T::zero(); N],
            pow2: [T::zero(); N],
        };
        ret.init_factorial(n * 2 + 1);
        ret.init_binomial(n * 2);
        ret.init_bernoulli(n);
        ret.init_euler(n);
        ret.init_pows();
        ret
    }

    #[inline]
    pub fn binom(&self, n: usize, m: usize) -> T { self.binomial[n][m] }
    #[inline]
    pub fn bernoulli(&self, n: usize) -> T { self.bernoulli[n] }
    #[inline]
    pub fn euler(&self, n: usize) -> T { self.euler[n] }
    #[inline]
    pub fn factorial(&self, n: usize) -> T { self.factorial[n] }
    #[inline]
    pub fn pow_pi(&self, n: usize) -> T { self.pow_pi[n] }
    #[inline]
    pub fn pow2(&self, n: usize) -> T { self.pow2[n] }

    #[inline]
    pub fn two(&self) -> T { T::mp(2.0) }
}

impl<T: MyReal> Context<T> {}

impl<T: MyReal> Context<T> {
    fn init_binomial(&mut self, n: usize) {
        info!("initialize binomial numbers up to {}", n);
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

    fn init_factorial(&mut self, n: usize) {
        info!("initialize factorial up to {}", n);
        let mut factorial = vec![T::one(); n + 1];
        for i in 1..=n {
            factorial[i] = factorial[i - 1] * T::mp(i as f64);
        }
        self.factorial = factorial;
    }

    /// initialize Bernoulli numbers
    /// Use a recurrence from [here](https://maths-people.anu.edu.au/~brent/pd/tangent.pdf) for a more stable recurrence.
    fn init_bernoulli(&mut self, n: usize) {
        info!("initialize Bernoulli numbers up to {}", n);
        let mut bernoulli = vec![T::zero(); n + 1];
        bernoulli[0] = T::one();
        bernoulli[1] = T::one() / 2.0;

        for i in 1..=n / 2 {
            let mut b = T::zero();
            for j in 0..i {
                b += self.factorial(2 * i) / self.factorial(2 * i + 1 - 2 * j) * bernoulli[2 * j];
            }
            bernoulli[2 * i] = (T::one() - b) / self.factorial(2 * i);
        }
        for i in 1..=n / 2 {
            bernoulli[i * 2] *= self.factorial(2 * i) / T::mp(4.0).pow(i as i32);
        }

        debug!("bernoulli numbers = {:?}", &bernoulli[..10]);

        self.bernoulli = bernoulli;
    }

    fn init_euler(&mut self, n: usize) {
        info!("initialize Euler numbers up to {}", n);
        let mut euler = vec![T::zero(); n + 1];
        euler[0] = T::one();
        for i in 1..=n {
            let mut s = T::zero();
            for j in 0..i {
                s += (if (i + j) % 2 == 0 { T::one() } else { -T::one() })
                    * euler[j]
                    * self.binom(2 * i, 2 * j);
            }
            euler[i] = -s;
        }
        self.euler = euler;
    }

    fn init_pows(&mut self) {
        self.pow_pi[0] = T::one();
        self.pow2[0] = T::one();
        for i in 1..N {
            self.pow_pi[i] = self.pow_pi[i - 1] * T::PI();
            self.pow2[i] = self.pow2[i - 1] * 2.0;
        }
    }
}
