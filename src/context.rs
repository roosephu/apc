use crate::{
    traits::{ComplexFunctions, MyReal},
    unchecked_cast::UncheckedCast,
};
use log::{debug, info};
use num_complex::Complex;

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
    pub fn two(&self) -> T { 2.0f64.unchecked_cast::<T>() }
}

impl<T: MyReal> Context<T> {
    /// error estimation
    /// See [here](https://www.wikiwand.com/en/Stirling%27s_approximation#Error_bounds) for more details.
    fn loggamma_err(&self, ln_z: Complex<f64>, n: usize) -> f64 {
        let arg = ln_z.im;
        let norm = ln_z.re.exp();
        let err_coef = if arg < std::f64::consts::FRAC_PI_4 {
            1.0
        } else if arg < std::f64::consts::FRAC_PI_2 {
            1.0 / arg.sin().abs()
        } else {
            panic!("you should normalize z first!, z = {:?}", ln_z);
            1.0 / (arg / 2.0).cos().pow(2 * n as i32)
        };
        self.bernoulli(n * 2).unchecked_cast::<f64>().abs()
            / ((2 * n) * (2 * n - 1)) as f64
            / norm.powi((2 * n - 1) as i32)
            * err_coef
    }

    /// log Gamma function by Stirling series
    ///
    pub fn loggamma(&self, mut z: Complex<T>, eps: f64) -> Complex<T> {
        const N: usize = 20;

        assert!(z.re > (-20.0f64).unchecked_cast(), "beyond impl {:?}", z);
        let mut result = Complex::zero();
        while z.re < (N as i64).unchecked_into() {
            result -= z.ln();
            z += T::one();
        }

        let ln_z = z.ln();

        result += (z - 0.5f64.unchecked_cast::<T>()) * ln_z - z + (T::PI() * 2.0).ln() / 2.0;
        let z2 = z * z;
        let mut zpow = z;
        for i in 1..N {
            let contrib = self.bernoulli(i * 2) / ((2 * i) * (2 * i - 1)) as f64 / zpow;
            result += contrib;

            zpow *= z2;
            if contrib.l1_norm().unchecked_cast::<f64>() * 10.0 < eps {
                break;
            }
        }
        let err = self.loggamma_err(ln_z.approx(), N);
        assert!(err <= eps);

        result
    }
}

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
            factorial[i] = factorial[i - 1] * (i as i32).unchecked_cast::<T>();
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
            bernoulli[i * 2] *= self.factorial(2 * i) / 4i32.unchecked_cast::<T>().pow(i as i32);
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
