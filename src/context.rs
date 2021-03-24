use crate::traits::{ComplexFunctions, GenericFloat};
use log::{debug, info};
use num_complex::Complex;
use num_traits::AsPrimitive;

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
        ret.init_factorial(n * 2 + 1);
        ret.init_binomial(n * 2);
        ret.init_bernoulli(n);
        ret.init_euler(n);
        ret
    }

    #[inline]
    pub fn binom(&self, n: usize, m: usize) -> T {
        self.binomial[n][m]
    }
    #[inline]
    pub fn bernoulli(&self, n: usize) -> T {
        self.bernoulli[n]
    }
    #[inline]
    pub fn euler(&self, n: usize) -> T {
        self.euler[n]
    }
    #[inline]
    pub fn factorial(&self, n: usize) -> T {
        self.factorial[n]
    }
}

impl<T: GenericFloat> Context<T> {
    /// error estimation
    /// See [here](https://www.wikiwand.com/en/Stirling%27s_approximation#Error_bounds) for more details.
    fn loggamma_err(&self, ln_z: Complex<f64>, n: usize) -> f64 {
        let arg = ln_z.im;
        let norm = ln_z.re.exp();
        let err_coef = if arg < std::f64::consts::PI / 4.0 {
            1.0
        } else if arg < std::f64::consts::PI / 2.0 {
            1.0 / arg.sin().abs()
        } else {
            1.0 / (arg / 2.0).cos().pow(2 * n as i32)
        };
        AsPrimitive::<f64>::as_(self.bernoulli(n * 2)) / ((2 * n) * (2 * n - 1)) as f64 * err_coef
    }

    /// log Gamma function by Stirling series
    ///
    pub fn loggamma(&self, mut z: Complex<T>, eps: f64) -> Complex<T> {
        const N: usize = 20;

        let mut result = Complex::zero();
        while z.re < T::from(N).unwrap() {
            result -= z.ln();
            z += T::one();
        }

        let ln_z = z.ln();

        result += (z - T::from(0.5).unwrap()) * ln_z - z
            + (T::PI() * T::from(2).unwrap()).ln() / T::from(2).unwrap();
        let z2 = z * z;
        let mut zpow = z;
        for i in 1..N {
            result += self.bernoulli(i * 2) / T::from((2 * i) * (2 * i - 1)).unwrap() / zpow;
            zpow *= z2;
        }
        let err = self.loggamma_err(ln_z.approx(), N);
        assert!(err <= eps);

        result
    }
}

impl<T: GenericFloat> Context<T> {
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
            factorial[i] = factorial[i - 1] * T::from(i).unwrap();
        }
        self.factorial = factorial;
    }

    /// initialize Bernoulli numbers
    /// Use a recurrence from [here](https://maths-people.anu.edu.au/~brent/pd/tangent.pdf) for a more stable recurrence.
    fn init_bernoulli(&mut self, n: usize) {
        info!("initialize Bernoulli numbers up to {}", n);
        let mut bernoulli = vec![T::zero(); n + 1];
        bernoulli[0] = T::one();
        bernoulli[1] = -T::from(2).unwrap().recip();

        for i in 1..=n / 2 {
            let mut b = T::zero();
            for j in 0..i {
                b += self.factorial(2 * i) / self.factorial(2 * i + 1 - 2 * j) * bernoulli[2 * j];
            }
            bernoulli[2 * i] = (T::one() - b) / self.factorial(2 * i);
        }
        for i in 1..=n / 2 {
            bernoulli[i * 2] *= self.factorial(2 * i) / T::from(4).unwrap().pow(i as i32);
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
                s += (if (i + j) % 2 == 0 {
                    T::one()
                } else {
                    -T::one()
                }) * euler[j]
                    * self.binom(2 * i, 2 * j);
            }
            euler[i] = -s;
        }
        self.euler = euler;
    }
}