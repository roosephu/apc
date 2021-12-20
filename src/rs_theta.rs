use crate::{traits::MyReal, contexts::Bernoulli};


pub struct RiemannSiegelTheta<T> {
    K: usize,
    coeffs: Vec<T>,
}

impl<T: MyReal + Bernoulli> RiemannSiegelTheta<T> {
    /// see https://arxiv.org/pdf/1609.03682.pdf for the "wrong" formula
    /// also see https://arblib.org/gamma.html
    pub fn new(K: usize) -> Self {
        let mut coeffs = vec![T::zero(); K + 1];
        for j in 1..=K {
            coeffs[j] = (T::one() - T::from_f64(2.0).unwrap().powi(1 - 2 * j as i32))
                * T::bernoulli(2 * j).abs()
                / (4 * j * (2 * j - 1)) as f64;
        }
        Self { K, coeffs }
    }

    // See [Sec 3.11, Pugh].
    pub fn theta(&self, t: T, eps: f64) -> T {
        // as it's typically used with RiemannSiegelZ, we hope it's not too small.
        assert!(t.to_f64().unwrap() >= 200.0 && eps > 1e-33);
        const K: usize = 7;

        // needs high precision base computation here.
        let mut ret = t / 2.0 * (t / 2.0 / T::PI() / T::E()).ln() - T::FRAC_PI_8();
        let mut tpow = t;
        let tsqr = t * t;
        for i in 1..=K {
            ret += self.coeffs[i] / tpow;
            tpow *= tsqr;
        }
        ret
    }
}
