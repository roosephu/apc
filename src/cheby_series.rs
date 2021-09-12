use std::ops::{Add, AddAssign, DivAssign, Mul, MulAssign, SubAssign};

use num::{FromPrimitive, ToPrimitive};
use num_traits::NumAssign;

use crate::traits::MyReal;

/// n = working order
/// n = max order
#[derive(Debug, Clone)]
pub struct ChebySeries<T> {
    pub n: usize,
    pub coeffs: Vec<T>,
}

impl<T: Copy + Clone + NumAssign> ChebySeries<T> {
    pub fn new<X: MyReal, F: Fn(X) -> T>(n: usize, f: F) -> Self
    where
        T: Mul<X, Output = T> + Add<T, Output = T>,
    {
        let cos: Vec<X> = (0..=2 * n).map(|k| (X::PI() * (k as f64) / (n as f64)).cos()).collect();
        let values: Vec<T> = cos[0..=n].iter().map(|&x| f(x)).collect();

        let mut coeffs = vec![T::zero(); n + 1];
        for a in 0..=n {
            let mut coeff = T::zero();
            for b in 0..=n {
                let c = values[b] * cos[a * b % (2 * n)];
                if b == 0 || b == n {
                    coeff += c;
                } else {
                    coeff += c * X::from_f64(2.0).unwrap();
                }
            }
            coeff = coeff * X::from_usize(n).unwrap().recip();
            coeffs[a] = coeff;
        }

        Self { n, coeffs }
    }

    pub fn eval<X: MyReal>(&self, x: X) -> T
    where
        T: Mul<X, Output = T> + Add<T, Output = T>,
    {
        let mut ret = self.coeffs[0] * X::from_f64(0.5).unwrap();
        let mut cur = x; // T_k(x)
        let mut pre = X::one(); // T_{k-1}(x)
        let x2 = x * 2.0;
        for k in 1..=self.n {
            ret += self.coeffs[k] * cur;

            let new = cur * x2 - pre;
            pre = cur;
            cur = new;
        }
        ret
    }
}

impl<T: NumAssign + Copy> ChebySeries<T> {
    pub fn differentiate_<X: FromPrimitive>(&mut self)
    where
        T: Mul<X, Output = T>,
    {
        let n = self.n;

        let mut coeffs = vec![T::zero(); n + 1];
        for k in (0..n).rev() {
            coeffs[k] = self.coeffs[k + 1] * X::from_usize(2 * (k + 1)).unwrap();
            if k + 2 < n {
                coeffs[k] = coeffs[k + 2] + coeffs[k];
            }
        }

        self.coeffs = coeffs;
    }
}

impl<T: Copy + NumAssign> AddAssign<&ChebySeries<T>> for ChebySeries<T> {
    fn add_assign(&mut self, rhs: &Self) {
        assert!(self.n == rhs.n);
        for i in 0..=self.n {
            self.coeffs[i] += rhs.coeffs[i];
        }
    }
}

impl<T: Copy + NumAssign> SubAssign<&ChebySeries<T>> for ChebySeries<T> {
    fn sub_assign(&mut self, rhs: &Self) {
        assert!(self.n == rhs.n);
        for i in 0..=self.n {
            self.coeffs[i] -= rhs.coeffs[i];
        }
    }
}

impl<X: Copy, T: MulAssign<X>> MulAssign<X> for ChebySeries<T> {
    fn mul_assign(&mut self, rhs: X) {
        for i in 0..=self.n {
            self.coeffs[i] *= rhs;
        }
    }
}

impl<X: Copy, T: DivAssign<X>> DivAssign<X> for ChebySeries<T> {
    fn div_assign(&mut self, rhs: X) {
        for i in 0..=self.n {
            self.coeffs[i] /= rhs;
        }
    }
}

impl<T: Copy + NumAssign> ChebySeries<T> {
    pub fn zero(n: usize) -> Self { Self { n, coeffs: vec![T::zero(); n + 1] } }
}

impl<T> ChebySeries<T> {
    pub fn truncate_(&mut self, n: usize) {
        assert!(n <= self.n);
        self.n = n;
        self.coeffs.truncate(n + 1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cheby() {
        // let f = |x: f64| x.exp() / (2.0 + x.cos());
        // let f_prime = |x: f64| x.exp() * (2.0 + x.cos() + x.sin()) / (2.0 + x.cos()).powi(2);
        // let f = |x: f64| x * x / (2.0 + x.cos());
        // let f_prime = |x: f64| x * (x * x.sin() + 4.0 + 2.0 * x.cos()) / (2.0 + x.cos()).powi(2);
        let f = |x: f64| (-x * x + x).exp() / (2.0 + x);
        let f_prime =
            |x: f64| (-2.0 * x * x - 3.0 * x + 1.0) * (-x * x + x).exp() / (2.0 + x).powi(2);

        // let f = |&x: &f64| { x.exp() };
        // let f_prime = |&x: &f64| { x.exp() };

        // let f = |&x: &f64| { (2.0 * x).exp() };
        // let f_prime = |&x: &f64| { 2.0 * (2.0 * x).exp() };

        let mut cheby_series = ChebySeries::new(30, f);
        let x = 0.3;

        // test new and query
        let gt = f(x);
        let result = cheby_series.eval(x);
        println!("gt = {}, result = {}, diff = {:.6e}", gt, result, gt - result);
        assert!((result - gt).abs() < 5e-16);

        cheby_series.differentiate_::<f64>();
        let gt = f_prime(x);
        let result = cheby_series.eval(x);
        println!("gt = {}, result = {}, diff = {:.6e}", gt, result, gt - result);
        assert!((result - gt).abs() < 5e-15);
    }
}
