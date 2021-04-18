use std::{
    cmp::{max, min},
    ops::{AddAssign, DivAssign, MulAssign, SubAssign},
};

use num::{Float, Num};
use num_traits::NumAssignOps;

use crate::unchecked_cast::UncheckedFrom;

#[derive(Debug, Clone)]
pub struct PowerSeries<T> {
    pub n: usize,
    pub N: usize,
    pub data: Vec<T>,
}

impl<T: Copy + Num> PowerSeries<T> {
    pub fn new(max_order: usize, x: T) -> Self {
        let mut data = vec![T::zero(); max_order];
        data[0] = x;
        data[1] = T::one();
        Self { n: 2, N: max_order, data }
    }

    pub fn from_vec(max_order: usize, mut data: Vec<T>) -> Self {
        data.resize(max_order, T::zero());
        Self { n: data.len(), N: max_order, data }
    }
}

impl<T: Copy + NumAssignOps> AddAssign<&PowerSeries<T>> for PowerSeries<T> {
    fn add_assign(&mut self, rhs: &Self) {
        self.n = max(self.n, rhs.n);
        for i in 0..self.n {
            self.data[i] += rhs.data[i];
        }
    }
}

impl<T: Copy + NumAssignOps> SubAssign<&PowerSeries<T>> for PowerSeries<T> {
    fn sub_assign(&mut self, rhs: &Self) {
        self.n = max(self.n, rhs.n);
        for i in 0..self.n {
            self.data[i] -= rhs.data[i];
        }
    }
}

impl<T: Copy + Num + NumAssignOps> MulAssign<T> for PowerSeries<T> {
    fn mul_assign(&mut self, rhs: T) {
        for i in 0..self.n {
            self.data[i] *= rhs;
        }
    }
}

impl<T: Copy + Num + NumAssignOps> MulAssign<&PowerSeries<T>> for PowerSeries<T> {
    fn mul_assign(&mut self, rhs: &Self) {
        let mut data = vec![T::zero(); self.N];
        let new_n = min(self.n + rhs.n, self.N);
        for i in 0..self.n {
            for j in 0..=min(new_n - i, rhs.n) {
                data[i + j] += self.data[i] * rhs.data[j];
            }
        }
        self.n = new_n;
        self.data.clone_from_slice(&data);
    }
}

impl<T: Copy + Num + NumAssignOps> DivAssign<&PowerSeries<T>> for PowerSeries<T> {
    fn div_assign(&mut self, rhs: &Self) {
        let mut data = self.data.clone();
        let N = self.N;
        for i in 0..N {
            let d = data[i] / rhs.data[0];
            self.data[i] = d;
            for j in 0..N - i {
                data[i + j] -= d * rhs.data[j];
            }
        }
    }
}

impl<T: Copy + Num + UncheckedFrom<i32> + NumAssignOps> PowerSeries<T> {
    /// function composition
    /// Assuming N is small, brute-forcing can typically be very fast.
    /// input = \sum_{i=0}^\infty a_i t^i
    /// let x = a_0, D = \sum_{i=1}^\infty a_i t^i
    /// then f(x + D) = \sum_{i=0}^\infty f^{(i)}(x) D^i / i!
    fn compose(&mut self, derivatives: &mut [T]) {
        let x = self.data[0];
        let mut series = vec![T::zero(); self.N];
        let mut factorial = T::one();
        for i in 0..self.N {
            derivatives[i] /= factorial;
            factorial *= T::unchecked_from(i as i32 + 1);
        }

        series[0] = self.data[0];
        self.data[0] = T::zero();
        for i in (0..self.N).rev() {
            // basically it's `ret = ret * D + f^{(i)}`
            // We truncate the series to order n.
            let n = self.N - i;
            for j in (0..n).rev() {
                let mut s = T::zero();
                for k in 1..min(j + 1, self.n) {
                    s += series[j - k] * self.data[k];
                }
                series[j] = s;
            }
            series[0] += derivatives[i];
        }
        self.n = self.N;
        self.data = series;
    }
}

impl<T: Copy + Float + UncheckedFrom<i32> + NumAssignOps> PowerSeries<T> {
    pub fn cos_(&mut self) {
        let mut derivatives = vec![T::zero(); self.N];
        let (sin_x, cos_x) = self.data[0].sin_cos();
        let table = [cos_x, -sin_x, -cos_x, sin_x];
        for i in 0..self.N {
            derivatives[i] = table[i % 4];
        }
        self.compose(&mut derivatives);
    }

    pub fn sin_(&mut self) {
        let mut derivatives = vec![T::zero(); self.N];
        let (sin_x, cos_x) = self.data[0].sin_cos();
        let table = [sin_x, cos_x, -sin_x, -cos_x];
        for i in 0..self.N {
            derivatives[i] = table[i % 4];
        }
        self.compose(&mut derivatives);
    }

    pub fn exp_(&mut self) {
        let exp = self.data[0].exp();
        let mut derivatives = vec![exp; self.N];
        self.compose(&mut derivatives);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn composition() {
        type PS = PowerSeries<f64>;
        let z = 0.2;
        let mut numer = PS::from_vec(20, vec![z * z * PI / 2.0 + 3.0 * PI / 8.0, PI * z, PI / 2.0]);
        numer.cos_();

        let mut denom = PS::from_vec(20, vec![z * PI, PI]);
        denom.cos_();

        // println!("numer  = {:?}", numer);
        // println!("denom  = {:?}", denom);

        numer /= &denom;
        // println!("result = {:?}", numer);

        let gt = [
            0.40038394798914,
            0.17910471843231912,
            0.4686598273731518,
            0.10347886444127957,
            0.1226725478937987,
            -0.022532170247690487,
            -0.02932815687717693,
            -0.023188194967356527,
            -0.016239785896345466,
            -0.0026737163058292426,
        ];
        for i in 0..10 {
            assert!((gt[i] - numer.data[i]).abs() < 1e-14);
        }
    }
}