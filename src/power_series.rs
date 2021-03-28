use std::ops::{AddAssign, DivAssign, MulAssign};

use num::{Float, Num};
use num_traits::NumAssignOps;

use crate::unchecked_cast::UncheckedFrom;

#[derive(Debug)]
pub struct PowerSeries<T, const N: usize> {
    pub data: [T; N],
}

impl<T: Copy + Num, const N: usize> PowerSeries<T, N> {
    fn new(x: T) -> Self {
        let mut ret = Self { data: [T::zero(); N] };
        ret.data[0] = x;
        ret.data[1] = T::one();
        ret
    }
}

impl<T: Copy + NumAssignOps, const N: usize> AddAssign<&PowerSeries<T, N>> for PowerSeries<T, N> {
    fn add_assign(&mut self, rhs: &Self) {
        for i in 0..N {
            self.data[i] += rhs.data[i];
        }
    }
}

impl<T: Copy + Num + NumAssignOps, const N: usize> MulAssign<&PowerSeries<T, N>>
    for PowerSeries<T, N>
{
    fn mul_assign(&mut self, rhs: &Self) {
        for i in (0..N).rev() {
            let mut s = T::zero();
            for j in 0..=i {
                s += self.data[i - j] * rhs.data[j];
            }
            self.data[i] = s;
        }
    }
}

impl<T: Copy + Num + NumAssignOps, const N: usize> DivAssign<&PowerSeries<T, N>>
    for PowerSeries<T, N>
{
    fn div_assign(&mut self, rhs: &Self) {
        for i in 0..N {
            let d = self.data[i] / rhs.data[0];
            for j in 0..N - i {
                self.data[i + j] -= d * rhs.data[j];
            }
            self.data[i] = d;
        }
    }
}

impl<T: Copy + Num + UncheckedFrom<i32> + NumAssignOps, const N: usize> PowerSeries<T, N> {
    /// function composition
    /// Assuming N is small, brute-forcing can typically be very fast.
    /// input = \sum_{i=0}^\infty a_i t^i
    /// let x = a_0, D = \sum_{i=1}^\infty a_i t^i
    /// then f(x + D) = \sum_{i=0}^\infty f^{(i)}(x) D^i / i!
    fn compose(&mut self, mut derivatives: [T; N]) {
        let x = self.data[0];
        let mut series = [T::zero(); N];
        let mut factorial = T::one();
        for i in 0..N {
            derivatives[i] /= factorial;
            factorial *= T::unchecked_from(i as i32 + 1);
        }

        series[0] = self.data[0];
        self.data[0] = T::zero();
        for i in (0..N).rev() {
            // basically it's `ret = ret * D + f^{(i)}`
            // We truncate the series to length n.
            let n = N - i;
            for j in (0..n).rev() {
                let mut s = T::zero();
                for k in 1..=j {
                    s += series[j - k] * self.data[k];
                }
                series[j] = s;
            }
            series[0] += derivatives[i];
        }
        self.data = series;
    }
}

impl<T: Copy + Float + UncheckedFrom<i32> + NumAssignOps, const N: usize> PowerSeries<T, N> {
    fn cos_(&mut self) {
        let mut derivatives = [T::zero(); N];
        let (sin_x, cos_x) = self.data[0].sin_cos();
        for i in 0..N {
            derivatives[i] = match i % 4 {
                0 => cos_x,
                1 => -sin_x,
                2 => -cos_x,
                3 => sin_x,
                _ => unreachable!(),
            }
        }
        self.compose(derivatives);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn composition() {
        type PS = PowerSeries<f64, 10>;
        let z = 0.2;
        let mut numer = PS::new(z * z * PI / 2.0 + 3.0 * PI / 8.0);
        numer.data[1] = PI * z;
        numer.data[2] = PI / 2.0;
        numer.cos_();

        let mut denom = PS::new(z * PI);
        denom.data[1] = PI;
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
