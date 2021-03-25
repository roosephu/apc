use log::debug;
use num::Complex;
use rustfft::{FftNum, FftPlanner};

use crate::traits::MyFloat;

pub trait ExpPolyApprox: Sized {
    type Output: Iterator<Item = (usize, Complex<Self>)>;

    fn get_poly_approx() -> Self::Output;
}

const COEFFS: [Complex<f64>; 18] = [
    Complex { re: 1.0, im: 1.0812144107799266e-23 },
    Complex { re: -4.479610786345225e-21, im: 1.0 },
    Complex { re: -0.5, im: 3.027415987199093e-19 },
    Complex { re: -8.413828845781633e-18, im: -0.16666666666666669 },
    Complex { re: 0.04166666666666679, im: 1.1656663448809618e-16 },
    Complex { re: -1.0603511152022404e-15, im: 0.008333333333332324 },
    Complex { re: -0.0013888888888827402, im: 5.789264486508273e-15 },
    Complex { re: -2.491995923872859e-14, im: -0.00019841269843586228 },
    Complex { re: 2.4801587374768556e-5, im: 6.704175576866034e-14 },
    Complex { re: -1.594987515102099e-13, im: 2.7557317787217356e-6 },
    Complex { re: -2.755729303110001e-7, im: 2.3127460502103687e-13 },
    Complex { re: -3.2663668749921504e-13, im: -2.5052389834713885e-8 },
    Complex { re: 2.087985316554709e-9, im: 2.5867211760028217e-13 },
    Complex { re: -2.2167241850689593e-13, im: 1.6041263496425594e-10 },
    Complex { re: -1.1352710114429515e-11, im: 8.943908448871146e-14 },
    Complex { re: -4.542339711641447e-14, im: -7.962911435347713e-13 },
    Complex { re: 5.979573239083729e-14, im: 7.185782517642856e-15 },
    Complex { re: -1.970149077208406e-15, im: 1.9701490772084063e-15 },
];

impl ExpPolyApprox for f64 {
    type Output = impl Iterator<Item = (usize, Complex<Self>)>;

    fn get_poly_approx() -> Self::Output { COEFFS.iter().enumerate().map(|(idx, &x)| (idx, x)) }
}

// compute \sum_{j=1}^n j^{-(s + i j delta)} for every 0 <= t < m
pub fn sum_trunc_dirichlet<T: ExpPolyApprox + MyFloat + FftNum>(
    s: Complex<T>,
    n: usize,
    m: usize,
    delta: T,
) -> Vec<Complex<T>> {
    debug!("[OS-FKBJ] s = {:.6}, n = {}, m = {}, delta = {:.6}", s, n, m, delta);
    let M2 = (m + m % 2) / 2;
    let TM2 = T::from(M2).unwrap();
    let s = s + Complex::new(T::zero(), TM2 * delta);
    let a: Vec<_> =
        (1..=n).map(|x| Complex::new(T::from(x).unwrap(), T::zero()).powc(-s)).collect();
    let g: Vec<T> = (1..=n).map(|x| T::from(x).unwrap().ln() * -delta).collect();
    let R = (m + 1).next_power_of_two() * 2;
    let div = T::TAU() / T::from(R).unwrap();
    let w: Vec<i64> = g.iter().map(|&x| ((x / div).round().as_())).collect();
    let d: Vec<T> = g.iter().zip(w.iter()).map(|(&x, &y)| x - T::from(y).unwrap() * div).collect();

    let mut ret = vec![Complex::zero(); m + 1];
    let mut f = vec![Complex::zero(); R];
    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_inverse(f.len());

    for (e, c) in T::get_poly_approx() {
        for x in f.iter_mut() {
            x.set_zero()
        }
        for j in 0..n {
            f[w[j].rem_euclid(R as i64) as usize] += a[j] * (d[j] * TM2).pow(e as i32);
        }
        fft.process(&mut f);
        for x in 0..=m {
            ret[x] +=
                c * f[(x + R - M2) % R] * (T::from(x).unwrap() / TM2 - T::one()).pow(e as i32);
        }
    }

    ret
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sum_trunc_dirichlet() {
        let N = 20;
        let s = Complex::new(1.3, 2.6);
        let M = 9;
        let delta = 1.1;
        let result = sum_trunc_dirichlet(s, N, M, delta);
        for i in 0..=M {
            let mut sum = Complex::<f64>::zero();
            let z = s + Complex::new(0.0, i as f64 * delta);
            for j in 1..=N {
                sum += Complex::new(j as f64, 0.0).powc(-z);
            }
            assert!((sum - result[i]).norm() < 2e-15);
        }
    }
}
