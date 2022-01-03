use log::debug;
use num::Complex;
use rustfft::{FftNum, FftPlanner};

use crate::{contexts::ExpPolyApprox, traits::MyReal};

/// compute $F(t) = \sum_{(a, g)} a \exp(i t g)$ for $t$ in $[-m, m]$
pub(crate) fn sum_weighted_exp<T: MyReal + FftNum>(
    a: &[Complex<T>],
    g: &[T],
    m: usize,
    poly_approx: impl Iterator<Item = (usize, Complex<T>)>,
) -> Vec<Complex<T>> {
    let TM = T::from_usize(m).unwrap();
    let n = a.len();
    let R = (m as usize + 1).next_power_of_two() * 2;
    let div = T::TAU() / (R as f64);
    let w: Vec<i64> = g.iter().map(|&x| ((x / div).round()).to_i64().unwrap()).collect();
    let d: Vec<T> =
        g.iter().zip(w.iter()).map(|(&x, &y)| x - T::from_i64(y).unwrap() * div).collect();

    let mut ret = vec![Complex::zero(); m * 2 + 1];
    let mut f = vec![Complex::zero(); R];
    let mut planner = FftPlanner::<T>::new();
    let fft = planner.plan_fft_inverse(R);

    for (e, c) in poly_approx {
        for x in f.iter_mut() {
            x.set_zero()
        }
        for j in 0..n {
            f[w[j].rem_euclid(R as i64) as usize] += a[j] * (d[j] * TM).pow(e as i32);
        }
        fft.process(f.as_mut_slice());
        for x in 0..=m * 2 {
            ret[x] += c * f[(x + R - m) % R] * (T::from_usize(x).unwrap() / TM - 1.0).pow(e as i32);
        }
    }

    ret
}

/// compute \sum_{j=n0}^n1 j^{-(s + i t delta)} for every 0 <= t <= m
pub(crate) fn sum_trunc_dirichlet<T: MyReal + FftNum + ExpPolyApprox>(
    s: Complex<T>,
    n0: usize,
    n1: usize,
    m: usize,
    delta: T,
) -> Vec<Complex<T>> {
    // debug!("[OS-FKBJ] s = {:.6}, n = [{}, {}], m = {}, delta = {:.6}", s, n0, n1, m, delta);
    let M2 = (m + m % 2) / 2;
    let s = s + Complex::new(T::zero(), delta * M2 as f64);

    // ans[t] = \sum_{j=0}^{n - 1} a[j] * exp(i gamma[j] (t - M2))
    let ln_x: Vec<_> = (n0..=n1).map(|x| T::from_usize(x).unwrap().ln()).collect();
    let a: Vec<_> = ln_x.iter().map(|&ln_x| (-s * ln_x).exp()).collect();
    let gamma: Vec<T> = ln_x.iter().map(|&ln_x| -delta * ln_x).collect();

    sum_weighted_exp(&a, &gamma, M2, T::get_poly_approx())
}

#[cfg(test)]
mod tests {
    use F64x2::f64x2;
    use F64x2::test_utils::*;

    use super::*;

    type T = f64x2;

    #[test]
    fn test_sum_trunc_dirichlet() {
        let N = 20;
        let s = Complex::new(T::mp(1.3), T::mp(260.0));
        let M = 9;
        let delta = T::mp(1.1);
        let result = sum_trunc_dirichlet(s, 1, N, M, delta);
        for t in 0..=M {
            let mut sum = Complex::<T>::zero();
            let z = s + Complex::new(T::zero(), delta * t as f64);
            for j in 1..=N {
                sum += (-z * T::mp(j as f64).ln()).exp();
            }
            println!("{:e} {:e}, diff = {:e}", sum, result[t], sum - result[t]);
            assert_complex_close(sum, result[t], 1e-29);
        }
    }

    /// can test: add, power for complex numbers.
    #[test]
    fn test_truncated_dirichlet_series() {
        let n = 1000;
        let s = Complex::new(f64x2::new(1.5, 0.0), f64x2::new(1.0, 0.0));
        let mut ans: Complex<f64x2> = Complex::zero();
        for i in 1..=n {
            let x = Complex::new(f64x2::new(i as f64, 0.0), f64x2::zero()).powc(s);
            ans += x;
        }
        println!("ans = {}", ans);
        let gt = Complex::new(
            f64x2::new(1.1409175704943191e7, 8.194110941294255e-10),
            f64x2::new(2.8472593984260815e6, 1.4500054775776998e-10),
        );
        assert_complex_close(ans, gt, 6e-32);
    }
}
