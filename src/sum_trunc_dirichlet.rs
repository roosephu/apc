use log::debug;
use num::Complex;
use rustfft::{FftNum, FftPlanner};

use crate::traits::{ExpPolyApprox, MyReal};

/// compute F(t) = \sum_{(a, g)} a exp(i t g) for t in [-m, m]
pub(crate) fn sum_weighted_exp<T: ExpPolyApprox + MyReal + FftNum>(
    a: &[Complex<T>],
    g: &[T],
    m: usize,
) -> Vec<Complex<T>> {
    let TM = (m as f64).unchecked_cast::<T>();
    // g[j] = w[j] * div + d[j];
    let n = a.len();
    let R = (m as usize + 1).next_power_of_two() * 2;
    let div = T::TAU() / (R as f64);
    let w: Vec<i64> = g.iter().map(|&x| ((x / div).round().unchecked_cast())).collect();
    let d: Vec<T> =
        g.iter().zip(w.iter()).map(|(&x, &y)| x - y.unchecked_cast::<T>() * div).collect();

    let mut ret = vec![Complex::zero(); m * 2 + 1];
    let mut f = vec![Complex::zero(); R];
    let mut planner = FftPlanner::<T>::new();
    let fft = planner.plan_fft_inverse(R);

    for (e, c) in T::get_poly_approx() {
        for x in f.iter_mut() {
            x.set_zero()
        }
        for j in 0..n {
            f[w[j].rem_euclid(R as i64) as usize] += a[j] * (d[j] * TM).pow(e as i32);
        }
        fft.process(f.as_mut_slice());
        for x in 0..=m * 2 {
            ret[x] += c
                * f[(x + R - m) % R]
                * ((x as i32).unchecked_cast::<T>() / TM - 1.0).pow(e as i32);
        }
    }

    ret
}

/// compute \sum_{j=n0}^n1 j^{-(s + i t delta)} for every 0 <= t <= m
pub(crate) fn sum_trunc_dirichlet<T: ExpPolyApprox + MyReal + FftNum>(
    s: Complex<T>,
    n0: usize,
    n1: usize,
    m: usize,
    delta: T,
) -> Vec<Complex<T>> {
    debug!("[OS-FKBJ] s = {:.6}, n = [{}, {}], m = {}, delta = {:.6}", s, n0, n1, m, delta);
    let M2 = (m + m % 2) / 2;
    let s = s + Complex::new(T::zero(), delta * M2 as f64);

    // ans[t] = \sum_{j=0}^{n - 1} a[j] * exp(i gamma[j] (t - M2))
    let ln_x: Vec<_> = (n0..=n1).map(|x| (x as f64).unchecked_cast::<T>().ln()).collect();
    let a: Vec<_> = ln_x.iter().map(|&ln_x| (-s * ln_x).exp()).collect();
    let gamma: Vec<T> = ln_x.iter().map(|&ln_x| -delta * ln_x).collect();

    sum_weighted_exp(&a, &gamma, M2)
}

#[cfg(test)]
mod tests {
    use crate::test_utils::*;
    use crate::unchecked_cast::UncheckedCast;
    use F64x2::f64x2;

    use super::*;

    type T = f64x2;

    #[test]
    fn test_sum_trunc_dirichlet() {
        let res: Option<()> = try {
            let N = 20;
            let s = Complex::new( T::from_f64(1.3)?, 260.unchecked_cast::<T>());
            let M = 9;
            let delta = T::from_f64(1.1)?;
            let result = sum_trunc_dirichlet(s, 1, N, M, delta);
            for t in 0..=M {
                let mut sum = Complex::<T>::zero();
                let z = s + Complex::new(T::zero(), delta * t as f64);
                for j in 1..=N {
                    sum += (-z * T::from_f64(j as f64)?.ln()).exp();
                }
                println!("{:e} {:e}, diff = {:e}", sum, result[t], sum - result[t]);
                assert_complex_close(sum, result[t], 1e-29);
            }
        };
        res.unwrap();
    }
}
