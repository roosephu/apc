use log::debug;
use num::Complex;
use rustfft::{FftNum, FftPlanner};

use crate::traits::{ExpPolyApprox, MyReal};

/// compute \sum_{j=1}^n j^{-(s + i t delta)} for every 0 <= t <= m
pub fn sum_trunc_dirichlet<T: ExpPolyApprox + MyReal + FftNum>(
    s: Complex<T>,
    n: usize,
    m: usize,
    delta: T,
) -> Vec<Complex<T>> {
    debug!("[OS-FKBJ] s = {:.6}, n = {}, m = {}, delta = {:.6}", s, n, m, delta);
    let M2 = (m + m % 2) / 2;
    let TM2 = (M2 as i64).unchecked_cast::<T>();
    let s = s + Complex::new(T::zero(), TM2 * delta);

    // ans[t] = \sum_{j=0}^{n - 1} a[j] * exp(i gamma[j] (t - M2))
    let a: Vec<_> = (1..=n).map(|x| (-s * (x as f64).unchecked_cast::<T>().ln()).exp()).collect();
    let gamma: Vec<T> = (1..=n).map(|x| -delta * (x as f64).unchecked_cast::<T>().ln()).collect();

    // gamma[j] = w[j] * div + d[j];
    let R = (m + 1).next_power_of_two() * 2;
    let div = T::TAU() / (R as f64);
    let w: Vec<i64> = gamma.iter().map(|&x| ((x / div).round().unchecked_cast())).collect();
    let d: Vec<T> =
        gamma.iter().zip(w.iter()).map(|(&x, &y)| x - y.unchecked_cast::<T>() * div).collect();

    let mut ret = vec![Complex::zero(); m + 1];
    let mut f = vec![Complex::zero(); R];
    let mut planner = FftPlanner::<T>::new();
    let fft = planner.plan_fft_inverse(R);

    // for t in 0..=m {
    //     for j in 0..n {
    //         // ret[t] += a[j] * Complex::new(T::zero(), gamma[j] * (t as i64 - M2 as i64) as f64).exp();
    //         // ret[t] += a[j] * Complex::new(T::zero(), (div * w[j] as f64 + d[j]) * (t as f64 - M2 as f64)).exp();
    //         let rem = d[j] * (t as f64 - M2 as f64);
    //         let mut exp = Complex::<T>::zero();
    //         for (e, c) in T::get_poly_approx() {
    //             exp += c * rem.pow(e as i32);
    //         }
    //         let gt = Complex::new(T::zero(), rem).exp();
    //         // println!("diff norm = {:.e}, exp = {:?}, gt = {:?}, rem = {:?}", (gt - exp).norm(), exp, gt, rem);
    //         ret[t] += a[j]
    //             * Complex::new(T::zero(), (div * w[j] as f64) * (t as f64 - M2 as f64)).exp()
    //             * exp;
    //     }
    // }
    // return ret;

    for (e, c) in T::get_poly_approx() {
        for x in f.iter_mut() {
            x.set_zero()
        }
        for j in 0..n {
            f[w[j].rem_euclid(R as i64) as usize] += a[j] * (d[j] * TM2).pow(e as i32);
        }
        fft.process(f.as_mut_slice());
        for x in 0..=m {
            ret[x] += c
                * f[(x + R - M2) % R]
                * ((x as i32).unchecked_cast::<T>() / TM2 - 1.0).pow(e as i32);
        }
    }

    ret
}

#[cfg(test)]
mod tests {
    use crate::test_utils::*;
    use crate::{f64x2, unchecked_cast::UncheckedCast};

    use super::*;

    type T = f64x2;

    #[test]
    fn test_sum_trunc_dirichlet() {
        let N = 20;
        let s = Complex::new(1.3.unchecked_cast::<T>(), 260.unchecked_cast::<T>());
        let M = 9;
        let delta = 1.1.unchecked_cast::<T>();
        let result = sum_trunc_dirichlet(s, N, M, delta);
        for t in 0..=M {
            let mut sum = Complex::<T>::zero();
            let z = s + Complex::new(T::zero(), delta * t as f64);
            for j in 1..=N {
                sum += (-z * (j as f64).unchecked_cast::<T>().ln()).exp();
            }
            println!("{:e} {:e}, diff = {:e}", sum, result[t], sum - result[t]);
            assert_complex_close(sum, result[t], 1e-29);
        }
    }
}
