use num::{ToPrimitive, Complex};

use crate::{traits::{MyReal, GabckeExpansion}, bandwidth_interp::BandwidthInterp, rs_theta::RiemannSiegelTheta, contexts::*};
use crate::types::T;

// TODO: determine K wisely
fn calc_gabcke_n_terms<T: MyReal>(t: T, eps: f64) -> usize {
    7
}

fn gabcke_series<T: MyReal + GabckeExpansion>(t: T, eps: f64) -> T {
    let a = (t / T::PI() / 2.0).sqrt();
    let n = a.floor();
    let K = calc_gabcke_n_terms(t, eps);
    T::expand(a, T::one() - (a - n) * 2.0, K, eps)
}

fn find_zeros<T: MyReal + Sinc + GabckeExpansion + Contexts>(n: usize) {
    let sigma = T::from_f64(0.5).unwrap();
    let dir = BandwidthInterp::new(n, sigma);
    let lo = T::PI() * (2 * n * n) as f64;
    let hi = T::PI() * (2 * (n + 1) * (n + 1)) as f64;
    let theta = RiemannSiegelTheta::new(20);
    println!("l = {}, r = {}", lo, hi);

    let eps = 1e-15;
    let mut n_points = ((theta.theta(hi, eps) - theta.theta(lo, eps)) / T::PI()).to_usize().unwrap() + 1;

    let z = |x| (dir.query(x, eps) * Complex::new(T::zero(), theta.theta(x, eps)).exp()).re * 2.0 + gabcke_series(x, eps) / (x / T::PI() / 2.0).sqrt().sqrt() * (if n % 2 == 0 { -1.0 } else { 1.0 });
    let mid = (lo + hi) * 0.5;
    println!("mid = {}, Z(mid) = {}, gabcke = {}", mid, z(mid), gabcke_series(mid, eps));

    let mut signs = vec![];
    for i in 0..=n_points {
        let ratio = T::from_usize(i).unwrap() / T::from_usize(n_points).unwrap();
        let x = lo + (hi - lo) * ratio;
        signs.push(z(x).is_sign_positive());
    }

    // loop {
    for _ in 0..5 {
        let mut new_signs = vec![];
        for i in 0..n_points {
            new_signs.push(signs[i]);
            let ratio = T::from_usize(i * 2 + 1).unwrap() / T::from_usize(n_points * 2).unwrap();
            let x = lo + (hi - lo) * ratio;
            new_signs.push(z(x).is_sign_positive());
        }
        new_signs.push(signs[n_points]);
        signs = new_signs;
        let n_roots: usize = signs.iter().zip(signs[1..].iter()).map(|(&a, &b)| (a != b) as usize).sum();
        println!("# points = {}, # roots = {}", n_points, n_roots);
        n_points *= 2;
        // if n_points >= 100 {
        //     break
        // }
    }
}

#[cfg(test)]
mod tests {
    use crate::types::T;
    use super::find_zeros;

    #[test]
    fn test_zeta_zeros() {
        crate::contexts::init();
        find_zeros::<T>(100);
    }
}
