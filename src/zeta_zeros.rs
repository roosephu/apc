use num::Complex;

use crate::{bandwidth_interp::BandwidthInterp, contexts::*, traits::MyReal};

fn isolate_zeros(n: usize) {
    use crate::traits::MyReal;
    let σ = 0.5;
    let lo = f64::PI() * 2.0 * (n as f64).powi(2);
    let hi = f64::PI() * 2.0 * (n as f64 + 1.0).powi(2);
    let eps = 1e-15;
    let dir = BandwidthInterp::new(n, σ);

    let mut n_points = ((hi.rs_theta(eps) - lo.rs_theta(eps)) / std::f64::consts::PI) as usize + 1;

    let z = |x| {
        (dir.query(x, eps) * Complex::new(0.0, x.rs_theta(eps)).exp()).re * 2.0
            + x.gabcke_series(eps)
    };
    let mid = (lo + hi) * 0.5;
    println!("mid = {}, Z(mid) = {}, gabcke = {}", mid, z(mid), mid.gabcke_series(eps));

    let mut signs = vec![];
    for i in 0..=n_points {
        let ratio = i as f64 / n_points as f64;
        let x = lo + (hi - lo) * ratio;
        signs.push(z(x).is_sign_positive());
    }

    for _ in 0..5 {
        let mut new_signs = vec![];
        for i in 0..n_points {
            new_signs.push(signs[i]);
            let ratio = (i * 2 + 1) as f64 / (n_points * 2) as f64;
            let x = lo + (hi - lo) * ratio;
            new_signs.push(z(x).is_sign_positive());
        }
        new_signs.push(signs[n_points]);
        signs = new_signs;
        let n_roots = signs.iter().zip(signs[1..].iter()).filter(|(&a, &b)| a != b).count();
        println!("# points = {}, # roots = {}", n_points, n_roots);
        n_points *= 2;
    }
}

fn find_zeros<T: MyReal + Sinc + GabckeSeries + Contexts + RiemannSiegelTheta>(n: usize) {
    let sigma = T::from_f64(0.5).unwrap();
    let dir = BandwidthInterp::new(n, sigma);
    let lo = T::PI() * (2 * n * n) as f64;
    let hi = T::PI() * (2 * (n + 1) * (n + 1)) as f64;
    println!("l = {}, r = {}", lo, hi);

    let eps = 1e-15;

    let z = |x| {
        (dir.query(x, eps) * Complex::new(T::zero(), x.rs_theta(eps)).exp()).re * 2.0
            + x.gabcke_series(eps)
    };
    let mid = (lo + hi) * 0.5;
    println!("mid = {}, Z(mid) = {}, gabcke = {}", mid, z(mid), mid.gabcke_series(eps));
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::T;

    #[test]
    fn test_zeta_zeros() {
        use crate::contexts::*;

        crate::contexts::init();
        isolate_zeros(100);

        let x = T { hi: 63463.313195167415, lo: 0.0 };
        let t = x.gabcke_series(1e-10);
        println!("series = {}", t);
        // find_zeros::<T>(100);
    }
}
