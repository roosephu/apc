use num::{Complex, Signed};

use crate::{
    bandwidth_interp::BandwidthInterp, contexts::*, sum_trunc_dirichlet::sum_trunc_dirichlet,
    traits::MyReal,
};

fn isolate_zeros(n: usize) {
    use crate::traits::MyReal;

    let σ = 0.5;
    let lo = f64::PI() * 2.0 * (n as f64).powi(2);
    let hi = f64::PI() * 2.0 * (n as f64 + 1.0).powi(2);
    let eps = 1e-13;
    let dir = BandwidthInterp::new(n, σ);
    println!("height from {:.0} to {:.0}", lo, hi);

    let mut n_points = ((hi.rs_theta(eps) - lo.rs_theta(eps)) / std::f64::consts::PI) as usize + 1;

    let z = |x| {
        (dir.query(x, eps) * Complex::from_polar(1.0, x.rs_theta(eps))).re * 2.0
            + x.gabcke_series(eps)
    };

    for _ in 0..6 {
        let s = Complex::new(0.5, lo);
        let dir_sum = sum_trunc_dirichlet(s, 1, n, n_points + 1, (hi - lo) / n_points as f64);

        let mut signs = vec![];
        for i in 0..=n_points {
            let ratio = i as f64 / n_points as f64;
            let x = ratio * hi + (1.0 - ratio) * lo;
            let z = (dir_sum[i] * Complex::from_polar(1.0, x.rs_theta(eps))).re * 2.0
                + x.gabcke_series(eps);
            let dir_query = dir.query(x, eps);
            assert!(
                (dir_query - dir_sum[i]).norm() <= x * eps,
                "x = {}, query = {}, func = {}, z = {}",
                x,
                dir_query,
                dir_sum[i],
                z
            );
            signs.push(z.is_sign_positive());
        }
        let n_roots = (0..n_points).filter(|&i| signs[i] != signs[i + 1]).count();
        println!("# points = {}, # roots = {}", n_points, n_roots);

        n_points *= 2;
    }

    // let mut signs = vec![];
    // for i in 0..=n_points {
    //     let ratio = i as f64 / n_points as f64;
    //     let x = lo + (hi - lo) * ratio;
    //     signs.push(z(x).is_sign_positive());
    // }

    // for _ in 0..1 {
    //     let mut new_signs = vec![];
    //     for i in 0..n_points {
    //         new_signs.push(signs[i]);
    //         let ratio = (i * 2 + 1) as f64 / (n_points * 2) as f64;
    //         let x = lo + (hi - lo) * ratio;
    //         new_signs.push(z(x).is_sign_positive());
    //     }
    //     new_signs.push(signs[n_points]);
    //     signs = new_signs;
    //     let n_roots = (0..n_points * 2).filter(|&i| signs[i] != signs[i + 1]).count();
    //     println!("# points = {}, # roots = {}", n_points, n_roots);
    //     n_points *= 2;
    // }
}

fn find_zeros<T: MyReal + Sinc + GabckeSeries + Contexts + RiemannSiegelTheta + Signed>(n: usize) {
    let sigma = T::mp(0.5);
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
        isolate_zeros(4110);

        let x = T::new(63463.313195167415, 0.0);
        let t = x.gabcke_series(1e-10);
        println!("series = {}", t);
        // find_zeros::<T>(100);
    }
}
