use log::debug;
use num::{Complex, Signed};
use F64x2::traits::FpOps;

use crate::{
    bandwidth_interp::BandwidthInterp,
    brentq::{brentq2, BrentqStatus},
    contexts::*,
    sum_trunc_dirichlet::sum_trunc_dirichlet,
    traits::MyReal,
};

/// Sketch of the algorithm
/// 1. We focus on Gram points
/// 2. Each time we find a Rosser block, and assume it satisfies Rosser's rule
/// 3. If it contains one Gram interval: done isolating
/// 4. Otherwise, we find one
fn isolate_zeros(n: usize) {
    use crate::traits::MyReal;

    let σ = 0.5;
    let lo = f64::PI() * 2.0 * (n as f64).powi(2);
    let hi = f64::PI() * 2.0 * (n as f64 + 1.0).powi(2);
    let eps = 1e-13;
    let dir = BandwidthInterp::new(n, σ);

    let mut n_points = ((hi.rs_theta(eps) - lo.rs_theta(eps)) / std::f64::consts::PI) as usize + 1;

    let _z = |x| {
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

use super::brentq::EvalPoint;

fn count_variations(block: &[EvalPoint<f64>]) -> usize {
    let mut variations = 0;
    for i in 1..block.len() {
        let a = &block[i];
        let b = &block[i - 1];
        if !a.has_same_sign(b) {
            variations += 1;
        }
    }
    variations
}

#[inline]
fn try_separate_same_sign(
    block: &mut Vec<EvalPoint<f64>>,
    p: &EvalPoint<f64>,
    q: &EvalPoint<f64>,
    rsz: &mut RiemannSiegelZ,
) {
    // We want to solve f(x) = 0 in [a, b] but now f(a) f(b) > 0.
    // So we model f(x) as $f(x) = p (x - q)^2$ and use $q$ here.
    let t = (q.f / p.f).sqrt();
    let g = EvalPoint::new((t * p.x + q.x) / (t + 1.0), |x| rsz.query(x, QueryMode::Separate));
    let should_dfs = g.has_same_sign(p) && g.f.abs() < p.f.abs() && g.f.abs() < q.f.abs();
    if should_dfs {
        try_separate_same_sign(block, p, &g, rsz);
    }
    block.push(g);
    if should_dfs {
        try_separate_same_sign(block, &g, q, rsz);
    }
}

#[inline]
fn try_separate_diff_sign(
    block: &mut Vec<EvalPoint<f64>>,
    p: &EvalPoint<f64>,
    q: &EvalPoint<f64>,
    rsz: &mut RiemannSiegelZ,
) {
    block.push(EvalPoint::new((p.x + q.x) / 2.0, |x| rsz.query(x, QueryMode::Separate)));
}

/// Given a Rosser block, try to separate all zeros in the block.
/// `n_zeros`: Expected number of zeros in the block
fn try_separate_rosser_block(
    mut block: Vec<EvalPoint<f64>>,
    n_zeros: usize,
    rsz: &mut RiemannSiegelZ,
) -> (bool, Vec<EvalPoint<f64>>) {
    for _ in 0..4 {
        if count_variations(&block) == n_zeros {
            return (true, block);
        }
        let mut new_block = vec![block[0]];
        for i in 1..block.len() {
            let p = &block[i - 1];
            let q = &block[i];
            if p.has_same_sign(q) {
                try_separate_same_sign(&mut new_block, p, q, rsz);
            } else {
                try_separate_diff_sign(&mut new_block, p, q, rsz);
            }
            new_block.push(*q);
        }
        block = new_block;
    }
    (false, block)
}

#[inline]
fn approx_lambert_w(x: f64) -> f64 {
    assert!(x >= f64::E());
    let r = x.ln();
    let l = r - r.ln();
    let f = |t: f64| t.exp() * t - x;
    let result = brentq2(f, EvalPoint::new(l, f), EvalPoint::new(r, f), 1e-11, 0.0, 20);
    assert!(result.status == BrentqStatus::Converged);
    result.x
}

#[inline]
fn calc_gram_point(n: usize, eps: f64) -> f64 {
    let gram = |x: f64| x.rs_theta(eps) - n as f64 * f64::PI();
    let t = (n as f64 + 0.125) / f64::E();
    assert!(t >= f64::E());
    let w = approx_lambert_w(t.fp());
    let x0 = f64::PI() * 2.0 * f64::E() * t / w;
    let a = x0 - 1.0;
    let b = x0 + 1.0;
    let fa = gram(a);
    let fb = gram(b);
    // debug!("n = {n}, a = {a:.6}, b = {b:.6}");
    let result = crate::brentq::brentq2(
        gram,
        EvalPoint { x: a, f: fa },
        EvalPoint { x: b, f: fb },
        eps,
        0.0,
        100,
    );
    // assert!(result.status == BrentqStatus::Converged, "{:?}", result); // TODO: find good stop criterion
    result.x
}

#[inline]
fn is_good_gram_point(n: usize, g: &EvalPoint<f64>) -> bool {
    g.f.is_sign_positive() == (n % 2 == 0)
}

/// returns whether if g_n is a good Gram point
fn next_rosser_block(
    mut n: usize,
    g: EvalPoint<f64>,
    rsz: &mut RiemannSiegelZ,
) -> Vec<EvalPoint<f64>> {
    assert!(is_good_gram_point(n, &g));
    // assume g_n is a good Gram point
    let mut ret = vec![g];
    loop {
        n += 1;
        let x = calc_gram_point(n, 1e-13);
        let g = EvalPoint::new(x, |x| rsz.query(x, QueryMode::Separate));
        ret.push(g);
        if is_good_gram_point(n, &g) {
            // (-1)^n Z(g_n) > 0
            return ret;
        }
    }
}

fn try_isolate(rsz: &mut RiemannSiegelZ, xtol: f64, _rtol: f64) -> Vec<f64> {
    let mut n = 10010;
    let x = calc_gram_point(n, 1e-13);
    let mut g = EvalPoint::new(x, |x| rsz.query(x, QueryMode::Separate));
    assert!(is_good_gram_point(n, &g));

    let mut roots = vec![];
    while n < 110000 {
        let block = next_rosser_block(n, g, rsz);
        let n_zeros = block.len() - 1;

        let (separated, block) = try_separate_rosser_block(block, n_zeros, rsz);
        if separated {
            // yeah!
            for i in 1..block.len() {
                let a = block[i - 1];
                let b = block[i];
                if !a.has_same_sign(&b) {
                    let f = |x| rsz.query(x, QueryMode::Locate);
                    let result = crate::brentq::brentq2(f, a, b, xtol, 0.0, 100);
                    assert!(result.status == BrentqStatus::Converged);
                    roots.push(result.x);
                }
            }
        } else {
            panic!(
                "Exception of Rosser's rule! block = {:?}, variations = {}",
                block,
                count_variations(&block)
            );
        }

        n += n_zeros;
        g = block[block.len() - 1];
        let n0 = (block[0].x.fp() / 2.0 / f64::PI()).sqrt().floor() as usize;
        let n1 = (g.x.fp() / 2.0 / f64::PI()).sqrt().floor() as usize;
        if n0 != n1 {
            rsz.drop(n0);
        }
    }
    roots
    // debug!("{} zeros found between {} and {}", zeta_zeros.len(), zeta_zeros[0], zeta_zeros.last().unwrap());
}

pub enum QueryMode {
    Separate,
    Locate,
}

pub struct RiemannSiegelZ {
    pub dirichlet: Vec<Option<BandwidthInterp<f64>>>,
    pub eps: f64,
    pub counts: [usize; 2],
    pub mode: QueryMode,
}

impl RiemannSiegelZ {
    pub fn new(max_height: f64, eps: f64) -> Self {
        let n = (max_height / 2.0 / f64::PI()).sqrt().ceil() as usize + 10;
        let dirichlet = (0..n).map(|_| None).collect::<Vec<_>>();
        Self { dirichlet, eps, counts: [0, 0], mode: QueryMode::Separate }
    }

    fn get(&mut self, n: usize) -> &BandwidthInterp<f64> {
        if self.dirichlet[n].is_none() {
            self.dirichlet[n] = Some(BandwidthInterp::new(n, 0.5));
        }
        self.dirichlet[n].as_ref().unwrap()
    }

    pub fn level(&self, x: f64) -> usize { (x / f64::PI() / 2.0).sqrt().floor() as usize }

    pub fn query(&mut self, x: f64, mode: QueryMode) -> f64 {
        self.counts[mode as usize] += 1;
        let n = self.level(x);
        let eps = self.eps;
        let dir = self.get(n);
        (dir.query(x, eps) * Complex::from_polar(1.0, x.rs_theta(eps))).re * 2.0
            + x.gabcke_series(eps)
    }

    pub fn drop(&mut self, n: usize) {
        for i in 0..n {
            self.dirichlet[i] = None;
        }
    }
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
    use F64x2::test_utils::assert_close;

    use super::*;
    use crate::types::T;

    #[test]
    fn test_zeta_zeros() {
        use crate::contexts::*;

        crate::init();
        isolate_zeros(4110);

        let x = T::new(63463.313195167415, 0.0);
        let t = x.gabcke_series(1e-10);
        println!("series = {}", t);
        // find_zeros::<T>(100);
    }

    #[test]
    fn test_try_isolate() {
        crate::init();
        let mut rsz = RiemannSiegelZ::new(1e5, 1e-12);

        let roots = try_isolate(&mut rsz, 1e-5, 1e-30);
        let n_calls_separate = rsz.counts[0];
        let n_calls_locate = rsz.counts[0];
        let n_zeros = roots.len();
        debug!(
            "To separate {} zeros: {:.3} calls to separate, {:.3} calls to locate",
            n_zeros,
            n_calls_separate as f64 / n_zeros as f64,
            n_calls_locate as f64 / n_zeros as f64
        );

        panic!();
    }

    #[test]
    fn test_gram() {
        crate::init();
        let eps1 = 1e-11;
        let eps2 = 1e-9;
        assert_close(calc_gram_point(100, eps1), 238.5825905145, eps2);
    }
}
