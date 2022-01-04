use super::brentq::EvalPoint;
use log::debug;
use num::{Complex, Signed};
use num_traits::FloatConst;
use F64x2::traits::FpOps;

use crate::{
    bandwidth_interp::BandwidthInterp,
    brentq::{brentq2, Brentq, BrentqStatus},
    contexts::*,
    traits::MyReal,
};

fn count_variations<T: Copy + Signed>(block: &[EvalPoint<T>]) -> usize {
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

// /// For root finding
// struct Lounge<T> {
//     intervals: Vec<T>,
// }

// impl<T: MyReal> Lounge<T> {
//     pub fn new() -> Self {
//         Self { intervals: vec![] }
//     }

//     pub fn add()
// }

/// We want to solve f(x) = 0 in [a, b] but now f(a) f(b) > 0.
/// So we model f(x) as $f(x) = p (x - q)^2$ and use $q$ here.
/// TODO: recursive?
#[inline]
fn try_separate_same_sign<T: MyReal>(
    block: &mut Vec<EvalPoint<T>>,
    p: &EvalPoint<T>,
    q: &EvalPoint<T>,
    f: impl FnMut(T) -> T,
) {
    let g;
    if rand::random::<u8>() % 2 == 0 {
        let t = (q.f / p.f).sqrt();
        g = EvalPoint::new((t * p.x + q.x) / (t + 1.0), f);
    } else {
        g = EvalPoint::new((p.x + q.x) * 0.5, f);
    }
    // let should_dfs = g.has_same_sign(p) && g.f.abs() < p.f.abs() && g.f.abs() < q.f.abs();
    // if should_dfs {
    //     try_separate_same_sign(block, p, &g, rsz);
    // }
    block.push(g);
    // if should_dfs {
    //     try_separate_same_sign(block, &g, q, rsz);
    // }
}

/// TODO: Illinois/Brentq?
#[inline]
fn try_separate_diff_sign<T: MyReal>(
    block: &mut Vec<EvalPoint<T>>,
    p: &EvalPoint<T>,
    q: &EvalPoint<T>,
    f: impl FnMut(T) -> T,
) {
    block.push(EvalPoint::new((p.x + q.x) / 2.0, f));
}

/// Given a Rosser block, try to separate all zeros in the block.
/// `n_zeros`: Expected number of zeros in the block
#[inline(never)]
fn try_separate_rosser_block<T: MyReal>(
    mut block: Vec<EvalPoint<T>>,
    n_zeros: usize,
    n_iters: usize,
    mut f: impl FnMut(T) -> T,
) -> (bool, Vec<EvalPoint<T>>) {
    for _ in 0..n_iters {
        if count_variations(&block) == n_zeros {
            return (true, block);
        }
        let mut new_block = vec![block[0]];
        for i in 1..block.len() {
            let p = &block[i - 1];
            let q = &block[i];
            if p.has_same_sign(q) {
                try_separate_same_sign(&mut new_block, p, q, &mut f);
            } else {
                try_separate_diff_sign(&mut new_block, p, q, &mut f);
            }
            new_block.push(*q);
        }
        block = new_block;
    }
    (false, block)
}

/// Compute Lambert W function $W(x)$, which is the solution to $t e^t = x$ for
/// positive $x$.
#[inline]
fn approx_lambert_w(x: f64) -> f64 {
    assert!(x >= f64::E());
    let r = x.ln();
    let l = r - r.ln();
    let f = |t: f64| t.exp() * t - x;
    let result = Brentq::new(EvalPoint::new(l, f), EvalPoint::new(r, f), 1e-11, 0.0).solve(f, 20);
    assert!(result.status == BrentqStatus::Converged);
    result.x
}

/// computes the Gram point $g_n$ with accuracy $\epsilon$.
#[inline]
fn gram_point<T: MyReal + RiemannSiegelTheta>(n: usize, eps: f64) -> T {
    assert!(n <= (1usize << 52), "`n as f64` has rounding error");

    let gram = |x: T| x.rs_theta(eps) - T::PI() * n as f64;
    let t = (T::mp(n as f64) + 0.125) / T::E();
    assert!(t >= T::E());
    let w = approx_lambert_w(t.fp());
    let x0 = T::PI() * 2.0 * T::E() * t / w;
    let a = x0 - 1.0;
    let b = x0 + 1.0;
    let result = Brentq::new(EvalPoint::new(a, gram), EvalPoint::new(b, gram), eps, 0.0)
        .solve(gram, 100);
    assert!(result.status != BrentqStatus::SignError);
    // assert!(result.status == BrentqStatus::Converged, "{:?}", result); // TODO: find good stop criterion
    result.x
}

/// (-1)^n Z(g_n) > 0
#[inline]
fn is_good_gram_point<T: Copy + Signed>(n: usize, g: &EvalPoint<T>) -> bool {
    g.f.is_positive() == (n % 2 == 0)
}

/// Returns whether if g_n is a good Gram point. Note that the
/// Gram points are not accurate: we don't need them to be accurate.
fn next_rosser_block<T: MyReal + RiemannSiegelTheta>(
    mut n: usize,
    g: EvalPoint<T>,
    mut f: impl FnMut(T) -> T,
) -> Vec<EvalPoint<T>> {
    assert!(is_good_gram_point(n, &g));
    // assume g_n is a good Gram point
    let mut ret = vec![g];
    loop {
        n += 1;
        let x = T::mp(gram_point::<f64>(n, 1e-13));
        let g = EvalPoint::new(x, &mut f);
        ret.push(g);
        if is_good_gram_point(n, &g) {
            return ret;
        }
    }
}

pub enum QueryMode {
    Separate,
    Locate,
}

pub trait RiemannSiegelZReq = MyReal + Sinc + ExpPolyApprox + RiemannSiegelTheta + GabckeSeries;

pub struct RiemannSiegelZ<T: RiemannSiegelZReq> {
    pub dirichlet: Vec<Option<BandwidthInterp<T>>>,
    pub eps: f64,
    pub counts: [usize; 2],
}

impl<T: RiemannSiegelZReq> RiemannSiegelZ<T> {
    pub fn new(max_height: f64, eps: f64) -> Self {
        let n = (max_height / 2.0 / f64::PI()).sqrt().ceil() as usize + 10;
        let dirichlet = (0..n).map(|_| None).collect::<Vec<_>>();
        Self { dirichlet, eps, counts: [0, 0] }
    }

    fn get(&mut self, n: usize) -> &BandwidthInterp<T> {
        if self.dirichlet[n].is_none() {
            self.dirichlet[n] = Some(BandwidthInterp::new(n, T::mp(0.5)));
        }
        self.dirichlet[n].as_ref().unwrap()
    }

    pub fn level(&self, x: T) -> usize { (x.fp() / f64::PI() / 2.0).sqrt().floor() as usize }

    pub fn query(&mut self, x: T, mode: QueryMode) -> T {
        self.counts[mode as usize] += 1;
        let n = self.level(x);
        let eps = self.eps;
        let dir = self.get(n);
        (Complex::from_polar(T::mp(1.0), x.rs_theta(eps)) * dir.query(x, eps)).re * 2.0
            + x.gabcke_series(eps)
    }

    pub fn purge(&mut self) { todo!() }
}

/// Sketch of the algorithm
/// 1. We focus on Gram points
/// 2. Each time we find a Rosser block, and assume it satisfies Rosser's rule
/// 3. If it contains one Gram interval: done isolating
/// 4. Otherwise, we find one
pub fn try_isolate<T: RiemannSiegelZReq>(
    rsz: &mut RiemannSiegelZ<T>,
    n0: usize,
    n1: usize,
    xtol: f64,
    _rtol: f64,
) -> Vec<T> {
    let mut n = n0;
    let x = gram_point(n, xtol);
    let mut g = EvalPoint::new(x, |x| rsz.query(x, QueryMode::Separate));
    assert!(is_good_gram_point(n, &g));

    let mut max_n_intervals = 0;

    let mut roots = vec![];
    while n < n1 {
        let block = next_rosser_block(n, g, |x| rsz.query(x, QueryMode::Separate));
        let initial_block = block.clone();
        let n_zeros = block.len() - 1;

        let (separated, block) =
            try_separate_rosser_block(block, n_zeros, 10, |x| rsz.query(x, QueryMode::Separate));
        if block.len() > max_n_intervals {
            max_n_intervals = block.len();
            debug!("{max_n_intervals} blocks!");
        }
        if separated {
            // yeah!
            for i in 1..block.len() {
                let a = block[i - 1];
                let b = block[i];
                if !a.has_same_sign(&b) {
                    let result = Brentq::new(a, b, xtol, 0.0)
                        .solve(|x| rsz.query(x, QueryMode::Locate), 100);
                    assert!(result.status == BrentqStatus::Converged);
                    roots.push(result.x);
                }
            }
        } else {
            panic!(
                "Exception of Rosser's rule! n = {n}, block = {:?}, expect variations = {}, variations = {}",
                initial_block,
                n_zeros,
                count_variations(&block),
            );
        }

        n += n_zeros;
        g = block[block.len() - 1];
    }
    // rsz.purge();
    roots
}

#[cfg(test)]
mod tests {
    use F64x2::test_utils::assert_close;

    use super::*;
    use crate::types::T;

    #[test]
    fn test_gram() {
        crate::init();
        let eps1 = 1e-11;
        let eps2 = 1e-9;
        assert_close(gram_point(100, eps1), 238.5825905145, eps2);
    }
}
