use std::borrow::Borrow;

use super::brentq::EvalPoint;
use log::debug;
use num::{Complex, Signed};
use num_traits::FloatConst;
use F64x2::traits::FpOps;

use crate::{
    bandwidth_interp::BandwidthInterp,
    brentq::{brentq2, Brentq, BrentqStatus},
    contexts::*,
    illinois::Illinois,
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

pub struct SameSign<T: MyReal> {
    xa: EvalPoint<T>,
    xb: EvalPoint<T>,
}

pub enum BracketMaintainer<T: MyReal> {
    DiffSign(Illinois<T>),
    SameSign(SameSign<T>),
}

pub struct Bracket<T: MyReal> {
    entry: BracketMaintainer<T>,
    priority: f64,
}

impl<T: MyReal> Bracket<T> {
    fn new_diffsign(entry: Illinois<T>) -> Bracket<T> {
        let priority = (entry.xa.x.fp() - entry.xb.x.fp()).abs();
        Self { entry: BracketMaintainer::DiffSign(entry), priority }
    }

    fn new_samesign(entry: SameSign<T>) -> Bracket<T> {
        let priority = (entry.xa.x.fp() - entry.xb.x.fp()).abs();
        Self { entry: BracketMaintainer::SameSign(entry), priority }
    }

    fn get_interval(&self) -> (EvalPoint<T>, EvalPoint<T>) {
        match &self.entry {
            BracketMaintainer::DiffSign(x) => (x.xa, x.xb),
            BracketMaintainer::SameSign(x) => (x.xa, x.xb),
        }
    }
}

/// A `Lounge` is a place where roots are waiting to be separated.
struct Lounge<T: MyReal> {
    brackets: Vec<Bracket<T>>,
    variations: usize,
    xtol: f64,
    rtol: f64,
}

impl<T: MyReal> Lounge<T> {
    pub fn new(xtol: f64, rtol: f64) -> Self {
        Self { brackets: vec![], variations: 0, xtol, rtol }
    }

    pub fn add(&mut self, xa: EvalPoint<T>, xb: EvalPoint<T>) {
        if xa.has_same_sign(&xb) {
            self.brackets.push(Bracket::new_samesign(SameSign { xa, xb }));
        } else {
            self.brackets.push(Bracket::new_diffsign(Illinois::new(xa, xb)));
            self.variations += 1;
        }
    }

    pub fn count_variations(&mut self) -> usize {
        let mut ret = 0;
        for bracket in &self.brackets {
            if let BracketMaintainer::DiffSign(_) = bracket.entry {
                ret += 1;
            }
        }
        ret
    }

    fn select(&mut self) -> Bracket<T> {
        let mut cur_idx = 0;
        let mut cur_priority = f64::MIN;
        for (i, bracket) in self.brackets.iter().enumerate() {
            if bracket.priority > cur_priority {
                cur_priority = bracket.priority;
                cur_idx = i;
            }
        }
        // debug!(
        //     "selected interval = {:?}, priority = {}",
        //     self.brackets[cur_idx].get_interval(),
        //     cur_priority
        // );
        self.brackets.remove(cur_idx)
    }

    pub fn try_isolate(
        &mut self,
        n_zeros: usize,
        n_iters: usize,
        mut f: impl FnMut(T) -> T,
    ) -> bool {
        for _ in 0..n_iters {
            if self.variations == n_zeros {
                break;
            }
            // select the one bracket
            let bracket = self.select();

            // split the bracket
            match bracket.entry {
                BracketMaintainer::DiffSign(mut x) => {
                    let xa = x.xa;
                    let xb = x.xb;
                    x.step(&mut f);
                    if x.xa.x != xa.x {
                        assert!(x.xb.x == xb.x);
                        self.add(xa, x.xa);
                    } else {
                        assert!(x.xa.x == xa.x);
                        self.add(xb, x.xb);
                    }
                    self.brackets.push(Bracket::new_diffsign(x));
                }
                BracketMaintainer::SameSign(x) => {
                    let (p, q) = (x.xa, x.xb);
                    let g = if rand::random::<u8>() % 2 == 0 {
                        let t = (q.f / p.f).sqrt();
                        EvalPoint::new((t * p.x + q.x) / (t + 1.0), &mut f)
                    } else {
                        EvalPoint::new((p.x + q.x) * 0.5, &mut f)
                    };
                    // let g = EvalPoint::new((p.x + q.x) * 0.5, &mut f);
                    self.add(x.xa, g);
                    self.add(g, x.xb);
                }
            }
        }
        self.variations == n_zeros
    }
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
#[inline(never)]
fn gram_point<T: MyReal + RiemannSiegelTheta>(n: usize, eps: f64) -> T {
    assert!(n <= (1usize << 52), "`n as f64` has rounding error");

    let gram = |x: T| x.rs_theta(eps) - T::PI() * n as f64;
    let t = (T::mp(n as f64) + 0.125) / T::E();
    assert!(t >= T::E());
    let w = approx_lambert_w(t.fp());
    let x0 = T::PI() * 2.0 * T::E() * t / w;
    let a = x0 - 1.0;
    let b = x0 + 1.0;
    let result =
        Brentq::new(EvalPoint::new(a, gram), EvalPoint::new(b, gram), eps, 0.0).solve(gram, 100);
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

pub trait RiemannSiegelZReq =
    MyReal + Sinc + ExpPolyApprox + RiemannSiegelTheta + GabckeSeries + Bernoulli + Factorial;

/// Apply Euler-Maclaurin formula to compute $\zeta(1/2 + i t)$.
/// See Eq (1.2) and (1.3) in Fast Algorithm For Multiple Evaluation Of
/// The Riemann Zeta Function, A.M. Odlyzko and A. Schonhage.
///
/// Although Euler-Maclaurin method requires O(t) time, the $n$ we used is
/// fixed as we only apply Euler-Maclaurin to small $t$ where Riemann-Siegel
/// can't provide enough precision, and larger $n$ always works better.

pub struct EulerMaclaurinMethod<T> {
    dirichlet: BandwidthInterp<T>,
    n: usize,
    atol: f64,
}

impl<T: RiemannSiegelZReq> EulerMaclaurinMethod<T> {
    pub fn new(n: usize, atol: f64) -> Self {
        let dirichlet =
            BandwidthInterp::<T>::new(n - 1, T::zero(), T::mp(n as f64), T::mp(0.5), atol);
        Self { dirichlet, n, atol }
    }

    pub fn query(&self, x: T) -> T {
        let atol = self.atol;

        let n = self.n; // or another?
        let s = Complex::new(T::mp(0.5), x);
        let n_t = T::mp(n as f64);
        let ln_n = n_t.ln();
        let n_pow_minus_s = (-s * ln_n).exp();
        let mut zeta = self.dirichlet.query(x)
            + n_pow_minus_s * T::mp(0.5)
            + n_pow_minus_s * n_t / (s - T::one()); // TODO: merge two

        let mut term = n_pow_minus_s / n_t * s;
        let n_sqr = T::mp((n * n) as f64);
        for k in 1..=n {
            let value = term * T::bernoulli(2 * k) / T::factorial(2 * k);
            let error = value * (s + T::mp((2 * k + 1) as f64)) / T::mp((2 * k + 1) as f64 + 0.5);
            if error.norm().fp() < atol {
                break;
            }
            zeta += value;
            term = term / n_sqr * (s + T::mp((2 * k - 1) as f64)) * (s + T::mp((2 * k) as f64));
        }
        (zeta * Complex::from_polar(T::mp(1.0), x.rs_theta(atol))).re
    }
}

/// Compute HardyZ function by different methods: Euler-Maclaurin method
/// or Riemann-Siegel method.
pub struct HardyZ<T: RiemannSiegelZReq> {
    pub dirichlet: Vec<Option<BandwidthInterp<T>>>,
    pub eps: f64,
    pub counts: [usize; 2],
    pub levels: Vec<(f64, usize)>,
    pub euler_maclaurin: EulerMaclaurinMethod<T>,
}

/// TODO: mixed precision
impl<T: RiemannSiegelZReq> HardyZ<T> {
    pub fn new(max_height: f64, max_order: usize, eps: f64) -> Self {
        let n = (max_height / 2.0 / f64::PI()).sqrt().ceil() as usize + 10;
        let dirichlet = (0..=n).map(|_| None).collect::<Vec<_>>();
        let (levels, em_height) = Self::get_plans(max_height, max_order, eps);
        let euler_maclaurin = EulerMaclaurinMethod::new(em_height as usize, eps);
        Self { dirichlet, eps, counts: [0, 0], levels, euler_maclaurin }
    }

    fn get_plans(max_height: f64, max_order: usize, eps: f64) -> (Vec<(f64, usize)>, f64) {
        const COEFFS: [f64; 11] =
            [0.127, 0.053, 0.011, 0.031, 0.017, 0.061, 0.661, 9.2, 130.0, 1837.0, 25966.0];
        assert!(
            max_order < COEFFS.len(),
            "don't have enough information about Riemann-Siegel formula"
        );
        let mut levels: Vec<(f64, usize)> = vec![];
        let mut cur_t = max_height;
        for k in 0..=max_order {
            let min_t = (COEFFS[k] / eps).powf(4.0 / (3.0 + k as f64 * 2.0));
            if min_t < cur_t {
                levels.push((min_t, k));
                cur_t = min_t;
            }
        }
        (levels, cur_t)
    }

    fn get(&mut self, n: usize) -> &BandwidthInterp<T> {
        if self.dirichlet[n].is_none() {
            let min_t = T::mp(f64::PI() * 2.0 * (n as f64).powi(2));
            let max_t = T::mp(f64::PI() * 2.0 * (n as f64 + 1.0).powi(2));
            self.dirichlet[n] = Some(BandwidthInterp::new(n, min_t, max_t, T::mp(0.5), self.eps));
        }
        self.dirichlet[n].as_ref().unwrap()
    }

    pub fn level(&self, x: T) -> usize { (x.fp() / f64::PI() / 2.0).sqrt().floor() as usize }

    pub fn riemann_siegel(&mut self, x: T, order: usize) -> T {
        let eps = self.eps;
        let n = self.level(x);
        let dir = self.get(n);
        (Complex::from_polar(T::mp(1.0), x.rs_theta(eps)) * dir.query(x)).re * 2.0
            + x.gabcke_series(order, eps)
    }

    pub fn query(&mut self, x: T, mode: QueryMode) -> T {
        self.counts[mode as usize] += 1;
        let x_fp = x.fp();
        for &(min_t, order) in &self.levels {
            if x_fp >= min_t {
                return self.riemann_siegel(x, order);
            }
        }
        self.euler_maclaurin.query(x)
    }

    pub fn purge(&mut self) { todo!() }
}

/// Sketch of the algorithm
/// 1. We focus on Gram points
/// 2. Each time we find a Rosser block, and assume it satisfies Rosser's rule
/// 3. If it contains one Gram interval: done isolating
/// 4. Otherwise, we find one
pub fn try_isolate<T: RiemannSiegelZReq>(
    rsz: &mut HardyZ<T>,
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

        let mut lounge = Lounge::new(xtol, 0.0);
        for i in 0..n_zeros {
            lounge.add(block[i], block[i + 1]);
        }
        let separated = lounge.try_isolate(n_zeros, 1000, |x| rsz.query(x, QueryMode::Separate));
        let n_blocks_evaled = lounge.brackets.len();

        if n_blocks_evaled > max_n_intervals {
            max_n_intervals = n_blocks_evaled;
            debug!("{max_n_intervals} blocks!");
        }
        if separated {
            // yeah!
            for bracket in lounge.brackets {
                if let BracketMaintainer::DiffSign(x) = bracket.entry {
                    let result = Brentq::new(x.xa, x.xb, xtol, 0.0)
                        .solve(|x| rsz.query(x, QueryMode::Locate), 100);
                    assert!(result.status == BrentqStatus::Converged);
                    roots.push(result.x);
                }
            }
        } else {
            panic!(
                "Exception of Rosser's rule! n = {n}, block = {:?}, expect variations = {}, found = {}",
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
    use F64x2::f64x2;

    #[test]
    fn test_gram() {
        crate::init();
        let eps1 = 1e-11;
        let eps2 = 1e-9;
        assert_close(gram_point(100, eps1), 238.5825905145, eps2, 0.0);
    }

    #[test]
    fn test_root_find() {
        crate::init();
        type T = f64x2;
        let mut rsz = HardyZ::<T>::new(1e8, 10, 1e-18);

        let mut f = |x: T| rsz.query(x, QueryMode::Locate);
        let center = T::mp(74929.812159119);
        let xa = EvalPoint::new(center - 0.4, &mut f);
        let xb = EvalPoint::new(center + 0.3, &mut f);

        let rtol = 1e-18;
        let result = Brentq::new(xa, xb, 0.0, rtol).solve(f, 100);
        println!("status = {:?}, x = {},f = {}", result.status, result.x, result.f);
        let x = result.x;
        assert_close(x, center, 0.0, rtol);
    }

    #[test]
    fn test_euler_maclaurin() {
        crate::init();
        type T = f64x2;
        let atol = 1e-18;
        let euler_maclaurin = EulerMaclaurinMethod::new(1000, atol);
        let x = T::mp(1000.0);
        let z = euler_maclaurin.query(x);
        let z_gt = f64x2::new(0.9977946375215866, 2.2475077894406e-17);
        assert_close(z, z_gt, atol, 0.0);
    }
}
