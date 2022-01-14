use std::{borrow::Borrow, cmp::Ordering, collections::VecDeque, ops::Range};

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
    n_zeros: usize,
    atol: f64,
    rtol: f64,
}

impl<T: MyReal> Lounge<T> {
    pub fn new(n_zeros: usize, atol: f64, rtol: f64) -> Self {
        Self { brackets: vec![], variations: 0, n_zeros, atol, rtol }
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
        self.brackets.remove(cur_idx)
    }

    pub fn try_isolate(&mut self, n_iters: usize, mut f: impl FnMut(T) -> T) {
        for _ in 0..n_iters {
            if self.separated() {
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
    }

    #[inline]
    pub fn separated(&self) -> bool { self.variations == self.n_zeros }

    pub fn merge(&mut self, l: Self) {
        self.brackets.extend(l.brackets);
        self.variations += l.variations;
        self.n_zeros += l.n_zeros;
    }

    pub fn bounds(&self) -> Range<f64> {
        self.brackets.iter().fold(f64::MAX..f64::MIN, |r, x| {
            let interval = x.get_interval();
            r.start.min(interval.0.x.fp())..r.end.max(interval.1.x.fp())
        })
    }
}

/// computes the Gram point $g_n$ with accuracy $\epsilon$.
#[inline(never)]
fn gram_point<T: MyReal + RiemannSiegelTheta>(n: i64, atol: f64) -> T {
    assert!(-1 <= n && n <= (1i64 << 52), "`n as f64` has rounding error");

    let gram = |x: T| x.rs_theta(atol) - T::PI() * n as f64;
    let t = (T::mp(n as f64) + 0.125) / T::E();
    let w = rgsl::lambert_w::lambert_W0(t.fp());
    let x0 = T::PI() * 2.0 * T::E() * t / w;
    let a = x0 - 1.0;
    let b = x0 + 1.0;
    let result =
        Brentq::new(EvalPoint::new(a, gram), EvalPoint::new(b, gram), atol, 0.0).solve(gram, 100);
    // assert!(result.status == BrentqStatus::Converged, "{:?}, n = {}, atol = {:?}", result, n, atol); // TODO: find good stop criterion
    result.x
}

/// (-1)^n Z(g_n) > 0
#[inline]
fn is_good_gram_point<T: Copy + Signed>(n: i64, g: &EvalPoint<T>) -> bool {
    g.f.is_negative() == (n.rem_euclid(2) == 1)
}

/// Returns whether if g_n is a good Gram point. Note that the
/// Gram points are not accurate: we don't need them to be accurate.
fn next_good_gram_point<T: MyReal + RiemannSiegelTheta>(
    mut n: i64,
    g: EvalPoint<T>,
    mut f: impl FnMut(T) -> T,
) -> Vec<EvalPoint<T>> {
    // assume g_n is a good Gram point
    let mut ret = vec![g];
    loop {
        n += 1;
        let x = T::mp(gram_point::<f64>(n, 1e-12));
        let g = EvalPoint::new(x, &mut f);
        ret.push(g);
        if is_good_gram_point(n, &g) {
            return ret;
        }
    }
}

pub trait HardyZDep =
    MyReal + Sinc + ExpPolyApprox + RiemannSiegelTheta + GabckeSeries + Bernoulli + Factorial;

/// Apply Euler-Maclaurin formula to compute $\zeta(1/2 + i t)$.
/// See Eq (1.2) and (1.3) in Fast Algorithm For Multiple Evaluation Of
/// The Riemann Zeta Function, A.M. Odlyzko and A. Schonhage.
///
/// Although Euler-Maclaurin method requires O(t) time, the $n$ we used is
/// fixed as we only apply Euler-Maclaurin to small $t$ where Riemann-Siegel
/// can't provide enough precision, and larger $n$ always works better.

pub struct EulerMaclaurinMethod<T: MyReal> {
    dirichlet: BandwidthInterp<T>,
    n: usize,
    atol: f64,
}

impl<T: HardyZDep> EulerMaclaurinMethod<T> {
    pub fn new(n: usize, atol: f64) -> Self {
        let dirichlet =
            BandwidthInterp::<T>::new(n - 1, T::zero(), T::mp(n as f64), T::mp(0.5), atol);
        Self { dirichlet, n, atol }
    }

    #[inline(never)]
    pub fn query(&self, x: T) -> T {
        let atol = self.atol;

        let n = self.n; // or another?
        let s = Complex::new(T::mp(0.5), x);
        let n_t = T::mp(n as f64);
        let ln_n = n_t.ln();
        let n_pow_minus_s = (-s * ln_n).exp();
        let mut zeta =
            self.dirichlet.query(x) + n_pow_minus_s * (T::mp(0.5) + n_t / (s - T::one()));

        let mut term = n_pow_minus_s / n_t * s;
        let n_sqr = T::mp(n as f64 * n as f64);
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
pub struct HardyZ<T: HardyZDep> {
    pub dirichlet: Vec<Option<BandwidthInterp<T>>>,
    pub eps: f64,
    pub levels: Vec<(f64, usize)>,
    pub euler_maclaurin: EulerMaclaurinMethod<T>,

    // for memory management
    timestep: usize,
    visit: Vec<usize>,
    cur_mem: usize,
    max_mem: usize,
}

impl<T: HardyZDep> HardyZ<T> {
    pub fn new(max_height: f64, max_order: usize, eps: f64) -> Self {
        let n = (max_height / 2.0 / f64::PI()).sqrt().ceil() as usize + 10;
        let dirichlet = (0..=n).map(|_| None).collect::<Vec<_>>();
        let (levels, em_height) = Self::get_plans(max_height, max_order, eps);
        debug!("levels = {:?}, Euler-Maclaurin height = {:.3e}", levels, em_height);
        let euler_maclaurin = EulerMaclaurinMethod::new(em_height as usize, eps);

        let max_mem = n * 10; // default: we save 10 sets of precomputed results.
        Self {
            dirichlet,
            eps,
            levels,
            euler_maclaurin,
            timestep: 0,
            visit: vec![0; n],
            cur_mem: 0,
            max_mem,
        }
    }

    fn get_plans(max_height: f64, max_order: usize, atol: f64) -> (Vec<(f64, usize)>, f64) {
        const COEFFS: [f64; 11] =
            [0.127, 0.053, 0.011, 0.031, 0.017, 0.061, 0.661, 9.2, 130.0, 1837.0, 25966.0];
        assert!(
            max_order < COEFFS.len(),
            "don't have enough information about Riemann-Siegel formula"
        );
        let mut levels: Vec<(f64, usize)> = vec![];
        let mut cur_t = max_height;
        for k in 0..=max_order {
            let min_t = (COEFFS[k] / atol).powf(4.0 / (3.0 + k as f64 * 2.0));
            if min_t < 200.0 {
                // Riemann-Siegel formula doesn't work for too small $t$.
                // Use Euler-Maclaurin instead.
                break;
            }
            if min_t < cur_t {
                levels.push((min_t, k));
                cur_t = min_t;
            }
        }
        (levels, cur_t)
    }

    fn get(&mut self, n: usize) -> &BandwidthInterp<T> {
        self.timestep += 1;
        self.visit[n] = self.timestep;

        if self.dirichlet[n].is_none() {
            let min_t = T::mp(f64::PI() * 2.0 * (n as f64).powi(2));
            let max_t = T::mp(f64::PI() * 2.0 * (n as f64 + 1.0).powi(2));
            self.dirichlet[n] = Some(BandwidthInterp::new(n, min_t, max_t, T::mp(0.5), self.eps));
            self.cur_mem += n;
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

    pub fn query(&mut self, x: T) -> T {
        let x_fp = x.fp();
        for &(min_t, order) in &self.levels {
            if x_fp >= min_t {
                return self.riemann_siegel(x, order);
            }
        }
        self.euler_maclaurin.query(x)
    }

    pub fn purge(&mut self) {
        if self.cur_mem > self.max_mem {
            let mut indices = (0..self.visit.len()).collect::<Vec<_>>();
            indices.sort_unstable_by_key(|&x| self.visit[x]);
            let goal = self.max_mem / 2;

            let mut released = vec![];
            for i in indices {
                if self.cur_mem <= goal {
                    break;
                }
                if self.dirichlet[i].is_some() {
                    self.dirichlet[i] = None;
                    self.cur_mem -= i;
                    released.push(i);
                }
            }
            debug!("[HardyZ] drop cache: {:?}", released);
        }
    }
}

pub struct HybridPrecHardyZ<T: HardyZDep> {
    lo: HardyZ<f64>,
    hi: HardyZ<T>,
    lo_eps: f64,
}

impl<T: HardyZDep> HybridPrecHardyZ<T> {
    pub fn new(max_height: f64, max_order: usize, eps: f64) -> Self {
        // TODO: need analysis
        let lo_eps = 1e-3;
        let lo = HardyZ::<f64>::new(max_height, max_order, lo_eps);
        let hi = HardyZ::<T>::new(max_height, max_order, eps);
        Self { lo, hi, lo_eps }
    }

    pub fn query(&mut self, x: T) -> T {
        let result = self.lo.query(x.fp());
        if result.abs() > self.lo_eps {
            T::mp(result)
        } else {
            self.hi.query(x)
        }
    }

    pub fn purge(&mut self) {
        self.lo.purge();
        self.hi.purge();
    }
}

#[derive(Default)]
pub struct IsolationStats<T> {
    pub count_locate: usize,
    pub count_separate: usize,
    pub roots: Vec<T>,
    pub gram_start: (i64, f64),
    pub gram_end: (i64, f64),
}

/// Given the current height, determine the number of good Rosser blocks to
/// ensure that all roots are accounted. More specifically, Theorem 3.2 in
/// [Brent79] and Corollary 2.3 in [Trudgan09] provides a lower bound and upper
/// bound for some $g_k$. By applying the theorem twices and trying to close the
/// gap between the lower bound and the upper bound, we can conclude that $g_k =
/// k + 1$.
fn calc_n_good_rosser_blocks_bound(t: f64) -> usize {
    let log_t = t.ln();
    let n1 = 0.0061 * log_t * log_t + 0.08 * log_t; // Theorem 3.2 in [Brent79]
    let n2 = 0.0031 * log_t * log_t + 0.11 * log_t; // Corollary 2.3 in [Trudgan09]
    n1.min(n2).ceil() as usize
}

macro_rules! make_closure {
    ($z: ident, $stats: ident, $mode: ident) => {
        |x| {
            $stats.$mode += 1;
            $z.query(x)
        }
    };
}

/// We first find a Rosser block and try to isolate the number of zeros with
/// moderate efforts. Then a block is put into a ``pending'' zone. We check the
/// pending zone frequently, and apply Turing's method to determine if all roots
/// are accounted. Once a Rosser leaves the pending zone, we have to make sure
/// that all roots inside the Rosser block have been identified.
///
/// The function finds all zeros in $[L, R]$ where $R \geq goal_height$ and $L$
/// is slightly larger than g(n).
///
/// Requirement: $N(g_{n_0}) = n_0 + 1.$
pub fn try_isolate<T: HardyZDep>(
    hardy_z: &mut HybridPrecHardyZ<T>,
    mut n: i64,
    goal_height: f64,
    atol: f64,
    rtol: f64,
    mut certified_n: bool, // whether g(n0) is a regular Gram point.
) -> IsolationStats<T> {
    let init_budget = 1000;
    let huge_budget = 10000;
    let n_iters_to_locate = 100;
    let mut stats = IsolationStats::default();

    if n == -1 && !certified_n {
        log::warn!("n = -1 can be certified");
        certified_n = true;
    }

    let mut g;

    if certified_n {
        g = EvalPoint::new(gram_point(n, atol), make_closure!(hardy_z, stats, count_separate));
        assert!(is_good_gram_point(n, &g));
        stats.gram_start = (n, g.x.fp());
    } else {
        // find the next good Gram point
        let mut block = next_good_gram_point(
            n - 1,
            EvalPoint { x: T::zero(), f: T::zero() },
            make_closure!(hardy_z, stats, count_separate),
        );
        if block.len() != 1 {
            let n0 = n;
            n += (block.len() - 1) as i64;
            debug!("g({n0}) is not good... the next good Gram point is g({})", n);
        }
        g = block.pop().unwrap();
    }

    let mut pending = VecDeque::new();
    let mut n_last_good_rooser_blocks = 0;
    while stats.gram_end.1 < goal_height {
        let block = next_good_gram_point(n, g, make_closure!(hardy_z, stats, count_separate));
        let n_gram_points = block.len() - 1;

        let mut lounge = Lounge::new(n_gram_points, atol, rtol);
        for i in 0..n_gram_points {
            lounge.add(block[i], block[i + 1]);
        }
        lounge.try_isolate(init_budget, make_closure!(hardy_z, stats, count_separate));
        if lounge.separated() {
            n_last_good_rooser_blocks += 1;
        } else {
            n_last_good_rooser_blocks = 0;
        }

        n += n_gram_points as i64;
        g = block[block.len() - 1];

        pending.push_back(lounge);
        let k = calc_n_good_rosser_blocks_bound(g.x.fp());
        if n_last_good_rooser_blocks >= 2 * k && pending.len() > k {
            // ready to identify zeros!

            // This time we have to identify all roots at all costs. The last
            // $k$ Rosser blocks can't be verified by Turing's method.
            let mut lounge = Lounge::new(0, atol, rtol);
            while pending.len() > k {
                lounge.merge(pending.pop_front().unwrap());
            }

            // the last Gram point in `lounge`
            let n_gram_end = n - pending.iter().map(|x| x.n_zeros).sum::<usize>() as i64;

            if !certified_n {
                // OK... We're not sure if all zeros below g(n0) are listed. So
                // we simply ignore all of them and start from the minimum
                // height of $t$ which we're able to certify via Turing's
                // method.
                certified_n = true;
                stats.gram_start = (n_gram_end, lounge.bounds().end);
                debug!(
                    "Gram point g({}) = {:.3} can be certified regular. ",
                    n_gram_end, stats.gram_start.1
                );
                continue;
            }
            stats.gram_end = (n_gram_end, lounge.bounds().end);

            // TODO: check if there is a violation of Rosser's rule here.

            lounge.try_isolate(huge_budget, make_closure!(hardy_z, stats, count_locate));
            assert!(lounge.separated(), "huge budget = {} is not enough!", huge_budget);

            // All roots are separated. Locate the roots with high precision now!
            let mut roots = vec![];
            for bracket in lounge.brackets {
                if let BracketMaintainer::DiffSign(x) = bracket.entry {
                    let result = Brentq::new(x.xa, x.xb, atol, rtol)
                        .solve(make_closure!(hardy_z, stats, count_locate), n_iters_to_locate);
                    assert!(result.status == BrentqStatus::Converged);
                    roots.push(result.x);
                }
            }
            roots.sort_by(|a, b| a.partial_cmp(b).unwrap());
            stats.roots.extend(roots);

            // Clean the memory. Probably we won't need to compute $Z(t)$ for a
            // small $t$ any more.
            hardy_z.purge();
        }
    }
    stats
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
        assert_close(gram_point(-1, eps1), 9.666908056130259, eps2, 0.0);
    }

    #[test]
    fn test_root_find() {
        crate::init();
        type T = f64x2;
        let mut hardy_z = HardyZ::<T>::new(1e8, 10, 1e-18);

        let mut f = |x: T| hardy_z.query(x);
        let center = T::mp(74929.812159119);
        let xa = EvalPoint::new(center - 0.4, &mut f);
        let xb = EvalPoint::new(center + 0.3, &mut f);

        let rtol = 1e-18;
        let result = Brentq::new(xa, xb, 0.0, rtol).solve(f, 100);
        println!("status = {:?}, x = {},f = {}", result.status, result.x, result.f);
        let x = result.x;
        assert_close(x, f64x2::new(74929.812159119, -4.328667942455271e-13), 0.0, rtol);

        // let x = f64x2::new(2669.25374488818, 1.3356510478335532e-13);
        // let z = hardy_z.query(x);
        // assert!
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
