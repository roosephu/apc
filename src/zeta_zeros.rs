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

/// For root finding
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
        // if n_zeros > 1 {
        //     debug!("initial block {:?}", block);
        // }

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

    #[test]
    fn test_gram() {
        crate::init();
        let eps1 = 1e-11;
        let eps2 = 1e-9;
        assert_close(gram_point(100, eps1), 238.5825905145, eps2);
    }
}
