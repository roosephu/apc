use num::Signed;

use crate::traits::MyReal;

// scipy defaults: xtol = 2e-12, rtol = 8.881784197001252e-16

pub fn brentq<T: MyReal, F>(f: F, xa: T, xb: T, xtol: T, rtol: T, iter: usize) -> Option<T>
where
    F: Fn(T) -> T,
{
    let mut xpre = xa;
    let mut xcur = xb;
    let mut xblk = T::zero();
    let mut fpre = f(xpre);
    let mut fcur = f(xcur);
    let mut fblk = T::zero();
    let mut spre = T::zero();
    let mut scur = T::zero();

    if (fpre * fcur).is_sign_positive() {
        return None;
    }
    if fpre.is_zero() {
        return Some(xpre);
    }
    if fcur.is_zero() {
        return Some(xcur);
    }

    for _ in 0..iter {
        if (fpre * fcur).is_sign_negative() {
            xblk = xpre;
            fblk = fpre;
            spre = xcur - xpre;
            scur = spre;
        }

        if fblk.abs() < fcur.abs() {
            xpre = xcur;
            xcur = xblk;
            xblk = xpre;

            fpre = fcur;
            fcur = fblk;
            fblk = fpre;
        }

        let delta = (xtol + rtol * xcur.abs()) / 2.0;
        let sbis = (xblk - xcur) / 2.0;

        if fcur.is_zero() || sbis.abs() < delta {
            return Some(xcur);
        }
        if spre.abs() > delta && fcur.abs() < fpre.abs() {
            let stry = if xpre == xblk {
                -fcur * (xcur - xpre) / (fcur - fpre)
            } else {
                let dpre = (fpre - fcur) / (xpre - xcur);
                let dblk = (fblk - fcur) / (xblk - xcur);
                -fcur * (fblk * dblk - fpre * dpre) / (dblk * dpre * (fblk - fpre))
            };
            if stry.abs() * 2.0 < spre.abs().min(sbis.abs() * 3.0 - delta) {
                spre = scur;
                scur = stry;
            } else {
                spre = sbis;
                scur = sbis;
            }
        } else {
            spre = sbis;
            scur = sbis;
        }

        xpre = xcur;
        fpre = fcur;

        if scur.abs() > delta {
            xcur += scur;
        } else {
            xcur += if sbis.is_sign_positive() { delta } else { -delta };
        }

        fcur = f(xcur);
    }
    Some(xcur)
}

#[derive(Clone, Copy, Debug)]
pub struct EvalPoint<T: Copy> {
    pub(crate) x: T,
    pub(crate) f: T,
}

impl<T: Copy + Signed> EvalPoint<T> {
    #[inline]
    pub fn sign(&self) -> bool { self.f.is_positive() }

    #[inline]
    pub fn has_same_sign(&self, other: &Self) -> bool { self.sign() == other.sign() }

    #[inline]
    pub fn slope(&self, other: &Self) -> T { (self.f - other.f) / (self.x - other.x) }

    #[inline]
    pub fn new(x: T, f: impl FnOnce(T) -> T) -> Self { Self { x, f: f(x) } }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum BrentqStatus {
    Converged,
    SignError,
    ConvError,
    InProgress,
}

#[derive(Debug, PartialEq, Eq)]
pub struct BrentqResult<T: Copy> {
    pub status: BrentqStatus,
    pub n_iters: usize,
    pub x: T,
    pub f: T,
}

impl<T: Copy> BrentqResult<T> {
    pub fn result(&self) -> EvalPoint<T> { EvalPoint { x: self.x, f: self.f } }
}

#[inline(never)]
pub fn brentq2<T: MyReal>(
    mut f: impl FnMut(T) -> T,
    a: EvalPoint<T>,
    b: EvalPoint<T>,
    xtol: f64,
    rtol: f64,
    iter: usize,
) -> BrentqResult<T> {
    let mut pre = a;
    let mut cur = b;
    let mut blk = EvalPoint { x: T::zero(), f: T::zero() };
    let mut spre = T::zero();
    let mut scur = T::zero();

    if pre.has_same_sign(&cur) {
        return BrentqResult { status: BrentqStatus::SignError, n_iters: 0, x: cur.x, f: cur.f };
    }
    if pre.f.is_zero() {
        return BrentqResult { status: BrentqStatus::Converged, n_iters: 0, x: pre.x, f: pre.f };
    }
    if cur.f.is_zero() {
        return BrentqResult { status: BrentqStatus::Converged, n_iters: 0, x: cur.x, f: cur.f };
    }

    for t in 1..=iter {
        if !pre.has_same_sign(&cur) {
            blk = pre;
            spre = cur.x - pre.x;
            scur = spre;
        }

        if blk.f.abs() < cur.f.abs() {
            pre = cur;
            cur = blk;
            blk = pre;
        }

        let delta = (cur.x.abs() * rtol + xtol) * 0.5;
        let sbis = (blk.x - cur.x) * 0.5;

        if cur.f.is_zero() || sbis.abs() < delta {
            return BrentqResult {
                status: BrentqStatus::Converged,
                n_iters: t,
                x: cur.x,
                f: cur.f,
            };
        }
        if spre.abs() > delta && cur.f.abs() < pre.f.abs() {
            let stry = if pre.x == blk.x {
                -cur.f * (cur.x - pre.x) / (cur.f - pre.f)
            } else {
                let dpre = cur.slope(&pre);
                let dblk = cur.slope(&blk);
                -cur.f * (blk.f * dblk - pre.f * dpre) / (dblk * dpre * (blk.f - pre.f))
            };
            if stry.abs() * 2.0 < spre.abs().min(sbis.abs() * 3.0 - delta) {
                spre = scur;
                scur = stry;
            } else {
                spre = sbis;
                scur = sbis;
            }
        } else {
            spre = sbis;
            scur = sbis;
        }

        pre = cur;

        if scur.abs() > delta {
            cur.x += scur;
        } else {
            cur.x += if sbis.is_sign_positive() { delta } else { -delta };
        }

        cur.f = f(cur.x);
    }
    BrentqResult { status: BrentqStatus::ConvError, n_iters: iter, x: cur.x, f: cur.f }
}

#[derive(Clone, Copy)]
pub struct Brentq<T: MyReal> {
    pub status: BrentqStatus,
    xtol: f64,
    rtol: f64,
    n_iters: usize,
    pub(crate) pre: EvalPoint<T>,
    pub(crate) cur: EvalPoint<T>,
    blk: EvalPoint<T>,
    spre: T,
    scur: T,
}

impl<T: MyReal> Brentq<T> {
    pub fn new(a: EvalPoint<T>, b: EvalPoint<T>, xtol: f64, rtol: f64) -> Self {
        let default = Self {
            status: BrentqStatus::InProgress,
            n_iters: 0,
            cur: b,
            pre: a,
            blk: a, // to always maintain a valid bracket [cur, blk]
            spre: T::zero(),
            scur: T::zero(),
            xtol,
            rtol,
        };
        if a.has_same_sign(&b) {
            Self { status: BrentqStatus::SignError, ..default }
        } else {
            default
        }
    }

    // return whether we should stop
    pub fn step(&mut self, f: impl FnOnce(T) -> T) -> bool {
        if self.status != BrentqStatus::InProgress {
            return true;
        }
        self.n_iters += 1;
        if !self.pre.has_same_sign(&self.cur) {
            self.blk = self.pre;
            self.spre = self.cur.x - self.pre.x;
            self.scur = self.spre;
        }

        if self.blk.f.abs() < self.cur.f.abs() {
            self.pre = self.cur;
            self.cur = self.blk;
            self.blk = self.pre;
        }

        let delta = (self.cur.x.abs() * self.rtol + self.xtol) * 0.5;
        let sbis = (self.blk.x - self.cur.x) * 0.5;

        if self.cur.f.is_zero() || sbis.abs() < delta {
            self.status = BrentqStatus::Converged;
            return true;
        }

        if self.spre.abs() > delta && self.cur.f.abs() < self.pre.f.abs() {
            let stry = if self.pre.x == self.blk.x {
                -self.cur.f * (self.cur.x - self.pre.x) / (self.cur.f - self.pre.f)
            } else {
                let dpre = self.cur.slope(&self.pre);
                let dblk = self.cur.slope(&self.blk);
                -self.cur.f * (self.blk.f * dblk - self.pre.f * dpre)
                    / (dblk * dpre * (self.blk.f - self.pre.f))
            };
            if stry.abs() * 2.0 < self.spre.abs().min(sbis.abs() * 3.0 - delta) {
                self.spre = self.scur;
                self.scur = stry;
            } else {
                self.spre = sbis;
                self.scur = sbis;
            }
        } else {
            self.spre = sbis;
            self.scur = sbis;
        }

        self.pre = self.cur;

        if self.scur.abs() > delta {
            self.cur.x += self.scur;
        } else {
            self.cur.x += if sbis.is_sign_positive() { delta } else { -delta };
        }

        self.cur.f = f(self.cur.x);
        false
    }

    pub fn result(&self) -> BrentqResult<T> {
        BrentqResult { status: self.status, n_iters: self.n_iters, x: self.cur.x, f: self.cur.f }
    }

    #[inline(never)]
    pub fn solve(&mut self, mut f: impl FnMut(T) -> T, n_iters: usize) -> BrentqResult<T> {
        for _ in 0..n_iters {
            if self.step(&mut f) {
                break;
            }
        }
        if self.status == BrentqStatus::InProgress {
            self.status = BrentqStatus::ConvError;
        }
        self.result()
    }
}

#[cfg(test)]
mod tests {
    use super::brentq;

    #[test]
    fn test_brentq() {
        let f = |x| x * x * x - 0.5;
        let result = brentq(f, 0.0, 1.0, 1e-8, 1e-8, 100);
        let solution = 0.5f64.powf(1.0 / 3.0);
        assert!((result.unwrap() - solution).abs() <= 1e-8);
    }
}
