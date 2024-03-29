use crate::{brentq::brentq, context::Context, traits::MyReal};
use crate::{
    sum_trunc_dirichlet::sum_trunc_dirichlet,
    traits::{ComplexFunctions, ExpPolyApprox},
};
use core::panic;
use log::debug;
use num::Complex;
use rustfft::FftNum;
use std::f64::consts::PI;

type Float = f64;
type Complex64 = num::Complex<Float>;

pub trait FnZeta<T> {
    fn zeta(&mut self, s: Complex<T>, eps: f64) -> Complex<T>;
    fn prepare_multi_eval(&mut self, h: T, eps: f64);
}

/// essentially the following, but avoid large exponent
/// z.powc(-s) * (z * z * Complex::new(0.0, PI)).exp()
fn g<T: MyReal>(z: Complex<T>, s: Complex<T>) -> Complex<T> {
    (-s * z.ln() + (z * z * T::PI()).mul_i()).exp()
}

/// z = O(n) + i O(ln(n)), s = O(1) + i t.
/// thus z.ln() = O(ln(n)) + i O(ln(n) / n)
/// maximum number in the computation is O(t ln(n) / n) = O(n ln(n)).
/// typically we return a number in O(n), so it's fine to use f64.
fn ln_g_norm(z: Complex<f64>, s: Complex<f64>) -> f64 {
    -(s * z.ln()).re - z.im * z.re * f64::TAU()
}

/// z = O(n) + i O(ln(n)), s = O(1) + i t.
fn f<T: MyReal>(z: Complex<T>, s: Complex<T>) -> Complex<T> {
    let z_im: f64 = z.im.fp();
    if z_im > 100.0 {
        let a = z.scale(T::PI()).mul_i().exp();
        g(z, s) / a
    } else if z_im < -100.0 {
        let a = z.scale(-T::PI()).mul_i().exp();
        g(z, s) / (-a)
    } else {
        let a = z.scale(T::PI()).mul_i().exp();
        g(z, s) / (a - a.inv())
    }
}

/// z = O(n) + i O(ln(n)), s = O(1) + i t.
fn ln_f_norm(z: Complex<f64>, s: Complex<f64>) -> f64 {
    if z.im.abs() > 30.0 {
        ln_g_norm(z, s) - f64::PI() * z.im.abs()
    } else {
        let a = z.mul_i().exp();
        ln_g_norm(z, s) - (a - a.inv()).norm().ln()
    }
}

fn H<T: MyReal>(w: Complex<T>) -> Complex<T> {
    if w.im < T::mp(-30.0f64) {
        Complex::<T>::zero()
    } else {
        T::one() / (T::one() - w.mul_i().scale(T::TAU()).exp())
    }
}

fn is_close(a: f64, b: f64) -> bool { ((a - b) / a).abs() < 1e-9 }

type Plan<T> = (i64, i64, i64, i64, Complex<f64>, Complex<f64>, Option<Complex<T>>);

#[derive(Default)]
struct ZetaGalwayPlanner<T: MyReal + ExpPolyApprox + FftNum> {
    eps: f64,
    ln_eps: f64,
    replan: bool,
    sigma: f64,
    delta: f64,
    h: Option<T>,
    s0: Complex<T>,
    sum_trunc_dirichlet: Vec<Complex<T>>,

    n: f64,
    m: f64,
    alpha_1: f64,
    alpha_2: f64,
}

impl<T: MyReal + ExpPolyApprox + FftNum> ZetaGalwayPlanner<T> {
    pub fn new() -> Self { Self { replan: true, ..Default::default() } }

    fn plan_sum_trunc_dirichlet(&mut self, s: Complex<T>, n: f64) {
        let mut m = 0;
        let h = self.h.unwrap();
        loop {
            let s_approx = Complex::new(s.re, s.im + T::from_usize(m).unwrap() * h).approx();
            let z_s = (s_approx / Complex64::new(0.0, PI * 2.0)).sqrt();
            let n2 = (z_s.re - z_s.im).floor().max(1.0) as f64;
            if n2 != n {
                break;
            }
            m += 1;
        }

        // add m by 1 to avoid offset-by-1 due to limited precision.
        self.sum_trunc_dirichlet = sum_trunc_dirichlet(s, 1, n as usize, m + 1, h);
        self.s0 = s;
    }

    fn test_m(
        m: f64,
        z_1: Complex<f64>,
        z_2: Complex<f64>,
        s: Complex<f64>,
        ln_eps: f64,
    ) -> Option<(Complex<f64>, i64, i64)> {
        let h = (z_2 - z_1) / m;
        let inv_2h = 0.5 / h;
        let z_l = (inv_2h * inv_2h - s.mul_i() / f64::TAU()).sqrt() + inv_2h;
        let z_r = (inv_2h * inv_2h - s.mul_i() / f64::TAU()).sqrt() - inv_2h;
        let ln_err1 = ln_g_norm(z_l, s) - ((z_1 - z_l) / h).im * f64::TAU();
        let ln_err2 = ln_g_norm(z_r, s) - ((z_r - z_1) / h).im * f64::TAU();
        if ln_err1 <= ln_eps && ln_err2 <= ln_eps {
            let n_l = (z_l.re - z_l.im).ceil().max(1.0) as i64;
            let n_r = (z_r.re - z_r.im).floor() as i64;
            Some((h, n_l, n_r))
        } else {
            None
        }
    }

    fn plan_from_scratch(&mut self, s: Complex<f64>, eps: f64) -> Plan<T> {
        let n = self.n;
        let sigma = s.re;
        let ln_eps = eps.ln();

        let delta = ((n.powf(-sigma) / eps).ln() / 2.0 / PI).max(0.0);
        let direction = Complex64::new(0.5f64.sqrt(), 0.5f64.sqrt());

        let f_alpha = |alpha: Float| ln_f_norm(n + 0.5 + direction * alpha * delta, s) - ln_eps;
        // println!("{} {} {}", f_alpha(0.0), f_alpha(2.0), f_alpha(-2.0));

        let alpha_1 = brentq(f_alpha, 0.0, 2.0, 1e-6, 1e-6, 20).unwrap();
        let alpha_2 = brentq(f_alpha, 0.0, -2.0, 1e-6, 1e-6, 20).unwrap();
        let z_1 = n + 0.5 + direction * alpha_1 * delta;
        let z_2 = n + 0.5 + direction * alpha_2 * delta;

        let mut m = 1.0;
        let mut result = None;

        for count in 0..=100 {
            if count == 100 {
                debug!("s = {}, eps = {}", s, eps);
                panic!();
            }

            result = Self::test_m(m, z_1, z_2, s, ln_eps);
            if result.is_some() {
                break;
            }
            m = (m * 1.1).ceil();
        }

        let (h, n_l, n_r) = result.unwrap();

        self.sigma = sigma;
        self.eps = eps;
        self.ln_eps = ln_eps;
        self.delta = delta;

        self.alpha_1 = alpha_1;
        self.alpha_2 = alpha_2;
        self.m = m;
        // debug!(
        //     "[Zeta plan from scratch] alpha_1 = {:.6}, alpha_2 = {:.6}",
        //     self.alpha_1, self.alpha_2
        // );
        // assert!((0.2..=0.8).contains(&alpha_1) && (-0.8..=-0.2).contains(&alpha_2));

        (n as i64, m as i64, n_l, n_r, h, z_1, None)
    }

    fn plan_incremental(&mut self, s: Complex64) -> Plan<T> {
        let n = self.n;
        let ln_eps = self.ln_eps;

        let delta = self.delta;
        let direction = Complex64::new(0.5f64.sqrt(), 0.5f64.sqrt());

        let f_alpha = |alpha: Float| ln_f_norm(n + 0.5 + direction * alpha * delta, s) - ln_eps;

        const ALPHA_DIFF: f64 = 0.03;
        let alpha_1 =
            brentq(f_alpha, self.alpha_1 - ALPHA_DIFF, self.alpha_1 + ALPHA_DIFF, 1e-6, 1e-6, 3)
                .unwrap();
        let alpha_2 =
            brentq(f_alpha, self.alpha_2 - ALPHA_DIFF, self.alpha_2 + ALPHA_DIFF, 1e-6, 1e-6, 3)
                .unwrap();
        let z_1 = n + 0.5 + direction * alpha_1 * delta;
        let z_2 = n + 0.5 + direction * alpha_2 * delta;

        let mut m = self.m;
        let mut result;

        loop {
            result = Self::test_m(m, z_1, z_2, s, ln_eps);
            if result.is_some() {
                break;
            }
            m += 1.0;
        }
        while m > 1.0 {
            let result2 = Self::test_m(m - 1.0, z_1, z_2, s, ln_eps);
            if result2.is_none() {
                break;
            }
            m -= 1.0;
            result = result2;
        }
        let (h, n_l, n_r) = result.unwrap();

        if rand::random::<f64>() < 0.001 {
            // debug!("[I] a1 = {:.6}, a2 = {:.6}, delta = {:.6}", alpha_1, alpha_2, delta);
            // debug!("m = {}", m);
        }
        // assert!((0.1..=0.9).contains(&alpha_1) && (-0.8..=-0.2).contains(&alpha_2));

        self.m = m;
        self.alpha_1 = alpha_1;
        self.alpha_2 = alpha_2;
        (n as i64, m as i64, n_l, n_r, h, z_1, None)
    }

    pub fn plan(&mut self, s: Complex<T>, eps: f64) -> Plan<T> {
        let old_n = self.n;

        let s_approx = s.approx();
        let z_s = (s_approx / Complex64::new(0.0, PI * 2.0)).sqrt();
        // can simply take sqrt(t).
        // let n = (s.im / PI / 2.0).sqrt().floor().max(1.0);
        self.n = (z_s.re - z_s.im).floor().max(1.0);

        let mut plan =
            if self.replan || old_n != self.n || s_approx.re != self.sigma || eps != self.eps {
                self.replan = false;
                if self.h.is_some() {
                    self.plan_sum_trunc_dirichlet(s, self.n);
                }
                self.plan_from_scratch(s_approx, eps)
            } else {
                self.plan_incremental(s_approx)
            };
        if let Some(h) = self.h {
            let idx = (s.im - self.s0.im) / h;
            let round_idx = idx.round();
            let int_idx = round_idx.to_usize().unwrap();
            assert!((idx - round_idx).abs().fp() <= 1e-8);

            plan.6 = Some(self.sum_trunc_dirichlet[int_idx]);
        }
        plan
    }
}

pub struct ZetaGalway<'a, T: MyReal> {
    ctx: &'a Context<T>,
    planners: [ZetaGalwayPlanner<T>; 2],
    pub complexity: i64,
}

impl<'a, T: MyReal> ZetaGalway<'a, T> {
    pub fn new(ctx: &'a Context<T>) -> Self {
        Self { ctx, planners: [ZetaGalwayPlanner::new(), ZetaGalwayPlanner::new()], complexity: 0 }
    }
}

impl<T: MyReal> ZetaGalway<'_, T> {
    fn I0(&self, s: Complex<T>, plan: Plan<T>) -> Complex<T> {
        let (n, m, n_l, n_r, h, z_1, plan_sum_trunc_dirichlet) = plan;
        // we don't care the precise value of z_1 and h, as it's only for correction.
        let z_1 = Complex::<T>::new(T::mp(z_1.re), T::mp(z_1.im));
        let h = Complex::<T>::new(T::mp(h.re), T::mp(h.im));

        let mut s0;
        match plan_sum_trunc_dirichlet {
            Some(x) => {
                s0 = x;
            }
            None => {
                s0 = Complex::<T>::zero();
                for i in 1..=n {
                    s0 += (-s * T::from_i64(i).unwrap().ln()).exp();
                }
                // panic!("why don't you use [FKBJ-OS]?");
            }
        }

        let mut s1 = Complex::<T>::zero();
        for i in 0..=m {
            s1 += f(z_1 + h * T::from_i64(i).unwrap(), s);
            // debug!("{:?}, {:?} {:?}", z_1 + h * i.unwrap(), s, f(z_1 + h * i.unwrap(), s));
        }

        let mut s2 = Complex::<T>::zero();
        for i in n_l..=n {
            let cast_i = T::from_i64(i).unwrap();
            s2 += Complex::<T>::new(cast_i, T::zero()).powc(-s) * H((cast_i - z_1) / h);
        }

        let mut s3 = Complex::<T>::zero();
        for i in n + 1..=n_r {
            let cast_i = T::from_i64(i).unwrap();
            s3 += Complex::<T>::new(cast_i, T::zero()).powc(-s) * H((z_1 - cast_i) / h);
            // println!("{:?} {:?}", (z_1 - cast_i) / h, H((z_1 - cast_i) / h));
        }
        // debug!("s = {}, plan = {:?}", s, plan);
        // debug!("s0 = {}, s1 = {}, s2 = {}, s3 = {}", s0, s1, s2, s3);
        // info!("ret = {}", s0 + s1 * h - s2 + s3);

        s0 + s1 * h - s2 + s3
    }

    fn test(&mut self, s: Complex<T>, eps: f64) {
        let log_chi = (s - T::mp(0.5)) * T::PI().ln()
            + self.ctx.loggamma((T::one() - s) * T::mp(0.5), eps)
            - self.ctx.loggamma(s * T::mp(0.5), eps);
        let a = self.ctx.loggamma((T::one() - s) * T::mp(0.5), eps);
        let b = self.ctx.loggamma(s * T::mp(0.5), eps);
        println!("a = {:?}, b = {:?}", a, b);
        println!("log chi = {:?}", log_chi);
    }
}

impl<T: MyReal> FnZeta<T> for ZetaGalway<'_, T> {
    fn zeta(&mut self, s: Complex<T>, eps: f64) -> Complex<T> {
        // println!("??? s = {}, eps = {}", s, eps);
        let log_chi = (s - T::mp(0.5)) * T::PI().ln()
            + self.ctx.loggamma((T::one() - s) * T::mp(0.5), eps)
            - self.ctx.loggamma(s * T::mp(0.5), eps);
        assert!(!log_chi.re.is_nan(), "{:?} {}", s, eps);
        let chi = log_chi.exp();

        let plan0 = self.planners[0].plan(s, eps);
        let plan1 = self.planners[1].plan(T::one() - s.conj(), eps);
        let ret = self.I0(s, plan0) + chi * self.I0(T::one() - s.conj(), plan1).conj();
        // println!("!!! ret = {}, {}, {}, {}", ret, self.I0(s, plan0), chi, self.I0(T::one() - s.conj(), plan1));
        ret
    }

    fn prepare_multi_eval(&mut self, h: T, _: f64) {
        self.planners[0].h = Some(h);
        self.planners[1].h = Some(h);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use F64x2::f64x2;

    #[test]
    fn test_zeta_galway() {
        let ctx2 = Context::<f64x2>::new(100);
        let mut zeta_galway2 = ZetaGalway::new(&ctx2);

        let s = Complex { re: f64x2 { hi: 1.5, lo: 0.0 }, im: f64x2 { hi: 10000.0, lo: 0.0 } };
        let eps = 1e-18;
        println!("s = {}, zeta(s) = {}", s, zeta_galway2.zeta(s, eps));
        // panic!();

        // zeta_galway2.test(s, eps);
    }
}
