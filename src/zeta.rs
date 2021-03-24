use crate::traits::ComplexFunctions;
use crate::{brentq::brentq, context::Context, traits::GenericFloat};
use log::{debug, info};
use num::Complex;
use std::f64::consts::PI;

type Float = f64;
type Complex64 = num::Complex<Float>;

pub trait FnZeta<T> {
    fn zeta(&mut self, s: Complex<T>, eps: f64) -> Complex<T>;
}

/// essentially the following, but avoid large exponent
/// z.powc(-s) * (z * z * Complex::new(0.0, PI)).exp()
fn g<T: GenericFloat>(z: Complex<T>, s: Complex<T>) -> Complex<T> {
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
fn f<T: GenericFloat>(z: Complex<T>, s: Complex<T>) -> Complex<T> {
    if z.im.as_() > 100.0 {
        let a = z.scale(T::PI()).mul_i().exp();
        g(z, s) / a
    } else if z.im.as_() < -100.0 {
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

fn H<T: GenericFloat>(w: Complex<T>) -> Complex<T> {
    T::one() / (T::one() - w.mul_i().scale(T::TAU()).exp())
}

fn is_close(a: f64, b: f64) -> bool {
    ((a - b) / a).abs() < 1e-9
}

type Plan<T> = (T, T, T, T, Complex<T>, Complex<T>);

#[derive(Default)]
struct ZetaGalwayPlanner {
    eps: f64,
    ln_eps: f64,
    replan: bool,
    sigma: f64,
    delta: f64,

    n: f64,
    m: f64,
    alpha_1: f64,
    alpha_2: f64,
}

impl ZetaGalwayPlanner {
    pub fn new() -> Self {
        Self {
            replan: true,
            ..Default::default()
        }
    }

    fn test_m(
        m: f64,
        z_1: Complex<f64>,
        z_2: Complex<f64>,
        s: Complex<f64>,
        ln_eps: f64,
    ) -> Option<(Complex<f64>, f64, f64)> {
        let h = (z_2 - z_1) / m;
        let inv_2h = 0.5 / h;
        let z_l = (inv_2h * inv_2h - s.mul_i() / f64::TAU()).sqrt() + inv_2h;
        let z_r = (inv_2h * inv_2h - s.mul_i() / f64::TAU()).sqrt() - inv_2h;
        let ln_err1 = ln_g_norm(z_l, s) - ((z_1 - z_l) / h).im * f64::TAU();
        let ln_err2 = ln_g_norm(z_r, s) - ((z_r - z_1) / h).im * f64::TAU();
        if ln_err1 <= ln_eps && ln_err2 <= ln_eps {
            let n_l = (z_l.re - z_l.im).ceil().max(1.0);
            let n_r = (z_r.re - z_r.im).floor();
            Some((h, n_l, n_r))
        } else {
            None
        }
    }

    fn plan_from_scratch(&mut self, s: Complex<f64>, eps: f64) -> Plan<f64> {
        let n = self.n;
        let sigma = s.re;
        let ln_eps = eps.ln();

        let delta = ((n.powf(-sigma) / eps).ln() / 2.0 / PI).max(0.0);
        let direction = Complex64::new(0.5f64.sqrt(), 0.5f64.sqrt());

        let f_alpha = |alpha: Float| ln_f_norm(n + 0.5 + direction * alpha * delta, s) - ln_eps;

        let alpha_1 = brentq(f_alpha, 0.0, 1.0, 1e-6, 1e-6, 10).unwrap();
        let alpha_2 = brentq(f_alpha, 0.0, -1.0, 1e-6, 1e-6, 10).unwrap();
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
        debug!(
            "[Plan from scratch] alpha_1 = {}, alpha_2 = {}",
            self.alpha_1, self.alpha_2
        );
        assert!(0.2 <= alpha_1 && alpha_1 <= 0.8 && -0.8 <= alpha_2 && alpha_2 <= -0.2);

        (n, m, n_l, n_r, h, z_1)
    }

    fn plan_incremental(&mut self, s: Complex64) -> Plan<f64> {
        let eps = self.eps;
        let n = self.n;
        let sigma = s.re;
        let ln_eps = self.ln_eps;

        let delta = self.delta;
        let direction = Complex64::new(0.5f64.sqrt(), 0.5f64.sqrt());

        let f_alpha = |alpha: Float| ln_f_norm(n + 0.5 + direction * alpha * delta, s) - ln_eps;

        const ALPHA_DIFF: f64 = 0.03;
        let alpha_1 = brentq(
            f_alpha,
            self.alpha_1 - ALPHA_DIFF,
            self.alpha_1 + ALPHA_DIFF,
            1e-6,
            1e-6,
            5,
        )
        .unwrap();
        let alpha_2 = brentq(
            f_alpha,
            self.alpha_2 - ALPHA_DIFF,
            self.alpha_2 + ALPHA_DIFF,
            1e-6,
            1e-6,
            5,
        )
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

        if rand::random::<f64>() < 0.01 {
            // debug!("[I] a1 = {:.6}, a2 = {:.6}, delta = {:.6}", alpha_1, alpha_2, delta);
            // debug!("t = {}", s.im);
        }
        assert!(0.2 <= alpha_1 && alpha_1 <= 0.8 && -0.8 <= alpha_2 && alpha_2 <= -0.2);

        self.m = m;
        self.alpha_1 = alpha_1;
        self.alpha_2 = alpha_2;
        (n, m, n_l, n_r, h, z_1)
    }

    pub fn plan(&mut self, s: Complex64, eps: f64) -> Plan<f64> {
        let old_n = self.n;

        let z_s = (s / Complex64::new(0.0, PI * 2.0)).sqrt();
        // can simply take sqrt(t).
        // let n = (s.im / PI / 2.0).sqrt().floor().max(1.0);
        self.n = (z_s.re - z_s.im).floor().max(1.0);

        if self.replan || old_n != self.n || s.re != self.sigma || eps != self.eps {
            self.replan = false;
            self.plan_from_scratch(s, eps)
        } else {
            self.plan_incremental(s)
        }
    }
}

pub struct ZetaGalway<'a> {
    ctx: &'a Context<f64>,
    planners: [ZetaGalwayPlanner; 2],
    pub complexity: i64,
}

impl<'a> ZetaGalway<'a> {
    pub fn new(ctx: &'a Context<f64>) -> Self {
        Self {
            ctx,
            planners: [ZetaGalwayPlanner::new(), ZetaGalwayPlanner::new()],
            complexity: 0,
        }
    }
}

impl ZetaGalway<'_> {
    fn I0(&self, s: Complex64, plan: Plan<f64>) -> Complex64 {
        let (n, m, n_l, n_r, h, z_1) = plan;
        // info!("n = {}, m = {}, h = {}, n_l = {}, n_r = {}", n, m, h, n_l, n_r);

        let mut s0 = Complex64::zero();
        for i in 1..=n as i64 {
            s0 += Complex64::new(i as Float, 0.0).powc(-s);
        }

        let mut s1 = Complex64::zero();
        for i in 0..=m as i64 {
            s1 += f(z_1 + i as Float * h, s);
        }

        let mut s2 = Complex64::zero();
        for i in n_l as i64..=n as i64 {
            s2 += Complex64::new(i as Float, 0.0).powc(-s) * H((i as Float - z_1) / h);
        }

        let mut s3 = Complex64::zero();
        for i in n as i64 + 1..=n_r as i64 {
            s3 += Complex64::new(i as Float, 0.0).powc(-s) * H((z_1 - i as Float) / h);
        }
        // info!("s0 = {}, s1 = {}, s2 = {}, s3 = {}", s0, s1, s2, s3);
        // info!("ret = {}", s0 + s1 * h - s2 + s3);

        s0 + s1 * h - s2 + s3
    }
}

impl FnZeta<f64> for ZetaGalway<'_> {
    fn zeta(&mut self, s: Complex64, eps: f64) -> Complex64 {
        let log_chi = (s - 0.5) * PI.ln() + self.ctx.loggamma((1.0 - s) / 2.0, eps)
            - self.ctx.loggamma(s / 2.0, eps);
        let chi = log_chi.exp();

        let s_approx = Complex::<f64>::new(s.re.as_(), s.im.as_());

        let plan0 = self.planners[0].plan(s_approx, eps);
        let plan1 = self.planners[1].plan(1.0 - s_approx.conj(), eps);
        self.I0(s, plan0) + chi * self.I0(1.0 - s.conj(), plan1).conj()
    }
}
