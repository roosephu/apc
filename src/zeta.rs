use std::f64::consts::PI;
use crate::{brentq::brentq, context::Context, traits::GenericFloat};
use log::{debug, info};
use crate::traits::Rotate90;

type Float = f64;
type Complex = num::Complex<Float>;

pub trait FnZeta {
    fn zeta(&mut self, s: Complex, eps: Float) -> Complex;
}

fn g(z: Complex, s: Complex) -> Complex {
    // essentially the following, but avoid large exponent
    // z.powc(-s) * (z * z * Complex::new(0.0, PI)).exp()
    (-s * z.ln() + (z * z * PI).rotate90()).exp()
}

fn f(z: Complex, s: Complex) -> Complex {
    let a = z.scale(PI).rotate90().exp();
    g(z, s) / (a - a.inv())
}

fn H(w: Complex) -> Complex {
    1.0 / (1.0 - w.rotate90().scale(2.0 * PI).exp())
}

fn ln_g_norm(z: Complex, s: Complex) -> f64 { -(s * z.ln()).re - 2.0 * z.im * z.re * PI }
fn ln_f_norm(z: Complex, s: Complex) -> f64 {
    if z.im.abs() > 30.0 {
        ln_g_norm(z, s) - PI * z.im.abs()
    } else {
        let a = z.rotate90().exp();
        ln_g_norm(z, s) - (a - a.inv()).norm().ln()
    }
}

fn is_close(a: f64, b: f64) -> bool { ((a - b) / a).abs() < 1e-9 }

#[derive(Default)]
struct ZetaGalwayPlanner {
    eps: f64,
    ln_eps: f64,
    replan: bool,
    sigma: f64,
    n: f64,
    n_l: f64,
    n_r: f64,
    m: f64,
    alpha_1: f64,
    alpha_2: f64,

    h: Complex,
    z_1: Complex,
    z_2: Complex,
}

impl ZetaGalwayPlanner {
    pub fn new() -> Self {
        Self {
            replan: true,
            ..Default::default()
        }
    }

    fn plan_from_scratch(&mut self, s: Complex, eps: f64)  {
        let n = self.n;
        let sigma = s.re;
        let ln_eps = eps.ln();

        let delta = ((n.powf(-sigma) / eps).ln() / 2.0 / PI).max(0.0);
        let direction = Complex::new(0.5f64.sqrt(), 0.5f64.sqrt());

        // TODO
        let f_alpha = |alpha: Float| {
            ln_f_norm(n + 0.5 + direction * alpha * delta, s) - ln_eps
        };

        let alpha_1 = brentq(f_alpha, 0.0, 1.0, 1e-6, 1e-6, 10).unwrap();
        let alpha_2 = brentq(f_alpha, 0.0, -1.0, 1e-6, 1e-6, 10).unwrap();
        let z_1 = n + 0.5 + direction * alpha_1 * delta;
        let z_2 = n + 0.5 + direction * alpha_2 * delta;

        fn test_g_norm(z: Complex, s: Complex) {
            let ans = g(z, s).norm().ln();
            let value = ln_g_norm(z, s);
            assert!((ans - value) / ans.abs() < 1e-9, "ans = {}, value = {}", ans, value);
        }

        let mut m = 1.0;
        let mut h;
        let mut z_l;
        let mut z_r;
        let mut count = 0;

        loop {
            h = (z_2 - z_1) / m;
            let inv_2h = 0.5 / h;
            z_l = (inv_2h * inv_2h + s / Complex::new(0.0, 2.0 * PI)).sqrt() + inv_2h;
            z_r = (inv_2h * inv_2h + s / Complex::new(0.0, 2.0 * PI)).sqrt() - inv_2h;
            let ln_err1 = ln_g_norm(z_l, s) - ((z_1 - z_l) / h).im * 2.0 * PI;
            let ln_err2 = ln_g_norm(z_r, s) - ((z_r - z_1) / h).im * 2.0 * PI;
            if ln_err1 <= ln_eps && ln_err2 <= ln_eps { break }
            m = (m * 1.1).ceil();
            count += 1;
            if count > 100 {
                debug!("s = {}, eps = {}", s, eps);
                panic!();
            }
        }

        let n_l = (z_l.re - z_l.im).ceil().max(1.0);
        let n_r = (z_r.re - z_r.im).floor();

        self.sigma = sigma;
        self.eps = eps;
        self.ln_eps = ln_eps;
        self.n_l = n_l;
        self.n_r = n_r;
        self.m = m;
        self.alpha_1 = alpha_1;
        self.alpha_2 = alpha_2;
        self.z_1 = z_1;
        self.z_2 = z_2;
        self.h = h;
    }

    fn plan_incremental(&mut self, s: Complex) {
        unimplemented!()
    }

    pub fn plan(&mut self, s: Complex, eps: f64) {
        let old_n = self.n;

        let z_s = (s / Complex::new(0.0, PI * 2.0)).sqrt();
        // can simply take sqrt(t).
        // let n = (s.im / PI / 2.0).sqrt().floor().max(1.0);
        self.n = (z_s.re - z_s.im).floor().max(1.0);

        if self.replan || old_n != self.n || s.re != self.sigma || eps != self.eps {
            self.replan = false;
            self.plan_from_scratch(s, eps);
        } else {
            // self.plan_incremental(s)
            self.plan_from_scratch(s, eps);
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
            complexity: 0
        }
    }
}

impl ZetaGalway<'_> {
    fn I0(&self, s: Complex, planner: &ZetaGalwayPlanner) -> Complex {
        let n = planner.n;
        let m = planner.m;
        let h = planner.h;
        let z_1 = planner.z_1;
        let z_2 = planner.z_2;
        let n_l = planner.n_l;
        let n_r = planner.n_r;

        // info!("n = {}, m = {}, h = {}, n_l = {}, n_r = {}", n, m, h, n_l, n_r);

        let mut s0 = Complex::zero();
        for i in 1..=n as i64 {
            s0 += Complex::new(i as Float, 0.0).powc(-s);
        }

        let mut s1 = Complex::zero();
        for i in 0..=m as i64 {
            s1 += f(z_1 + i as Float * h, s);
        }

        let mut s2 = Complex::zero();
        for i in n_l as i64..=n as i64 {
            s2 += Complex::new(i as Float, 0.0).powc(-s) * H((i as Float - z_1) / h);
        }

        let mut s3 = Complex::zero();
        for i in n as i64 + 1 ..= n_r as i64 {
            s3 += Complex::new(i as Float, 0.0).powc(-s) * H((z_1 - i as Float) / h);
        }
        // info!("s0 = {}, s1 = {}, s2 = {}, s3 = {}", s0, s1, s2, s3);
        // info!("ret = {}", s0 + s1 * h - s2 + s3);

        s0 + s1 * h - s2 + s3
    }
}

impl FnZeta for ZetaGalway<'_> {
    fn zeta(&mut self, s: Complex, eps: f64) -> Complex {
        let log_chi = (s - 0.5) * PI.ln() + self.ctx.loggamma((1.0 - s) / 2.0, eps) - self.ctx.loggamma(s / 2.0, eps);
        let chi = log_chi.exp();

        self.planners[0].plan(s, eps);
        self.planners[1].plan(1.0 - s.conj(), eps);
        self.I0(s, &self.planners[0]) + chi * self.I0(1.0 - s.conj(), &self.planners[1]).conj()
    }
}
