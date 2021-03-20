use std::f64::consts::PI;
use crate::{brentq::brentq, context::Context};
use num::Zero;
use log::{debug};

type Float = f64;
type Complex = num::Complex<Float>;

pub trait FnZeta {
    fn zeta(&self, s: Complex, eps: Float) -> Complex;
}

pub struct ZetaGalway<'a> {
    ctx: &'a Context<f64>,
    pub complexity: i64,
}

impl<'a> ZetaGalway<'a> {
    pub fn new(ctx: &'a Context<f64>) -> Self {
        Self { ctx, complexity: 0 }
    }
}

impl ZetaGalway<'_> {
    fn plan_zeta_galway(&self, s: Complex, eps: Float) -> (Float, Float, Complex, Complex, Complex, Float, Float) {

        fn g(z: Complex, s: Complex) -> Complex { z.powc(-s) * (z * z * Complex::new(0.0, PI)).exp() }
        fn H(w: Complex) -> Complex { 1.0 / (1.0 - (Complex::new(0.0, 2.0 * PI) * w).exp()) }
        fn f(z: Complex, s: Complex) -> Complex {
            let a = (z * Complex::new(0.0, PI)).exp();
            g(z, s) / (a - a.inv())
        }
        fn ln_g_norm(z: Complex, s: Complex) -> f64 { -(s * z.ln()).re - 2.0 * z.im * z.re * PI }
        fn ln_f_norm(z: Complex, s: Complex) -> f64 {
            if z.im.abs() > 30.0 {
                ln_g_norm(z, s) - PI * z.im.abs()
            } else {
                let a = (z * Complex::new(0.0, PI)).exp();
                ln_g_norm(z, s) - (a - a.inv()).norm().ln()
            }
        }
        fn is_close(a: f64, b: f64) -> bool { ((a - b) / a).abs() < 1e-9 }

        let sigma = s.re;

        // let z_s = (s / Complex::new(0.0, PI * 2.0)).sqrt();
        // let n = (z_s.re - z_s.im).floor().max(1.0);
        let n = (s.im / PI / 2.0).sqrt().floor().max(1.0);

        let delta = ((n.powf(-sigma) / eps).ln() / 2.0 / PI).max(0.0);
        let direction = Complex::new(0.5f64.sqrt(), 0.5f64.sqrt());

        let ln_eps = eps.ln();
        // TODO
        let f_alpha = |alpha: Float| {
            // assert!(
            //     is_close(f(n + 0.5 + direction * alpha * delta, s).norm().ln(), ln_f_norm(n + 0.5 + direction * alpha * delta, s)),
            //     "a = {}, b = {}, z = {}",
            //     f(n + 0.5 + direction * alpha * delta, s).norm(),
            //     ln_f_norm(n + 0.5 + direction * alpha * delta, s),
            //     n + 0.5 + direction * alpha * delta,
            // );
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
            m *= 1.1;
            count += 1;
            if count > 100 {
                debug!("s = {}, eps = {}", s, eps);
                assert!(false);
            }
        }

        let n_l = (z_l.re - z_l.im).ceil().max(1.0);
        let n_r = (z_r.re - z_r.im).floor();
        // self.complexity += (m + n_r - n_l) as i64;
        // if rand::random::<f64>() < 0.01 {
        //     debug!("[ZetaGalway plan] m = {}, n_r - n_l = {}", m, n_r - n_l);
        // }

        (n, m, h, z_1, z_2, n_l, n_r)
    }

    fn I0(&self, s: Complex, eps: Float) -> Complex {
        fn g(z: Complex, s: Complex) -> Complex {
            // essentially the following, but avoid large exponent
            // z.powc(-s) * (z * z * Complex::new(0.0, PI)).exp()
            (-s * z.ln() + z * z * Complex::new(0.0, PI)).exp()
        }
        fn H(w: Complex) -> Complex { 1.0 / (1.0 - (Complex::new(0.0, 2.0 * PI) * w).exp()) }
        fn f(z: Complex, s: Complex) -> Complex {
            let a = (z * Complex::new(0.0, PI)).exp();
            g(z, s) / (a - a.inv())
        }

        let (n, m, h, z_1, z_2, n_l, n_r) = self.plan_zeta_galway(s, eps);
        // debug!("{} {} {} {} {} {} {}", n, m, h, z_1, z_2, n_l, n_r);

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
        // debug!("s0 = {}, s1 = {}, s2 = {}, s3 = {}", s0, s1, s2, s3);
        // debug!("??? {} {}", f(z_1, s), g(z_1, s));
        // dbg!(n, m, n_l, n_r, z_1, z_2, (z_2 - z_1).norm());
        // debug!("n_l = {}, n_r = {}, z_l = {}, z_r = {}", n_l, n_r, z_l, z_r);
        // debug!("ret = {}", s0 + s1 * h - s2 + s3);

        s0 + s1 * h - s2 + s3
    }
}

impl FnZeta for ZetaGalway<'_> {
    fn zeta(&self, s: Complex, eps: Float) -> Complex {
        // zeta Galway
        let log_chi = (s - 0.5) * PI.ln() + self.ctx.loggamma((1.0 - s) / 2.0, eps) - self.ctx.loggamma(s / 2.0, eps);
        let chi = log_chi.exp();
        let zeta = self.I0(s, eps) + chi * self.I0(1.0 - s.conj(), eps).conj();
        debug!("s = {}, chi = {}, zeta = {}, {}", s, chi, zeta, self.I0(1.0 - s.conj(), eps));
        zeta
    }
}
