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
}

impl<'a> ZetaGalway<'a> {
    pub fn new(ctx: &'a Context<f64>) -> Self {
        Self { ctx }
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

        let sigma = s.re;

        // let z_s = (s / Complex::new(0.0, PI * 2.0)).sqrt();
        // let n = (z_s.re - z_s.im).floor().max(1.0);
        let n = (s.im / PI / 2.0).sqrt().floor().max(1.0);

        let delta = ((n.powf(-sigma) / eps).ln() / 2.0 / PI).max(0.0);
        let direction = Complex::new(0.5f64.sqrt(), 0.5f64.sqrt());

        // TODO
        let f_alpha = |alpha: Float| f(n + 0.5 + direction * alpha * delta, s).norm() - eps;

        let alpha_1 = brentq(f_alpha, 0.0, 1.0, 1e-6, 1e-6, 20).unwrap();
        let alpha_2 = brentq(f_alpha, 0.0, -1.0, 1e-6, 1e-6, 20).unwrap();
        let z_1 = n + 0.5 + direction * alpha_1 * delta;
        let z_2 = n + 0.5 + direction * alpha_2 * delta;
        // if sigma > 0.0 {
        //     // debug!("z1 = {}, z2 = {}, alpha_1 = {}, alpha_2 = {}", z_1, z_2, alpha_1, alpha_2);
        //     debug!("n = {}, alpha_1 = {}, alpha_2 = {}, o = {}", n, alpha_1, alpha_2, (s.im / 2.0 / PI).sqrt());
        // }
        // debug!("delta = {}, f_alpha(0) = {}, f_alpha(xx) = {}", delta, f_alpha(0.0), g(n + 0.5 + direction * 1.0 * delta, s));

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
            let err1 = g(z_l, s) * ((z_1 - z_l) / h * Complex::new(0.0, 2.0 * PI)).exp();
            let err2 = g(z_r, s) * ((z_r - z_1) / h * Complex::new(0.0, 2.0 * PI)).exp();
            if err1.norm() <= eps && err2.norm() <= eps { break }
            m *= 1.1;
            count += 1;
            if count > 100 {
                debug!("s = {}, eps = {}", s, eps);
                assert!(false);
            }
        }

        let n_l = (z_l.re - z_l.im).ceil().max(1.0);
        let n_r = (z_r.re - z_r.im).floor();

        (n, m, h, z_1, z_2, n_l, n_r)
    }

    fn I0(&self, s: Complex, eps: Float) -> Complex {
        fn g(z: Complex, s: Complex) -> Complex { z.powc(-s) * (z * z * Complex::new(0.0, PI)).exp() }
        fn H(w: Complex) -> Complex { 1.0 / (1.0 - (Complex::new(0.0, 2.0 * PI) * w).exp()) }
        fn f(z: Complex, s: Complex) -> Complex {
            let a = (z * Complex::new(0.0, PI)).exp();
            g(z, s) / (a - a.inv())
        }

        let (n, m, h, z_1, z_2, n_l, n_r) = self.plan_zeta_galway(s, eps);

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
        // debug!("s = {}, chi = {}, zeta = {}", s, chi, zeta);
        zeta
    }
}
