use std::f64::consts::PI;
use crate::gamma::gamma;
use crate::brentq::brentq;

type Float = f64;
type Complex = num::Complex<Float>;


fn I0(s: Complex, eps: Float) -> Complex {
    fn g(z: Complex, s: Complex) -> Complex {
        z.powc(-s) * (z * z * Complex::i() * PI).exp()
    }

    fn f(z: Complex, s: Complex) -> Complex {
        let a = (z * Complex::i() * PI).exp();
        let denom = a - a.inv();
        g(z, s) / denom
    }

    fn H(w: Complex) -> Complex {
        1.0 / (1.0 - (2.0 * PI * Complex::i() * w).exp())
    }

    let sigma = s.re;

    let z_s = (s / Complex::new(0.0, PI * 2.0)).sqrt();
    let n = (z_s.re - z_s.im).floor();

    let delta = ((n.powf(-sigma) / eps).ln() / 2.0 / PI).max(0.0);
    let direction = Complex::new(0.5f64.sqrt(), 0.5f64.sqrt());

    let f_alpha = |alpha: Float| f(n + 0.5 + direction * alpha * delta, s).norm() - eps;

    let alpha_1 = brentq(f_alpha, 0.0, 2.0, 1e-6, 1e-6, 20);
    let alpha_2 = brentq(f_alpha, 0.0, -2.0, 1e-6, 1e-6, 20);
    let z_1 = n + 0.5 + direction * alpha_1.unwrap() * delta;
    let z_2 = n + 0.5 + direction * alpha_2.unwrap() * delta;
    // dbg!(f_alpha(alpha_1.unwrap()), alpha_1);
    // dbg!(f_alpha(alpha_2.unwrap()), alpha_2);

    let mut m = 1.0;
    let mut h;
    let mut z_l;
    let mut z_r;

    loop {
        h = (z_2 - z_1) / m; // ???
        z_l = ((2. * h).powi(2).inv() + s / Complex::i() / 2.0 / PI).sqrt() + (2. * h).inv();
        z_r = ((2. * h).powi(2).inv() + s / Complex::i() / 2.0 / PI).sqrt() - (2. * h).inv();
        let err1 = g(z_l, s) * (-2.0 * (z_l - z_1) / h * Complex::i() * PI).exp();
        let err2 = g(z_r, s) * (2.0 * (z_r - z_1) / h * Complex::i() * PI).exp();
        if err1.norm() <= eps && err2.norm() <= eps { break }
        m *= 2.0;
    }

    let n_l = (z_l.re - z_l.im).ceil();
    let n_r = (z_r.re - z_r.im).floor();

    let mut s0 = Complex::new(0.0, 0.0);
    for i in 1..=n as i64 {
        s0 += Complex::new(i as Float, 0.0).powc(-s);
    }

    let mut s1 = Complex::new(0.0, 0.0);
    for i in 0..=m as i64 {
        s1 += f(z_1 + i as Float * h, s);
    }

    let mut s2 = Complex::new(0.0, 0.0);
    for i in n_l as i64..=n as i64 {
        s2 += Complex::new(i as Float, 0.0).powc(-s) * H((i as Float - z_1) / h);
    }

    let mut s3 = Complex::new(0.0, 0.0);
    for i in n as i64 + 1 ..= n_r as i64 {
        s3 += Complex::new(i as Float, 0.0).powc(-s) * H((z_1 - i as Float) / h);
    }
    // dbg!(s0, s1, s2, s3);

    s0 + s1 * h - s2 + s3
    // dbg!(n, m, n_l, n_r, z_1, z_2, (z_2 - z_1).norm());
    // dbg!(h);
}

pub fn zeta(s: Complex, eps: Float) -> Complex {
    let chi = Complex::new(2.0 * PI, 0.0).powc(s) / (PI * s / 2.0).cos() / gamma(s, eps) / 2.0;
    // dbg!(chi, gamma(s, eps));
    // dbg!(I0(ctx, s, eps / 2.0), I0(ctx, 1.0 - s.conj(), eps / 2.0));
    I0(s, eps / 2.0) + chi * I0(1.0 - s.conj(), eps / 2.0).conj()
}


#[cfg(test)]
mod tests {
    use super::{zeta, Complex};

    #[test]
    fn test_zeta() {
        let eps = 1e-10;
        let s = Complex::new(1.11, 100.0);
        let result = zeta(s, eps);
        assert!((result - Complex::new(1.528298620331271, -0.07024178360)).norm() < eps);
    }
}
