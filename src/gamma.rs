use std::f64::consts::PI;

use num::traits::*;
use rustfft::FftPlanner;

type Float = f64;
type Complex = num::Complex<Float>;

// https://epubs.siam.org/doi/pdf/10.1137/0731050
pub fn gamma(z: Complex, eps: Float) -> Complex {
    let z = z - 1.0;
    let mut a = (-z.re).max(2.0).floor() + 0.5;
    loop {
        let err = a.sqrt() / (2.0 * PI).powf(a + 0.5) / (z.re + a);
        if err < eps { break }
        a += 1.0;
    }
    let mut coef = Complex::new(1.0, 0.0);
    let mut k = 1.0;
    let mut fac = 1.0;
    while k < a {
        let c_k = (a - k).powf(k - 0.5) * (a - k).exp() * (-1.0).pow(k as i32 - 1) / (2.0 * PI).sqrt() / fac;
        fac *= k;
        coef += c_k / (z + k);
        k += 1.0;
    }
    // dbg!(coef, (z + a).powc(z + 0.5) / (z + a).exp() * (2.0 * PI).sqrt());
    coef * (z + a).powc(z + 0.5) / (z + a).exp() * (2.0 * PI).sqrt()
}

pub struct Context {
    bernoulli: Vec<Float>,
}

impl Context {
    fn new() -> Self {
        Self { bernoulli: vec![] }
    }

    const EXP_APPROX: [Complex; 18] = [
        Complex { re: 1.0, im: 1.0812144107799266e-23},
        Complex { re: -4.479610786345225e-21, im: 1.0},
        Complex { re: -0.5, im: 3.027415987199093e-19},
        Complex { re: -8.413828845781633e-18, im: -0.16666666666666669},
        Complex { re: 0.04166666666666679, im: 1.1656663448809618e-16},
        Complex { re: -1.0603511152022404e-15, im: 0.008333333333332324},
        Complex { re: -0.0013888888888827402, im: 5.789264486508273e-15},
        Complex { re: -2.491995923872859e-14, im: -0.00019841269843586228},
        Complex { re: 2.4801587374768556e-5, im: 6.704175576866034e-14},
        Complex { re: -1.594987515102099e-13, im: 2.7557317787217356e-6},
        Complex { re: -2.755729303110001e-7, im: 2.3127460502103687e-13},
        Complex { re: -3.2663668749921504e-13, im: -2.5052389834713885e-8},
        Complex { re: 2.087985316554709e-9, im: 2.5867211760028217e-13},
        Complex { re: -2.2167241850689593e-13, im: 1.6041263496425594e-10},
        Complex { re: -1.1352710114429515e-11, im: 8.943908448871146e-14},
        Complex { re: -4.542339711641447e-14, im: -7.962911435347713e-13},
        Complex { re: 5.979573239083729e-14, im: 7.185782517642856e-15},
        Complex { re: -1.970149077208406e-15, im: 1.9701490772084063e-15},
    ];

    #[inline]
    pub fn bernoulli(&self, idx: usize) -> Float {
        self.bernoulli[idx]
    }

    fn init_bernoulli(&mut self) {
        unimplemented!()
    }

    pub fn loggamma(&self, mut z: Complex, eps: Float) -> Complex {
        let mut result = Complex::new(0.0, 0.0);
        while z.re < 20.0 {
            result -= z.ln();
            z += 1.0;
        }

        result += (z - 0.5) * z.ln() - z + (PI * 2.0).ln() / 2.0;
        let z2 = z * z;
        let mut zpow = z;
        for i in 1..20 {
            result += self.bernoulli(i * 2) / ((2 * i) * (2 * i - 1)) as f64 / zpow;
            zpow *= z2;
        }
        result
    }

    fn I0(&self, s: Complex, eps: Float) -> Complex {
        unimplemented!()
    }

    pub fn zeta(&self, s: Complex, eps: Float) -> Complex {
        // zeta Galway
        let log_chi = (s - 0.5) * PI.ln() + self.loggamma((1.0 - s) / 2.0, eps) - self.loggamma(s / 2.0, eps);
        let chi = log_chi.exp();
        self.I0(s, eps) + chi * self.I0(1.0 - s.conj(), eps).conj()
    }

    fn ifft(&self, f: &mut Vec<Complex>) {
        let mut planner = FftPlanner::new();
        let fft = planner.plan_fft_inverse(f.len());
        fft.process(f)
    }

    // compute \sum_{n=1}^N n^{-(s + i n t)} for every 0 <= t < M
    pub fn sum_trunc_dirichlet(&self, N: usize, s: Complex, M: usize, delta: Float) -> Vec<Complex> {
        let M2 = (M + M % 2) / 2;
        let s = s + Complex::new(0.0, M2 as f64 * delta);
        let a: Vec<_> = (1..=N).map(|x| Complex::new(x as Float, 0.0).powc(-s)).collect();
        let g: Vec<_> = (1..=N).map(|x| (x as Float).ln() * -delta).collect();
        let R = (M + 1).next_power_of_two() * 2;
        let div = 2.0 * PI / R as Float;
        let w: Vec<_> = g.iter().map(|x| ((x / div).round() as i64)).collect();
        let d: Vec<_> = g.iter().zip(w.iter()).map(|(x, &y)| x - y as f64 * div).collect();

        let mut ret = vec![Complex::zero(); M + 1];
        let mut f= vec![Complex::zero(); R];
        for (e, c) in Self::EXP_APPROX.iter().enumerate() {
            for x in f.iter_mut() { x.set_zero() }
            for j in 0..N {
                f[w[j].rem_euclid(R as i64) as usize] += a[j] * (d[j] * M2 as f64).pow(e as i32);
            }
            self.ifft(&mut f);
            for x in 0..=M {
                ret[x] += c * f[(x + R - M2) % R] * (x as f64 / M2 as f64 - 1.0).pow(e as i32);
            }
        }

        ret
    }
}

#[cfg(test)]
mod tests {
    use super::{gamma, Complex, Float, Context};
    use num::Zero;

    fn _test_gamma(z: Complex, eps: Float, gt: Complex) {
        let result = gamma(z, eps);
        let diff = (result - gt).norm();
        println!("{:?} {:?}", result, gt);
        assert!(diff <= eps);
    }
    #[test]
    fn test_gamma() {
        _test_gamma(Complex::new(4.0, 10.0), 1e-10, Complex::new(0.0007715342942399662, -0.0010190827990417));
        _test_gamma(Complex::new(4.0, 0.0), 1e-10, Complex::new(6.0, 0.0));
    }

    #[test]
    fn test_sum_trunc_dirichlet() {
        let ctx = Context::new();
        let N = 20;
        let s = Complex::new(1.3, 2.6);
        let M = 9;
        let delta = 1.1;
        let result = ctx.sum_trunc_dirichlet(N, s, M, delta);
        for i in 0..=M {
            let mut sum = Complex::zero();
            let z = s + Complex::new(0.0, i as f64 * delta);
            for j in 1..=N {
                sum += Complex::new(j as f64, 0.0).powc(-z);
            }
            assert!((sum - result[i]).norm() < 2e-15);
        }
    }

    #[test]
    fn bench_sum_trunc_dirichlet() {
        let ctx = Context::new();
        let N = 20;
        let s = Complex::new(1.3, 2.6);
        let M = 100_000;
        let delta = 1.1;
        let result = ctx.sum_trunc_dirichlet(N, s, M, delta);
        for x in result {
            assert!(x.norm() != 0.0);
        }
    }
}
