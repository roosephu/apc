use num::{Complex, Float, Signed};

use crate::contexts::{ExpPolyApprox, Sinc};
use crate::sum_trunc_dirichlet::sum_trunc_dirichlet;
use crate::traits::MyReal;
use log::debug;

/// a data structure for querying \sum_{t=1}^k n^{-sigma - i t}
pub struct BandwidthInterp<T> {
    k0: usize,
    k1: usize,
    tau: T,
    sigma: T,
    t0: T,
    gap: T,

    alpha: T,
    beta: T,
    delta: T,
    data: Vec<Complex<T>>,
    h: SinhKernel<T>,
}

/// determine c: c / sinh(c) = eps
fn find_c(eps: f64) -> f64 {
    let c = 1.0;
    let c = (2.0 / eps * c).ln();
    let c = (2.0 / eps * c).ln();
    c
}

fn interpolate<T: MyReal>(data: &[Complex<T>], x: T, beta: T, h: &impl Kernel<T>) -> Complex<T> {
    let mut ret = Complex::<T>::zero();
    let delta = T::PI() / beta;
    let w = h.window();
    let window = w / delta.fp();
    let center = x.fp() / delta.fp();
    let u = beta * x;
    let center_T = u * T::FRAC_1_PI();

    if (center - center.round()).abs() < 1e-11 {
        // x is very close to a grid point.
        // Assuming G is 1-Lipschitz, we can simply use the data in the grid.
        ret = data[center as usize];
    } else {
        let sin_u = u.sin();
        assert!(x.fp() >= w);
        let l = (center - window).ceil() as usize;
        let r = (center + window).floor() as usize;

        for n in l..=r {
            ret += data[n] * h.query(x, n)
                / ((center_T - n as f64) * (if n % 2 == 1 { -1.0f64 } else { 1.0f64 }));
        }
        ret = ret * sin_u * T::FRAC_1_PI();
    }
    ret
}

/// Proposition 3 in [Gourdon].
/// We'd like to interpolate a function G(x) satisfying $|G(z) = O(e^{\tau y})|$.
/// We have a kernel function $h(z)$ such that $h(0) = 1$ and $h(z) = o(e^{
/// \gamma y})$. Let $\beta >= \gamma + \tau$. We have evaluations of $G(x)$ on
/// a grid $\{x = \delta n\}_{0 \leq n < N}$ where $\delta = \pi / \beta$.
///
/// $w$ is the window size of $h$: $|h(w)|$ is very small so that we can
/// safely ignore it.
///
/// Here we use the Sinc kernel in [Gourdon] for its efficency (Sec 2.4.3 in
/// [Gourdon]). That is, the kernel is
/// $h(t) = (sin(\gamma t / m) / (\gamma t/m))^m$. We choose $M$ to be a power
/// of 2.
struct BandwithInterploation<T, K: Kernel<T>> {
    β: T,
    δ: T,
    γ: T,
    m: usize,
    data: Vec<Complex<T>>,
    sin_cos: Vec<(T, T)>,
    kernel: K,
}

trait Kernel<T> {
    fn init(&mut self, δ: T, γ: T, eps: f64);
    fn window(&self) -> f64;
    fn query(&self, x: T, n: usize) -> T;
}

// struct SincKernel<T>();

struct SinhKernel<T> {
    δ: T,
    γ: T,
    c: T,
    c_over_c_sinh: T,
}

impl<T: MyReal> SinhKernel<T> {
    fn new() -> Self {
        Self { δ: T::zero(), γ: T::zero(), c: T::zero(), c_over_c_sinh: T::zero() }
    }
}

impl<T: MyReal> Kernel<T> for SinhKernel<T> {
    #[inline]
    fn init(&mut self, δ: T, γ: T, eps: f64) {
        self.δ = δ;
        self.γ = γ;
        self.c = T::mp(find_c(eps / 10.0));
        self.c_over_c_sinh = self.c / self.c.sinh();
    }

    #[inline]
    fn window(&self) -> f64 { self.c.fp() / (self.γ.fp() / 2.0) }

    fn query(&self, x: T, n: usize) -> T {
        let t = x - self.δ * n as f64;
        let gap = self.γ / 2.0;
        let w = (self.c * self.c - gap * gap * t * t).sqrt();
        if w.is_zero() {
            self.c_over_c_sinh
        } else {
            w.sinh() / w * self.c_over_c_sinh
        }
    }
}

/// TODO: refactor: this should only include bandwithinterp, no how parameters
/// should  be selected.
impl<T: MyReal + Sinc + ExpPolyApprox + Signed> BandwidthInterp<T> {
    /// tau = ln(k1/k0)
    /// precompute = O((c + k1 eps)(tau/gap + 2))
    /// query = O(k0 + c (tau/gap + 1))
    #[inline(never)]
    pub fn new(k: usize, min_t: T, max_t: T, sigma: T, eps: f64) -> Self {
        // let k0_int = std::cmp::min(4, k);
        let k0_int = 1;
        let k0 = T::mp(k0_int as f64);
        let k1 = T::mp(k as f64);
        let tau = (k1 / k0).ln();
        let gap = tau * 2.0;
        let beta = tau + gap * 2.0;
        let alpha = (k0 * k1).ln() / 2.0;
        let max_c = T::mp(80.0);
        let t0 = min_t - max_c / gap;
        let t1 = max_t + max_c / gap;
        let delta = T::PI() / beta;
        let m = ((t1 - t0) / delta).to_usize().unwrap();
        let data = sum_trunc_dirichlet(Complex::new(sigma, t0), k0_int, k, m, delta);
        let data = data
            .iter()
            .enumerate()
            .map(|(i, &x)| Complex::from_polar(T::one(), (t0 + delta * i as f64) * alpha) * x)
            .collect();
        // debug!("precompute {} terms", m);
        let mut h = SinhKernel::new();
        h.init(delta, gap * 2.0, eps);

        Self { k0: k0_int, k1: k, tau, sigma, alpha, beta, data, gap, t0, delta, h }
    }

    #[inline]
    fn h(&self, c: T, t: T) -> T {
        let w = (c * c - self.gap * self.gap * t * t).sqrt();
        if w.is_zero() {
            T::one()
        } else {
            w.sinh() / w
        }
    }

    /// we have all precomputed data of
    /// G1(dt) = G(dt + t0) = e^{i alpha (t0 + dt)} \sum_{n=k0}^{k1} n^{-sigma - i(t0 + dt)}
    /// for dt multipliers of delta = pi/beta
    pub fn query(&self, t: T) -> Complex<T> {
        let dt = t - self.t0;

        let mut ret = interpolate(&self.data, dt, self.beta, &self.h);
        ret *= Complex::from_polar(T::one(), -self.alpha * t);

        let s = Complex::new(self.sigma, t);
        for n in 1..self.k0 {
            ret += (-s * T::mp(n as f64).ln()).exp();
        }

        ret
    }
}

#[cfg(test)]
mod tests {
    use super::BandwidthInterp;
    use crate::traits::MyReal;
    use num::Complex;
    use F64x2::f64x2;
    use F64x2::test_utils::*;

    #[test]
    fn test_bandwidth_limited_interp() {
        type T = f64x2;

        let k = 100;
        let sigma = T::mp(0.5);
        let eps = 1e-26;
        let min_t = T::mp(2.0 * f64::PI() * (k as f64).powi(2));
        let max_t = T::mp(2.0 * f64::PI() * (k as f64 + 1.0).powi(2));
        let ds = BandwidthInterp::<T>::new(k, min_t, max_t, sigma, eps);
        let t = T::PI() * 2.0 * 10100.0;
        // let t = ds.t0;
        let mut gt = Complex::<T>::zero();
        let s = Complex::new(sigma, t);
        for i in 1..=k {
            gt += (-s * T::mp(i as f64).ln()).exp();
        }
        println!("ds params: alpha = {}, beta = {}, tau = {}", ds.alpha, ds.beta, ds.tau);

        let output = ds.query(t);
        // let output = ds.data[0];
        // gt *= Complex::new(T::zero(), ds.alpha * t).exp();
        println!("gt = {}, output = {}, diff = {:.e}", gt, output, gt - output);
        assert_complex_close(output, gt, eps);
    }

    #[test]
    fn test_bandwidth_limited_interp_f64() {
        type T = f64;

        let k = 100;
        let sigma = T::mp(0.5);
        let eps = 1e-10;
        let min_t = T::mp(2.0 * f64::PI() * (k as f64).powi(2));
        let max_t = T::mp(2.0 * f64::PI() * (k as f64 + 1.0).powi(2));
        let ds = BandwidthInterp::<T>::new(k, min_t, max_t, sigma, eps);
        let t = T::PI() * 2.0 * 10100.0;
        let mut gt = Complex::<T>::zero();
        let s = Complex::new(sigma, t);
        for i in 1..=k {
            gt += (-s * T::mp(i as f64).ln()).exp();
        }
        println!("ds params: alpha = {}, beta = {}, tau = {}", ds.alpha, ds.beta, ds.tau);

        let output = ds.query(t);
        println!("gt = {}, output = {}, diff = {:.e}", gt, output, gt - output);
        // (a - b).norm().approx() / b.norm().approx();
        let diff = (output - gt).norm() / gt.norm();
        assert!(diff <= eps);
        // assert_complex_close(output, gt, eps);
    }
}
