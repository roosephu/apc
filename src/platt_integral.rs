use crate::cache_stat::CacheStat;
use crate::traits::ComplexFunctions;
use crate::{power_series::PowerSeries, traits::MyReal};
use log::{debug, info};
use num::Complex;

/// We approximate Φ(σ + it) ≈ poly(t - t0) * C * ϕ(σ + it)
struct ExpansionIntegrator<T> {
    t: T, // center
    radius: T,
    w: T,
    coeff: Complex<T>,
    poly: PowerSeries<Complex<T>>,
}

impl<T: MyReal> ExpansionIntegrator<T> {
    /// expand exp(a z^2 + c z) at z = 0
    fn expand_exp_term(order: usize, a: Complex<T>, c: Complex<T>) -> PowerSeries<Complex<T>> {
        let mut poly = PowerSeries::<Complex<T>>::from_vec(order, vec![Complex::<T>::zero(), c, a]);
        let mut derivatives = vec![Complex::<T>::one(); order];
        for i in 1..order {
            derivatives[i] = derivatives[i - 1] / T::from_u32(i as u32).unwrap();
        }
        poly.compose_(&derivatives);

        poly
    }

    /// expand 1/(1 - bz)  at z = 0: 1 + b z + b^2 z^2 + ...
    fn expand_reciprocal_term(order: usize, b: Complex<T>) -> PowerSeries<Complex<T>> {
        let mut data = vec![Complex::zero(); order];
        let mut coeff = Complex::one();
        for i in 0..order {
            data[i] = coeff;
            coeff *= b;
        }
        PowerSeries::from_vec(order, data)
    }

    /// See Readme
    /// Expand f(w h) = exp(-λ^2 h^2 / 2 + i h λ^2 s_0) / (1 + i h / s_0) with w = i ln x.
    /// another option is to expand f(w h) = exp(-λ^2 h^2 / 2) / (1 + i h / s_0) with w = i (λ^2 s_0 + ln x).
    fn expand_at(s0: Complex<T>, ln_x: T, λ: T, order: usize) -> (T, PowerSeries<Complex<T>>) {
        let w = ln_x;
        let a = λ * λ / 2.0 / w / w;
        let b = -T::one() / s0 / ln_x;
        let c = s0 * (λ * λ / ln_x);
        let mut ps_exp = Self::expand_exp_term(order, Complex::<T>::new(a, T::zero()), c);
        let ps_recip = Self::expand_reciprocal_term(order, b);
        ps_exp *= &ps_recip;

        (w, ps_exp)
    }

    /// we want to evaluate $$ \sum_i poly_i (w h)^i$$ for a real $h$, however, poly_i c^i might grow too fast...
    /// Complex operations are slow.
    /// trick: we prepare poly'_i = poly_i (c/|c|)^i, so it becomes $$\sum_i poly'_i (|c| h)^i$$
    /// so it's numerically more stable.
    fn normalize_(poly: &mut PowerSeries<Complex<T>>, w: Complex<T>) -> (T, Complex<T>) {
        let w_norm = w.norm();
        let w_dir = w / w_norm;
        let mut w_dir_pow = Complex::<T>::one();
        for i in 0..poly.N {
            poly.coeffs[i] *= w_dir_pow;
            w_dir_pow *= w_dir;
        }
        (w_norm, w_dir)
    }

    fn new(t: T, σ: T, x: T, λ: T, eps: f64, order: usize) -> Self {
        let ln_x = x.ln();
        let s = Complex::new(σ, t);
        let (w, mut poly) = Self::expand_at(s, ln_x, λ, order);
        let hat_phi = (λ * λ / 2.0 * s * s + s * ln_x).exp() / s;
        let coeff = hat_phi / Complex::<T>::new(T::zero(), w);

        for i in 0..order {
            let mut signed_falling_factorial = T::one();
            for j in 0..i {
                signed_falling_factorial *= -T::from_usize(i - j).unwrap();
                let delta = poly.coeffs[i] * signed_falling_factorial;
                poly.coeffs[i - j - 1] += delta;
            }
        }

        let poly_eps = T::from_f64(eps).unwrap() / coeff.norm() / T::from_usize(order).unwrap();
        let radius =
            (poly_eps / poly.coeffs[order - 1].norm()).powf(T::one() / (order as f64 - 1.0)) / w;
        debug!("[Expansion] prepare {}, radius = {}", t, radius);

        // we've expanded in a new t0, and should clear cache.
        Self::normalize_(&mut poly, Complex::<T>::new(T::zero(), w));
        Self { t, radius, w, coeff, poly }
    }

    pub fn query(&self, t: T) -> Complex<T> {
        let diff = t - self.t;
        assert!(diff.abs() <= self.radius, "t = {}, t0 = {}", t, self.t);

        let z = self.w * diff;
        let mut poly = Complex::<T>::zero();
        let mut z_pow = T::one();
        for i in 0..self.poly.n {
            poly += self.poly.coeffs[i] * z_pow;
            z_pow *= z;
        }
        let (sin_z, cos_z) = z.sin_cos();
        let exp_i_z = Complex::new(cos_z, sin_z);

        self.coeff * (poly * exp_i_z - self.poly.coeffs[0])
    }
}

fn get_expansions<T: MyReal>(
    x: T,
    σ: T,
    λ: T,
    max_order: usize,
    limit: T,
    eps: f64,
) -> Vec<ExpansionIntegrator<T>> {
    let mut expansions = vec![];

    let mut t = T::zero();
    while t < limit {
        let expansion = ExpansionIntegrator::new(t, σ, x, λ, eps, max_order);
        t += expansion.radius / 1.2;
        expansions.push(expansion);
    }
    expansions
}

/// to compute $$\int x^s exp(λ^2 s^2 / 2) / s dh$$ for $s = \σ + ih$
pub struct PlattIntegrator<T> {
    expansions: Vec<ExpansionIntegrator<T>>,
    to_inf: Vec<Complex<T>>,
}

impl<T: MyReal> PlattIntegrator<T> {
    pub fn new(x: T, σ: T, λ: T, max_order: usize, eps: f64) -> Self {
        let limit = (x.ln() + x.ln().ln()).sqrt() / λ * 2.0;
        let expansions = get_expansions(x, σ, λ, max_order, limit, eps);

        let mut to_inf = vec![Complex::<T>::zero(); expansions.len()];
        let mut t = limit;
        for (i, expansion) in expansions.iter().enumerate().rev() {
            debug!("add to inf: t = {}, t0 = {}, radius = {}", t, expansion.t, expansion.radius);
            if i + 1 != expansions.len() {
                to_inf[i] = to_inf[i + 1] + expansion.query(t);
            }
            t = expansion.t;
        }
        info!(
            "[PlattIntegrator] up to {}, {} segments, integral = {}",
            limit,
            expansions.len(),
            to_inf[0],
        );

        Self { expansions, to_inf }
    }

    pub fn query(&mut self, t: T) -> Complex<T> {
        let index = self.expansions.partition_point(|key| key.t <= t);
        assert!(index > 0);
        let index = index - 1;
        self.to_inf[index] - self.expansions[index].query(t)
    }
}

pub struct HybridPrecIntegrator<T> {
    expansions: Vec<ExpansionIntegrator<f64>>,
    to_inf: Vec<Complex<T>>,
    pub max_err: f64,
}

impl<T: MyReal> HybridPrecIntegrator<T> {
    pub fn new(x: T, σ: T, λ: T, max_order: usize, eps: f64) -> Self {
        let high_prec = PlattIntegrator::new(x, σ, λ, max_order, eps);
        let mut expansions = vec![];

        let (x_, σ_, λ_) = (x.to_f64().unwrap(), σ.to_f64().unwrap(), λ.to_f64().unwrap());
        for high_prec_expansion in high_prec.expansions.iter() {
            let t = high_prec_expansion.t.to_f64().unwrap();
            let low_prec_expansion = ExpansionIntegrator::new(t, σ_, x_, λ_, eps, max_order);
            expansions.push(low_prec_expansion)
        }

        Self { expansions, to_inf: high_prec.to_inf, max_err: 0.0 }
    }

    pub fn query(&mut self, t: T) -> Complex<T> {
        let t = t.to_f64().unwrap();
        let index = self.expansions.partition_point(|key| key.t <= t);
        assert!(index > 0);
        let index = index - 1;
        let a = self.to_inf[index];
        let b = self.expansions[index].query(t);
        self.max_err += b.norm();
        Complex::new(a.re - b.re, a.im - b.im)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use F64x2::f64x2;
    use F64x2::test_utils::*;

    #[test]
    fn test_platt_integral() {
        // env_logger::init();

        type T = f64x2;
        let σ = T::from_f64(0.5).unwrap();
        let x = T::from_f64(1e6).unwrap();
        let λ = T::from_f64(0.003).unwrap();
        let mut integrator = PlattIntegrator::new(x, σ, λ, 20, 1e-20);

        let t1 = T::zero() + 14.0;
        let t2 = T::from_f64(8e3).unwrap();
        println!("{}", integrator.query(T::zero()));

        // panic!();
    }
}
