use crate::cache_stat::CacheStat;
use crate::traits::ComplexFunctions;
use crate::{power_series::PowerSeries, traits::MyReal};
use log::debug;
use num::Complex;

type Expansion<T> = (T, T, T, Complex<T>, PowerSeries<Complex<T>>);

/// to compute $$\int x^s exp(λ^2 s^2 / 2) / s dh$$ for $s = \σ + ih$
pub struct PlattIntegrator<T> {
    order: usize,
    eps: T,
    σ: T,
    λ_sqr: T,
    ln_x: T,
    expansion: Option<Expansion<T>>,
    cache: Option<(T, Complex<T>)>, // save the last query
    pub stat: CacheStat,
}

impl<T: MyReal> PlattIntegrator<T> {
    pub fn new(x: T, σ: T, λ: T, max_order: usize, eps: f64) -> Self {
        Self {
            ln_x: x.ln(),
            σ,
            eps: T::from_f64(eps).unwrap(),
            order: max_order,
            λ_sqr: λ * λ,
            expansion: None,
            cache: None,
            stat: CacheStat::new(),
        }
    }

    pub fn hat_phi(&self, s: Complex<T>) -> Complex<T> {
        (self.λ_sqr / T::from_f64(2.0).unwrap() * s * s + s * self.ln_x).exp() / s
    }

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
    fn expand_at(&self, s0: Complex<T>, order: usize) -> (T, PowerSeries<Complex<T>>) {
        let w = self.ln_x;
        let a = self.λ_sqr / 2.0 / w / w;
        let b = -T::one() / s0 / self.ln_x;
        let c = s0 * (self.λ_sqr / self.ln_x);
        let mut ps_exp =
            PlattIntegrator::expand_exp_term(order, Complex::<T>::new(a, T::zero()), c);
        let ps_recip = PlattIntegrator::expand_reciprocal_term(order, b);
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

    #[inline(never)]
    fn prepare(&mut self, t: T) {
        self.stat.miss();

        let s0 = Complex::new(self.σ, t);
        let (w, mut poly) = self.expand_at(s0, self.order);
        let mul_coeff = self.hat_phi(s0) / Complex::<T>::new(T::zero(), w);

        for i in 0..self.order {
            let mut signed_falling_factorial = T::one();
            for j in 0..i {
                signed_falling_factorial *= -T::from_usize(i - j).unwrap();
                let delta = poly.coeffs[i] * signed_falling_factorial;
                poly.coeffs[i - j - 1] += delta;
            }
        }

        let poly_eps = self.eps / mul_coeff.norm() / T::from_usize(self.order).unwrap();
        let radius = (poly_eps / poly.coeffs[self.order - 1].norm())
            .powf(T::one() / (self.order as f64 - 1.0))
            / w;
        debug!("[Integral] prepare {}, radius = {}", t, radius);

        // we've expanded in a new t0, and should clear cache.
        self.cache = None;

        Self::normalize_(&mut poly, Complex::<T>::new(T::zero(), w));
        self.expansion = Some((t, radius, w, mul_coeff, poly));
    }

    /// assuming the power series converges well at given point s.
    #[inline(never)]
    fn _query(&mut self, t: T) -> Complex<T> {
        if let Some((k, v)) = self.cache {
            if k == t {
                return v;
            }
        }

        self.stat.hit();
        let (t0, _, w, mul_coeff, ps) = self.expansion.as_ref().unwrap();
        let &w = w;

        let z = w * (t - *t0);
        let mut poly = Complex::<T>::zero();
        let mut z_pow = T::one();
        for i in 0..self.order {
            poly += ps.coeffs[i] * z_pow;
            z_pow *= z;
        }
        let (sin_z, cos_z) = z.sin_cos();
        let exp_i_z = Complex::new(cos_z, sin_z);
        mul_coeff * (poly * exp_i_z - ps.coeffs[0])
    }

    #[inline(never)]
    pub fn query(&mut self, mut t1: T, t2: T) -> Complex<T> {
        // debug!("[integral] query [{}, {}]", t1, t2);
        if let Some((t0, radius, _, _, _)) = self.expansion.as_ref() {
            if (t1 - *t0).abs() < *radius && (t2 - t1).abs() < *radius {
                // debug!("[integral] easy case.");

                let f_t1 = self._query(t1);
                let f_t2 = self._query(t2);
                self.cache = Some((t2, f_t2));
                return f_t2 - f_t1;
            }
        }

        let mut result = Complex::zero();
        loop {
            self.prepare(t1);
            let (_, radius, _, _, _) = self.expansion.as_ref().unwrap();
            let radius = *radius;

            let max_t = t1 + radius;
            if max_t > t2 {
                result += self._query(t2);
                break;
            } else {
                result += self._query(max_t);
                t1 = max_t;
            }
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::*;
    use F64x2::f64x2;

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
        println!("{}", integrator.query(T::zero(), t2));

        // panic!();
    }
}
