use crate::traits::ComplexFunctions;
use crate::{
    power_series::PowerSeries,
    traits::MyReal,
    unchecked_cast::{UncheckedCast, UncheckedInto},
};
use log::debug;
use num::Complex;

type Expansion<T> = (T, T, (T, Complex<T>), Complex<T>, PowerSeries<Complex<T>>);

/// to compute $$\int x^s exp(lambda^2 s^2 / 2) / s dh$$ for $s = \sigma + ih$
pub struct PlattIntegrator<T> {
    order: usize,
    eps: T,
    sigma: T,
    lambda_sqr: T,
    ln_x: T,
    expansion: Option<Expansion<T>>,
    cache: Option<(T, Complex<T>)>, // save the last query
}

impl<T: MyReal> PlattIntegrator<T> {
    pub fn new(x: T, sigma: T, lambda: T, max_order: usize, eps: f64) -> Self {
        Self {
            ln_x: x.ln(),
            sigma,
            eps: eps.unchecked_cast::<T>(),
            order: max_order,
            lambda_sqr: lambda * lambda,
            expansion: None,
            cache: None,
        }
    }

    pub fn hat_phi(&self, s: Complex<T>) -> Complex<T> {
        (self.lambda_sqr / 2.0.unchecked_cast::<T>() * s * s + s * self.ln_x).exp() / s
    }

    /// expand exp(a z^2) at z = 0: 1 + (a / 1!) z^2 + (a^2 / 2!) z^4 + ..
    fn expand_exp_term(order: usize, a: Complex<T>) -> PowerSeries<Complex<T>> {
        let mut data = vec![Complex::zero(); order];
        let mut coeff = Complex::one();
        for i in 0..(order + 1) / 2 {
            data[i * 2] = coeff;
            coeff *= a / (i as f64 + 1.0).unchecked_cast::<T>();
        }
        PowerSeries::from_vec(order, data)
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

    /// expand $\hat_\phi(s)$ at s0
    /// more specifically, $\hat_\phi(s_0 + ih) = \hat_\phi(s_0) exp(c h) exp(a (c h)^2) / (1 - b (c h))$
    // for a = ..., b = ..., c = ...
    /// so the integration becomes $\hat_\phi(s_0)/c \int exp(c h) exp(a (c h)^2) / (1 - b (c h)) d (c h)$
    /// which is $\hat_\phi(s_0) / c \int exp(z) exp(a z^2) / (1 - b z) d z$ for z = c h
    /// now we expand of exp(a z^2) and 1/(1 - b z) separately.
    fn expand_at(
        &self,
        s0: Complex<T>,
        order: usize,
    ) -> (Complex<T>, Complex<T>, Complex<T>, PowerSeries<Complex<T>>) {
        let c = (self.lambda_sqr * s0 + self.ln_x).mul_i();
        let a = self.lambda_sqr / -2.0 / c / c;
        let b = -(T::one() / s0 / c).mul_i();
        let mut ps_exp = PlattIntegrator::expand_exp_term(order, a);
        let ps_recip = PlattIntegrator::expand_reciprocal_term(order, b);
        ps_exp *= &ps_recip;

        (a, b, c, ps_exp)
    }

    /// we want to evaluate $$ \sum_i poly_i (c h)^i$$ for a real $h$, however, poly_i c^i might grow too fast...
    /// Complex operations are slow.
    /// trick: we prepare poly'_i = poly_i (c/|c|)^i, so it becomes $$\sum_i poly'_i (|c| h)^i$$
    /// so it's numerically more stable.
    fn normalize_(poly: &mut PowerSeries<Complex<T>>, c: Complex<T>) -> (T, Complex<T>) {
        let c_norm = c.norm();
        let c_dir = c / c_norm;
        let mut c_dir_pow = Complex::<T>::one();
        for i in 0..poly.N {
            poly.data[i] *= c_dir_pow;
            c_dir_pow *= c_dir;
        }
        (c_norm, c_dir)
    }

    #[inline(never)]
    fn prepare(&mut self, t: T) {
        let s0 = Complex::new(self.sigma, t);
        let (a, b, c, mut poly) = self.expand_at(s0, self.order);
        let mul_coeff = self.hat_phi(s0) / c;

        // now the integral becomes $\hat_\phi(s_0) / c \int exp(z) exp(a z^2) / (1 - b z) d z$
        // and we have $exp(a z^2) / (1 - bz) \approx poly(z)$
        // by integrating by parts: \int exp(z) poly(z) d z = exp(z) poly(z) - \int poly'(z) exp(z) d z
        // so the integration is $\hat_\phi(s_0) / c \times exp(z) (poly(z) - poly'(z) + poly''(z)...)$
        // now we're initializing the last term: f(z) = poly(z) - poly'(z) + poly''(z)...
        for i in 0..self.order {
            let mut signed_falling_factorial = T::one();
            for j in 0..i {
                signed_falling_factorial *= -((i - j) as f64).unchecked_cast::<T>();
                let delta = poly.data[i] * signed_falling_factorial;
                poly.data[i - j - 1] += delta;
            }
        }

        let poly_eps = self.eps / mul_coeff.norm() / (self.order as f64).unchecked_cast::<T>();
        let radius = (poly_eps / poly.data[self.order - 1].norm())
            .powf(T::one() / (self.order as f64 - 1.0))
            / c.norm();
        debug!("[Integral] prepare {}, radius = {}", t, radius);

        // we've expanded in a new t0, and should clear cache.
        self.cache = None;

        let c = Self::normalize_(&mut poly, c);
        self.expansion = Some((t, radius, c, mul_coeff, poly));
    }

    /// assuming the power series converges well at given point s.
    fn _query(&self, t: T) -> Complex<T> {
        if let Some((k, v)) = self.cache {
            if k == t {
                return v;
            }
        }
        let (t0, _, c, mul_coeff, ps) = self.expansion.as_ref().unwrap();
        let &(c_norm, c_dir) = c;

        let z = c_norm * (t - *t0);
        let mut poly = Complex::<T>::zero();
        let mut z_pow = T::one();
        for i in 0..self.order {
            poly += ps.data[i] * z_pow;
            z_pow *= z;
        }
        mul_coeff * (poly * (z * c_dir).exp() - ps.data[0])
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
    use crate::f64x2;
    use crate::test_utils::*;

    #[test]
    fn test_platt_integral() {
        env_logger::init();

        type T = f64x2;
        let sigma = 0.5.unchecked_cast::<T>();
        let x = 1e6.unchecked_cast::<T>();
        let lambda = 0.003.unchecked_cast::<T>();
        let mut integrator = PlattIntegrator::new(x, sigma, lambda, 20, 1e-20);

        let t1 = T::zero() + 14.0;
        let t2 = 8e3.unchecked_cast::<T>();
        println!("{}", integrator.query(T::zero(), t2));

        panic!();
    }
}
