// Riemann-Siegel formula in the critical line.

use num::Complex;

use crate::context::Context;
use crate::{
    power_series::PowerSeries,
    traits::{ComplexFunctions, MyReal},
    unchecked_cast::UncheckedCast,
};

/// Gabcke's function
// fn gen_gabcke_power_series<T: MyFloat, const N: usize>(z: T) -> PowerSeries<T, N> {
//     let mut numer = PowerSeries::<T, N>::new(
//         z * z * T::FRAC_PI_2() + 3.0f64.unchecked_cast::<T>() * T::FRAC_PI_8(),
//     );
//     numer.data[1] = T::PI() * z;
//     numer.data[2] = T::FRAC_PI_2();
//     numer.cos_();

//     let mut denom = PowerSeries::<T, N>::new(z * T::PI());
//     denom.data[1] = T::PI();
//     denom.cos_();

//     numer /= &denom;
//     numer
// }

/// Juan de Reyna's function
/// F(z) = F_re(z) + i F_im(z), where F_re(z) is the function used by Gabcke (then divide by 2):
/// F_re(z) = cos(PI (z^2 / 2 + 3/8)) / cos(PI z) / 2
/// F_im(z) = (sin(PI (z^2 / 2 + 3/8)) + \sqrt(2) cos(PI z / 2)) / cos(PI z) / 2
/// I guess computing real part and image part separately can be faster, as all operations are real.
fn gen_rs_power_series<T: MyReal>(z: T, N: usize) -> PowerSeries<Complex<T>> {
    // (z^2 / 2 + 3/8) PI
    let zpoly = PowerSeries::<T>::from_vec(
        N,
        vec![
            z * z * T::FRAC_PI_2() + 3.0f64.unchecked_cast::<T>() * T::FRAC_PI_8(),
            T::PI() * z,
            T::FRAC_PI_2(),
        ],
    );

    let mut denom = PowerSeries::<T>::from_vec(N, vec![z * T::PI(), T::PI()]);
    denom.cos_();
    denom *= 2.0f64.unchecked_cast::<T>();

    let mut real_series = zpoly.clone();
    real_series.cos_();
    real_series /= &denom;

    let mut imag_series = zpoly;
    imag_series.sin_();
    let mut imag_series2 = PowerSeries::<T>::from_vec(N, vec![z * T::FRAC_PI_2(), T::FRAC_PI_2()]);
    imag_series2.cos_();
    imag_series2 *= T::SQRT_2();
    imag_series -= &imag_series2;

    imag_series /= &denom;

    let data: Vec<_> = real_series
        .data
        .iter()
        .zip(imag_series.data.iter())
        .map(|(&re, &im)| Complex::new(re, im))
        .collect();
    PowerSeries::<Complex<T>>::from_vec(N, data)
}

pub struct RiemannSiegelZeta<'a, T> {
    ctx: &'a Context<T>,
    sigma: T,
    K: usize,
    coeffs1: Vec<Vec<T>>,
    coeffs2: Vec<Vec<T>>,
    planners: [RiemannSiegelPlanner<T>; 2],
}

pub struct RiemannSiegelPlanner<T> {
    sigma: f64,
    sum_trunc_dirichlet: Vec<Complex<T>>,
}

type Plan<T> = (usize, Option<Complex<T>>);

impl<T: MyReal> RiemannSiegelPlanner<T> {
    pub fn new(sigma: f64) -> Self { Self { sigma, sum_trunc_dirichlet: vec![] } }

    pub fn plan(&self, t: f64, eps: f64) -> Option<Plan<T>> {
        let a = (t.unchecked_cast::<f64>() / f64::PI() / 2.0).sqrt();
        let n = a.floor();
        let sigma: f64 = self.sigma.unchecked_cast();

        let mut k;
        let c1;
        if sigma >= 0.0 {
            k = 0usize;
            c1 = 2.0f64.powf(3.0 * sigma / 2.0) / 7.0;
        } else {
            k = (-sigma).ceil() as usize;
            c1 = 0.9f64.powi(k as i32) / 2.0;
        }

        let eps = eps * a.powf(sigma);
        let b = 10.0 / 11.0 * a;

        while k <= 50 {
            let err =
                c1 * rgsl::gamma_beta::gamma::gamma((k as f64 + 1.0) / 2.0) / b.powi(k as i32 + 1);
            if err < eps {
                return Some((k, None));
            }
            k += 1;
        }

        None
    }
}

impl<'a, T: MyReal> RiemannSiegelZeta<'a, T> {
    fn gen_coeffs(ctx: &'a Context<T>, sigma: T, K: usize) -> Vec<Vec<T>> {
        let mut coeffs: Vec<Vec<T>> = vec![];
        let w = T::one() - 2.0f64.unchecked_cast::<T>() * sigma;
        for k in 0..=K {
            let J = 3 * k / 2;
            let mut d = vec![T::zero(); J + 1];
            for j in 0..=J {
                let r = 3 * k - 2 * j;
                if k == 0 && j == 0 {
                    d[j] = T::one();
                } else if r == 0 {
                    let mut s = T::zero();
                    for i in 0..J {
                        if (J - i) % 2 == 0 {
                            s -= d[i] * ctx.factorial(3 * k - 2 * i) / ctx.factorial(J - i);
                        } else {
                            s += d[i] * ctx.factorial(3 * k - 2 * i) / ctx.factorial(J - i);
                        }
                    }

                    d[j] = s;
                } else {
                    d[j] = T::zero();
                    if r >= 3 {
                        // 2 * j <= 3 * (k - 1)  <=>  r >= 2
                        d[j] += coeffs[k - 1][j] * 0.5f64.unchecked_cast::<T>();
                    }
                    if 1 <= j && r >= 1 {
                        // 2 * (j - 1) <= 3 * (k - 1)  <=>  r >= 1
                        d[j] += w * coeffs[k - 1][j - 1];
                    }
                    if 2 <= j {
                        // 3(k-1) - 2(j-2) = r + 1 > 0
                        let c = ((2 * r * (r + 1)) as i32).unchecked_cast::<T>();
                        d[j] -= c * coeffs[k - 1][j - 2];
                    }
                    d[j] /= ((2 * r) as i32).unchecked_cast::<T>();
                }
            }
            coeffs.push(d);
        }
        for k in 0..=K {
            for j in 0..=3 * k / 2 {
                coeffs[k][j] *= ctx.factorial(3 * k - 2 * j) / ctx.pow_pi(2 * k - j) / ctx.pow2(j);
            }
        }

        coeffs
    }

    pub fn new(ctx: &'a Context<T>, sigma: T, K: usize) -> Self {
        let coeffs1 = RiemannSiegelZeta::gen_coeffs(ctx, sigma, K);
        let coeffs2 = RiemannSiegelZeta::gen_coeffs(ctx, T::one() - sigma, K);
        Self {
            ctx,
            sigma,
            K,
            coeffs1,
            coeffs2,
            planners: [
                RiemannSiegelPlanner::new(sigma.unchecked_cast::<f64>()),
                RiemannSiegelPlanner::new(1.0 - sigma.unchecked_cast::<f64>()),
            ],
        }
    }

    pub fn correct(
        &self,
        a: T,
        K: usize,
        ps: &PowerSeries<Complex<T>>,
        reflect: bool,
    ) -> Complex<T> {
        let mut ret = Complex::<T>::zero();

        let coeffs = if reflect { &self.coeffs2 } else { &self.coeffs1 };
        for k in 0..=K {
            let mut s = Complex::<T>::zero();
            for j in 0..=3 * k / 2 {
                s += ps.data[3 * k - 2 * j].mul_pow_i(4 - j % 4) * coeffs[k][j];
            }
            ret += s / a.powi(k as i32);
        }
        ret
    }

    pub fn solve(
        &self,
        t: T,
        ps: &PowerSeries<Complex<T>>,
        plan: Plan<T>,
        reflect: bool,
    ) -> Complex<T> {
        let ctx = self.ctx;
        let two = ctx.two();
        let a = (t / T::PI() / two).sqrt();
        let n = a.floor();
        let sigma = if reflect { T::one() - self.sigma } else { self.sigma };

        let (k, plan_sum_trunc_dirichlet) = plan;
        let mut sum_trunc_dirichlet;
        match plan_sum_trunc_dirichlet {
            Some(x) => {
                sum_trunc_dirichlet = x;
            }
            None => {
                sum_trunc_dirichlet = Complex::<T>::zero();
                let s = Complex::new(sigma, t);
                let n = n.unchecked_cast::<i32>();
                for i in 1..=n {
                    sum_trunc_dirichlet += (-s * i.unchecked_cast::<T>().ln()).exp();
                }
            }
        }
        let logU = t / two * (t / two / T::PI()).ln() - t / two - T::FRAC_PI_8();
        let U = Complex::new(T::zero(), -logU).exp();

        sum_trunc_dirichlet
            - self.correct(a, k, ps, reflect)
                * a.powf(-sigma)
                * U.mul_pow_i(2 * (n.unchecked_cast::<i32>() as usize))
    }

    pub fn zeta(&self, z: Complex<T>, eps: f64) -> Option<Complex<T>> {
        let log_chi = (z - 0.5f64.unchecked_cast::<T>()) * T::PI().ln()
            + self.ctx.loggamma((T::one() - z) * 0.5f64.unchecked_cast::<T>(), eps / 3.0)
            - self.ctx.loggamma(z * 0.5f64.unchecked_cast::<T>(), eps / 3.0);
        assert!(!log_chi.re.is_nan(), "{:?} {}", z, eps);
        let chi = log_chi.exp();

        let t = z.im;
        let plan0 = self.planners[0].plan(t.unchecked_cast(), eps / 3.0)?;
        let plan1 = self.planners[1].plan(t.unchecked_cast(), eps / 3.0)?;

        let K = std::cmp::max(plan0.0, plan1.0);
        let two = self.ctx.two();
        let a = (t / T::PI() / two).sqrt();
        let n = a.floor();
        let p = T::one() - two * (a - n);
        let ps = gen_rs_power_series(p, 3 * K + 1);

        println!("power series: K = {}", K);
        // println!("plan: {:?} {:?}", plan0, plan1);
        let R0 = self.solve(t, &ps, plan0, false);
        let R1 = self.solve(t, &ps, plan1, true);
        let ret = R0 + chi * R1.conj();
        Some(ret)
    }
}

#[cfg(test)]
mod tests {
    use crate::f64x2;
    use crate::*;
    use num::Complex;

    #[test]
    fn rszeta() {
        type T = f64;
        let ctx = Context::<T>::new(100);
        // let mut zeta_galway = ZetaGalway::new(&ctx);
        let zeta_rs = RiemannSiegelZeta::new(&ctx, T::from(1.5), 50);

        let s = Complex::new(T::from(1.5), f64::from(1000000.0));
        println!("s = {}, zeta(s) = {}", s, zeta_rs.zeta(s, 1e-17).unwrap());
        panic!();
    }
}
