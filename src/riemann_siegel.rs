// Riemann-Siegel formula in the critical line.

use num::Complex;

use crate::constants::RS_GABCKE_GAMMA;
use crate::{context::Context, traits::GabckeExpansion};
use crate::{
    power_series::PowerSeries,
    traits::{ComplexFunctions, MyReal},
};

fn riemann_siegel_chi<T: MyReal>(ctx: &Context<T>, z: Complex<T>, eps: f64) -> Complex<T> {
    let half = T::from_f64(0.5).unwrap();
    let log_chi = (z - half) * T::PI().ln() + ctx.loggamma((T::one() - z) * half, eps)
        - ctx.loggamma(z * half, eps);
    assert!(!log_chi.re.is_nan(), "{:?} {}", z, eps);
    let chi = log_chi.exp();
    chi
}

/// Juan de Reyna's function
/// F(z) = F_re(z) + i F_im(z), where F_re(z) is the function used by Gabcke (then divided by 2):
/// F_re(z) = cos(PI (z^2 / 2 + 3/8)) / cos(PI z) / 2
/// F_im(z) = (sin(PI (z^2 / 2 + 3/8)) + \sqrt(2) cos(PI z / 2)) / cos(PI z) / 2
/// I guess computing real part and image part separately can be faster, as all operations are real.
fn gen_rs_power_series<T: MyReal>(z: T, N: usize) -> PowerSeries<Complex<T>> {
    // (z^2 / 2 + 3/8) PI
    let zpoly = PowerSeries::<T>::from_vec(
        N,
        vec![z * z * T::FRAC_PI_2() + T::FRAC_PI_8() * 3.0, T::PI() * z, T::FRAC_PI_2()],
    );

    let mut denom = PowerSeries::<T>::from_vec(N, vec![z * T::PI(), T::PI()]);
    denom.cos_();
    denom *= T::from_f64(2.0).unwrap();

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
        .coeffs
        .iter()
        .zip(imag_series.coeffs.iter())
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
    K: usize,
    sum_trunc_dirichlet: Vec<Complex<T>>,

    /// \approx coeffs[k][0]
    max_coeff: Vec<f64>,

    /// gamma((k + 1) / 2)
    gamma_half: Vec<f64>,
}

type Plan<T> = (usize, Option<Complex<T>>);

impl<T: MyReal> RiemannSiegelPlanner<T> {
    pub fn new(sigma: f64, K: usize) -> Self {
        let gamma_half =
            (0..=K).map(|k| rgsl::gamma_beta::gamma::gamma((k as f64 + 1.0) / 2.0)).collect();
        let max_coeff = (0..=K)
            .map(|k| {
                if k == 0 {
                    1.0
                } else {
                    let k = k as f64;
                    (2.0 * k * k.ln() - 2.0 * k * (2.0 * f64::PI() * f64::E() / 3.0).ln()).exp()
                }
            })
            .collect();
        Self { K, max_coeff, sigma, gamma_half, sum_trunc_dirichlet: vec![] }
    }

    pub fn plan(&self, t: f64, eps: f64) -> Option<Plan<T>> {
        let a = (t.to_f64().unwrap() / f64::PI() / 2.0).sqrt();
        let n = a.floor();
        let sigma = self.sigma.to_f64().unwrap();

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
        let mut bpow = b.powi(k as i32 + 1);

        let coeff_bound = eps / T::epsilon().to_f64().unwrap() / 10.0f64;
        while k <= self.K {
            // bpow = b.powi(k + 1)
            let err = c1 * self.gamma_half[k] / bpow;

            if self.max_coeff[k] > coeff_bound {
                println!(
                    "k = {}, max_coef = {}, coeff bound = {}, err = {:.e}, eps = {:.e}",
                    k, self.max_coeff[k], coeff_bound, err, eps
                );
                return None;
            }

            if err < eps {
                return Some((k, None));
            }
            k += 1;
            bpow *= b;
        }

        None
    }
}

impl<'a, T: MyReal> RiemannSiegelZeta<'a, T> {
    fn gen_coeffs(ctx: &'a Context<T>, sigma: T, K: usize) -> Vec<Vec<T>> {
        let mut coeffs: Vec<Vec<T>> = vec![];
        let w = T::one() - T::from_f64(2.0).unwrap() * sigma;
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
                        d[j] += coeffs[k - 1][j] * T::from_f64(0.5).unwrap();
                    }
                    if 1 <= j && r >= 1 {
                        // 2 * (j - 1) <= 3 * (k - 1)  <=>  r >= 1
                        d[j] += w * coeffs[k - 1][j - 1];
                    }
                    if 2 <= j {
                        // 3(k-1) - 2(j-2) = r + 1 > 0
                        let c = T::from_usize(2 * r * (r + 1)).unwrap();
                        d[j] -= c * coeffs[k - 1][j - 2];
                    }
                    d[j] /= T::from_usize(2 * r).unwrap();
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
                RiemannSiegelPlanner::new(sigma.to_f64().unwrap(), K),
                RiemannSiegelPlanner::new(1.0 - sigma.to_f64().unwrap(), K),
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
                s += ps.coeffs[3 * k - 2 * j].mul_pow_i(4 - j % 4) * coeffs[k][j];
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
        let a = (t / T::PI() / 2.0).sqrt();
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
                let n = n.to_i32().unwrap();
                for i in 1..=n {
                    sum_trunc_dirichlet += (-s * T::from_i32(i).unwrap().ln()).exp();
                }
            }
        }
        let logU = t / 2.0 * (t / 2.0 / T::PI()).ln() - t / 2.0 - T::FRAC_PI_8();
        let U = Complex::new(T::zero(), -logU).exp();

        sum_trunc_dirichlet
            - self.correct(a, k, ps, reflect)
                * a.powf(-sigma)
                * U.mul_pow_i(2 * (n.to_i32().unwrap() as usize))
    }

    pub fn zeta(&self, z: Complex<T>, eps: f64) -> Option<Complex<T>> {
        let chi = riemann_siegel_chi(self.ctx, z, eps / 3.0);

        let t = z.im;
        let plan0 = self.planners[0].plan(t.to_f64().unwrap(), eps / 3.0)?;
        let plan1 = self.planners[1].plan(t.to_f64().unwrap(), eps / 3.0)?;

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
        println!("R0 = {}, R1 = {}", R0, R1);
        let ret = R0 + chi * R1.conj();
        Some(ret)
    }
}

pub struct RiemannSiegelTheta<'a, T> {
    ctx: &'a Context<T>,
    K: usize,
    coeffs: Vec<T>,
}

impl<'a, T: MyReal> RiemannSiegelTheta<'a, T> {
    /// see https://arxiv.org/pdf/1609.03682.pdf for the "wrong" formula
    pub fn new(ctx: &'a Context<T>, K: usize) -> Self {
        // (1 - T(2)^(1 - 2 * j)) * abs(T(_bernoulli[j + 1])) / (4 * j * (2 * j - 1))
        let mut coeffs = vec![T::zero(); K + 1];
        for j in 1..=K {
            coeffs[j] = (T::one() - T::from_f64(2.0).unwrap().powi(1 - 2 * j as i32))
                * ctx.bernoulli(2 * j).abs()
                / (4 * j * (2 * j - 1)) as f64;
        }
        Self { ctx, K, coeffs }
    }

    // See [Sec 3.11, Pugh].
    pub fn theta(&self, t: T, eps: f64) -> T {
        // as it's typically used with RiemannSiegelZ, we hope it's not too small.
        assert!(t.to_f64().unwrap() >= 200.0 && eps > 1e-33);
        const K: usize = 7;

        // needs high precision base computation here.
        let mut ret = t / 2.0 * (t / 2.0 / T::PI() / T::E()).ln() - T::FRAC_PI_8();
        let mut tpow = t;
        let tsqr = t * t;
        for i in 1..=K {
            ret += self.coeffs[i] / tpow;
            tpow *= tsqr;
        }
        ret
    }
}

/// Gabcke's
pub struct RiemannSiegelZ<'a, T> {
    ctx: &'a Context<T>,
    theta: RiemannSiegelTheta<'a, T>,
    K: usize,
    coeffs: Vec<Vec<T>>,
}

impl<'a, T: MyReal + GabckeExpansion> RiemannSiegelZ<'a, T> {
    pub fn new(ctx: &'a Context<T>, K: usize) -> Self {
        let coeffs = RiemannSiegelZ::gen_coeffs(ctx, K);
        Self { ctx, K, coeffs, theta: RiemannSiegelTheta::new(ctx, K) }
    }
}

impl<T: MyReal + GabckeExpansion> RiemannSiegelZ<'_, T> {
    fn gen_coeffs(ctx: &Context<T>, K: usize) -> Vec<Vec<T>> {
        let mut coeffs: Vec<Vec<T>> = vec![];
        let mut lambda = vec![T::zero(); K + 1];
        lambda[0] = T::one();
        for n in 0..=K / 4 {
            for k in 0..=n {
                let d = ctx.pow2(4 * k + 1) * ctx.euler(k + 1).abs() * lambda[n - k];
                lambda[n + 1] += d;
            }
            lambda[n + 1] = lambda[n + 1] / (n + 1) as f64;
        }

        for k in 0..=K {
            let J = 3 * k / 4;
            let mut d = vec![T::zero(); J + 1];
            for j in 0..=J {
                let r = 3 * k - 4 * j;
                if r == 0 {
                    d[j] = lambda[k / 4];
                } else {
                    if r > 2 {
                        d[j] += coeffs[k - 1][j] * ((r - 2) * (r - 1)) as f64;
                    }
                    if j >= 1 {
                        d[j] += coeffs[k - 1][j - 1];
                    }
                }
            }
            coeffs.push(d);
        }
        for k in 0..=K {
            for j in 0..=3 * k / 4 {
                coeffs[k][j] /= ctx.pow_pi(2 * k - 2 * j) * ctx.pow2(2 * k);
            }
        }

        coeffs
    }

    fn plan(&self, t: f64, eps: f64) -> Option<(usize, Option<T>)> {
        if t <= 200.0 {
            return None;
        }

        const COEFFS: [f64; 10] =
            [0.127, 0.053, 0.011, 0.031, 0.017, 0.061, 0.661, 9.2, 130.0, 1837.0];
        let mut pow = t.powf(-0.75);
        let sqrt_t = t.sqrt();
        let coeff_bound = eps / T::epsilon().to_f64().unwrap() / 10.0f64;

        for (k, c) in COEFFS.iter().enumerate() {
            if self.coeffs[k][0].to_f64().unwrap() > coeff_bound {
                return None;
            }
            if pow * c <= eps {
                return Some((k, None));
            }
            pow /= sqrt_t;
        }

        None
    }

    fn solve(&self, t: T, plan: (usize, Option<T>), eps: f64) -> T {
        let two = self.ctx.two();
        let a = (t / T::PI() / 2.0).sqrt();
        let n = a.floor();
        let (K, plan_sum_trunc_dirichlet) = plan;
        let z = Complex::new(T::from_f64(0.25).unwrap(), t * 0.5);
        // let theta = self.ctx.loggamma(z, eps).im - t * 0.5 * T::PI().ln();
        let theta = self.theta.theta(t, eps);

        let mut sum_trunc_dirichlet;
        match plan_sum_trunc_dirichlet {
            Some(x) => {
                sum_trunc_dirichlet = x;
            }
            None => {
                sum_trunc_dirichlet = T::zero();
                let n = n.to_i32().unwrap();
                for i in 1..=n {
                    let i = T::from_i32(i).unwrap();
                    sum_trunc_dirichlet += (theta - t * i.ln()).cos() / i.sqrt();
                }
            }
        }
        let sum_trunc_dirichlet = sum_trunc_dirichlet * 2.0;

        let correction = T::expand(a, T::one() - (a - n) * 2.0, K, eps);

        println!("sum trunc = {:?}, correction = {:?}", sum_trunc_dirichlet, correction);

        sum_trunc_dirichlet
            + correction / a.sqrt() * (if n.to_i32().unwrap() % 2 == 0 { -1.0 } else { 1.0 })
    }

    pub fn Z(&self, t: T, eps: f64) -> Option<T> {
        let plan = self.plan(t.to_f64().unwrap(), eps / 3.0)?;
        let K = plan.0;
        let a = (t / T::PI() / 2.0).sqrt();
        let n = a.floor();
        let p = T::one() - (a - n) * 2.0;
        // let ps = Self::gen_power_series(p, 3 * K + 1);
        // let ps = Self::gen_power_series_optimized(p, 3 * K + 1);

        println!("power series: K = {}, p = {:?}", K, p);
        let ret = self.solve(t, plan, eps);
        Some(ret)
    }
}

#[cfg(test)]
mod tests {
    use crate::*;
    use brentq::brentq;
    use num::Complex;
    use num::Zero;
    use riemann_siegel::RiemannSiegelTheta;
    use F64x2::f64x2;
    use F64x2::test_utils::*;

    #[test]
    fn rszeta() {
        type T = f64x2;
        let ctx = Context::<T>::new(100);
        let zeta_rs = RiemannSiegelZeta::new(&ctx, T::from(1.5), 20);

        let eps = 1e-10;
        let s = Complex::new(T::from(1.5), T::from(1000.0));
        let z = zeta_rs.zeta(s, eps).unwrap();
        let gt = Complex::new(
            f64x2 { hi: 0.9555445813034115, lo: 2.3151859668257056e-17 },
            f64x2 { hi: -0.09613241765159551, lo: 7.727321635166056e-19 },
        );
        println!("s = {}, zeta(s) = {}", s, z);
        assert_complex_close(z, gt, eps);

        let eps = 1e-24;
        let s = Complex::new(T::from(1.5), T::from(1000000.0));
        let z = zeta_rs.zeta(s, eps).unwrap();
        let gt = Complex::new(
            f64x2 { hi: 0.9380775942197226, lo: 1.0691936973721035e-17 },
            f64x2 { hi: 0.38907944826016855, lo: -2.7797363425251712e-18 },
        );
        println!("s = {}, zeta(s) = {}", s, z);
        assert_complex_close(z, gt, eps);
        // panic!();
    }

    #[test]
    fn test_rs_z() {
        type T = f64x2;

        let ctx = Context::<T>::new(100);
        let rs_z = RiemannSiegelZ::new(&ctx, 20);

        let eps = 1e-10;
        let t = T::from(1000.0);
        let z = rs_z.Z(t, eps).unwrap();
        let gt = f64x2 { hi: 0.9977946375215866, lo: 2.2475077894406e-17 };
        println!("t = {}, Z(t) = {}", t, z);
        assert_close(z, gt, eps);

        let eps = 1e-24;
        let t = T::from(1000000.0);
        let z = rs_z.Z(t, eps).unwrap();
        let gt = -f64x2 { hi: 2.8061338784306984, lo: 8.690744258391294e-17 };
        println!("t = {}, Z(t) = {}", t, z);
        assert_close(z, gt, eps);

        let eps = 1e-20;
        let t = f64x2 { hi: 74961.6886666817, lo: -0.000000000001845641837833102 };
        let z = rs_z.Z(t, eps).unwrap();
        let gt = f64x2 { hi: -2.502383979630123e-14, lo: -6.697961745086965e-31 };
        println!("t = {}, Z(t) = {:.e}", t, z);
        // assert_close(z, gt, eps);
        // panic!();
    }

    #[test]
    fn test_bisect() {
        type T = f64x2;
        let a = f64x2 { hi: 74961.0, lo: 0.0 };
        let b = f64x2 { hi: 74962.0, lo: 0.0 };

        let eps = 1e-20;
        let ctx = Context::<T>::new(100);
        let rs_z = RiemannSiegelZ::new(&ctx, 30);
        println!("{} {}", rs_z.Z(a, eps).unwrap(), rs_z.Z(b, eps).unwrap());

        let result = brentq(|x| rs_z.Z(x, eps).unwrap(), a, b, T::zero(), T::zero(), 30);

        // panic!();
    }

    #[test]
    fn test_rs_theta() {
        type T = f64x2;
        let t = f64x2 { hi: 1000.0, lo: 0.0 };

        let eps = 1e-30;
        let ctx = Context::<T>::new(100);
        let rs_theta = RiemannSiegelTheta::new(&ctx, 20);
        let gt = f64x2 { hi: 2034.5464280380315, lo: 7.28690383001782e-14 };
        let output = rs_theta.theta(t, eps);

        assert_close(output, gt, eps);
    }
}
