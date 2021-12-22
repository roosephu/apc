use crate::{traits::{MyReal, GabckeExpansion}, context::Context, rs_theta::RiemannSiegelTheta};


/// Gabcke's
pub struct RiemannSiegelZ<'a, T> {
    ctx: &'a Context<T>,
    theta: RiemannSiegelTheta<'a, T>,
    K: usize,
}

impl<'a, T: MyReal + GabckeExpansion> RiemannSiegelZ<'a, T> {
    pub fn new(ctx: &'a Context<T>, K: usize) -> Self {
        Self { ctx, K, theta: RiemannSiegelTheta::new(ctx, K) }
    }
}

impl<T: MyReal + GabckeExpansion> RiemannSiegelZ<'_, T> {
    fn plan(&self, t: f64, eps: f64) -> Option<(usize, Option<T>)> {
        if t <= 200.0 {
            return None;
        }

        const COEFFS: [f64; 10] =
            [0.127, 0.053, 0.011, 0.031, 0.017, 0.061, 0.661, 9.2, 130.0, 1837.0];
        let mut pow = t.powf(-0.75);
        let sqrt_t = t.sqrt();
        let coeff_bound = eps / T::epsilon().fp() / 10.0f64;

        for (k, c) in COEFFS.iter().enumerate() {
            if self.coeffs[k][0].fp() > coeff_bound {
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
        let a = (t / T::PI() / 2.0).sqrt();
        let n = a.floor();
        let (K, plan_sum_trunc_dirichlet) = plan;
        // let z = Complex::new(T::mp(0.25), t * 0.5);
        // let theta = self.ctx.loggamma(z, eps).im - t * 0.5 * T::PI().ln();

        let mut sum_trunc_dirichlet;
        match plan_sum_trunc_dirichlet {
            Some(x) => {
                sum_trunc_dirichlet = x;
            }
            None => {
                let theta = self.theta.theta(t, eps);
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
        let plan = self.plan(t.fp(), eps / 3.0)?;
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
