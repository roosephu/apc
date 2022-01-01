#![allow(non_snake_case)]
#![allow(dead_code)]
use rug::{
    ops::{CompleteRound, Pow},
    Assign, Float, Integer,
};

use super::brentq::find_zeros;

fn dot(A: &Vec<Vec<Float>>, b: &Vec<Float>, prec: u32) -> Vec<Float> {
    let n = A.len();
    let m = b.len();
    let mut ret = vec![Float::new(prec); n];
    for i in 0..n {
        for j in 0..m {
            ret[i] += (&A[i][j] * &b[j]).complete(prec);
        }
    }
    ret
}

/// Given $A$ in $\mathbb{R}^{n \times n}$ and $b$ in $R^{n}$, compute x such that $Ax = b$
fn solve(A: &Vec<Vec<Float>>, b: &Vec<Float>, prec: u32) -> Vec<Float> {
    let mut A = A.clone();
    let mut b = b.clone();

    let n = b.len();
    assert!(A.len() == n);
    for row in &A {
        assert!(row.len() == n);
    }

    let mut idx = (0..n).collect::<Vec<_>>();
    for i_ in 0..n {
        let mut ok = false;
        for j_ in i_..n {
            if !A[idx[i_]][idx[j_]].is_zero() {
                idx.swap(i_, j_);
                ok = true;
                break;
            }
        }
        assert!(ok, "singular matrix!");
        let i = idx[i_];
        for &j in &idx[i_ + 1..n] {
            let coeff = (&A[j][i] / &A[i][i]).complete(prec);
            for k in i + 1..n {
                let diff = (&coeff * &A[i][k]).complete(prec);
                A[j][k] -= diff;
            }
            b[j] = &b[j] - &b[i] * coeff;
            // A[j][i].assign(Float::new(prec));
        }
    }
    for &i in idx.iter().rev() {
        b[i] /= &A[i][i];
        for &j in &idx[0..i] {
            A[j][i] *= &b[i];
            b[j] -= &A[j][i];
        }
    }

    b
}

pub fn minimize(
    f: impl Fn(&Float) -> Float,
    xa: &Float,
    xb: &Float,
    prec: u32,
    eps: f64,
    iter: usize,
) -> Option<Float> {
    assert!(xa <= xb);

    let ratio = (3.0f64 - Float::with_val(prec, 5.0f64).sqrt()) * 0.5f64;
    let mut fulc = xa + ratio.clone() * (xb - xa).complete(prec);
    let mut nfc = fulc.clone();
    let mut xf = fulc.clone();
    let mut xa = xa.clone();
    let mut xb = xb.clone();
    let mut rat = Float::new(prec);
    let mut e = Float::new(prec);
    let mut x = xf.clone();
    let mut fx = f(&x);
    let (mut ffulc, mut fnfc) = (fx.clone(), fx.clone());
    let mut xm = (&xa + &xb).complete(prec) * 0.5f64;
    let mut tol1 = xf.to_f64().abs() * eps + eps / 3.0;
    let mut tol2 = tol1 * 2.0;
    let mut fu;

    for t in 0..iter {
        if (&xf - &xm).complete(prec).abs() <= tol2 - (&xb - &xa).complete(prec) * 0.5f64 {
            // println!("early stop {t}");
            return Some(x);
        }
        let mut use_golden = true;
        if e.clone().abs() > tol1 {
            let mut r = (&xf - &nfc).complete(prec) * (&fx - &ffulc).complete(prec);
            let q = (&xf - &fulc).complete(prec) * (&fx - &fnfc).complete(prec);
            let mut p = (&xf - &fulc).complete(prec) * &q - (&xf - &nfc).complete(prec) * &r;
            let mut q = (q - &r) * 2.0f64;
            if q > 0.0f64 {
                p = -p;
            } else {
                q = -q;
            }

            r.assign(&e);
            e.assign(&rat);
            if p.clone().abs() < (q.clone() * &r * 0.5f64).abs()
                && p > &q * (xa.clone() - &xf)
                && p < &q * (xb.clone() - &xf)
            {
                if &x - xa.clone() < tol2 || xb.clone() - &x < tol2 {
                    rat.assign(Float::with_val(prec, if xm >= xf { tol1 } else { -tol1 }));
                } else {
                    rat.assign(&p / &q);
                }
                use_golden = false;
            }
        }

        if use_golden {
            if xf >= xm {
                e.assign(&xa - &xf);
            } else {
                e.assign(&xb - &xf);
            }
            rat.assign(&ratio * &e);
        }

        let sign = if rat >= 0.0f64 { 1.0 } else { -1.0f64 };
        x = &xf + sign * rat.clone().abs().max(&Float::with_val(prec, tol1));
        fu = f(&x);

        if fu <= fx {
            if x >= xf {
                xa = xf.clone();
            } else {
                xb = xf.clone();
            }
            (fulc, ffulc) = (nfc, fnfc);
            (nfc, fnfc) = (xf, fx);
            (xf, fx) = (x.clone(), fu);
        } else {
            if x < xf {
                xa = x.clone();
            } else {
                xb = x.clone();
            }
            if fu <= fnfc || nfc == xf {
                (fulc, ffulc) = (nfc, fnfc);
                (nfc, fnfc) = (x.clone(), fu);
            } else if fu <= ffulc || fulc == xf || fulc == nfc {
                (fulc, ffulc) = (x.clone(), fu);
            }
        }

        xm.assign(&xa + &xb);
        xm /= 2.0f64;
        tol1 = xf.to_f64().abs() * eps + eps / 3.0;
        tol2 = tol1 * 2.0;
    }
    Some(x)
}

fn gen_coeff_remez(
    f: impl Fn(&Float) -> Float,
    exponents: &[usize],
    xa: f64,
    xb: f64,
    prec: u32,
) -> Vec<Float> {
    let pi = Float::with_val(prec, rug::float::Constant::Pi);
    let n = exponents.len();
    let mut points = vec![Float::new(prec); n + 1];
    for i in 0..=n {
        points[n - i] =
            (Float::with_val(prec, 2 * i + 1) / (2 * (n + 1)) as f64 * &pi).cos() * (xb - xa) / 2.0
                + (xa + xb) / 2.0;
        // println!("{:.6e}", points[n - i].to_f64());
    }
    let eps = 1e-200;
    for t in 0..100 {
        let mut A = vec![vec![Float::new(prec); n + 1]; n + 1];
        let b = points.iter().map(&f).collect::<Vec<_>>();
        for i in 0..n {
            for j in 0..=n {
                A[j][i] = points[j].clone().pow(Integer::from(exponents[i]));
            }
        }
        for i in 0..=n {
            A[i][n] = Float::with_val(prec, if i % 2 == 0 { 1i32 } else { -1i32 }) * &b[i];
        }
        // for i in 0..=n {
        //     for j in 0..=n {
        //         print!("{:.6e} ", A[i][j].to_f64());
        //     }
        //     println!(";");
        // }

        let coeffs = solve(&A, &b, prec);

        let poly = |x: &Float| {
            let mut ret = Float::new(prec);
            for i in 0..n {
                ret += &coeffs[i] * x.clone().pow(Integer::from(exponents[i]));
            }
            ret
        };
        // absolute error
        // let err = |x: Float| poly(x.clone()) - f(&x);

        // relative error
        let err = |x: &Float| {
            let my = poly(x);
            let gt = f(x);
            if gt == 0.0f64 {
                assert!(my == 0.0f64);
                Float::new(prec)
            } else {
                (my - &gt) / gt
            }
        };

        let mut roots = vec![Float::new(prec); n + 2];
        roots[0] = Float::with_val(prec, xa);
        for i in 1..=n {
            roots[i] = find_zeros(&err, &points[i - 1], &points[i], prec, eps, eps, 1000).unwrap();
        }
        roots[n + 1] = Float::with_val(prec, xb);
        for i in 0..=n {
            points[i] =
                minimize(|x| -err(x).abs(), &roots[i], &roots[i + 1], prec, eps, 1000).unwrap();
        }
        let errors = points.iter().map(|x| err(x).to_f64().abs()).collect::<Vec<_>>();

        let max_err = errors.iter().copied().reduce(f64::max).unwrap();
        let min_err = errors.iter().copied().reduce(f64::min).unwrap();
        let ratio = max_err / min_err;

        println!("iter {t}, err = {max_err:.6e} {min_err:.6e}, ratio = {ratio:.6e}");

        if ratio < 1.005 {
            break;
        }
    }
    unreachable!("Remez algorithm fails to converge!")
}

#[cfg(test)]
mod tests {
    use rug::{rand::RandState, Assign, Integer};

    use super::*;

    #[test]
    fn test_lineq_solver() {
        let n = 100;
        let prec = 1000;
        let mut rng = RandState::new();
        rng.seed(&Integer::from(12345u64));
        let mut A = vec![vec![Float::new(prec); n]; n];
        let mut b = vec![Float::new(prec); n];
        for i in 0..n {
            for j in 0..n {
                A[i][j].assign(Float::random_bits(&mut rng));
            }
        }
        for i in 0..n {
            b[i].assign(Float::random_bits(&mut rng));
        }
        let x = solve(&A, &b, prec);
        let b_ = dot(&A, &x, prec);
        // dbg!(A, &b, x, &b_);
        for i in 0..n {
            let diff = (&b[i] - &b_[i]).complete(prec).to_f64().abs();
            assert!(diff < 1e-290, "diff = {diff:.6e}");
        }
    }

    #[test]
    fn test_minimize() {
        let n = 100;
        let prec = 500;
        let f = |x: &Float| x.clone() * x * x - x;
        let eps = 1e-70;
        let iter = 1000;
        let result =
            minimize(f, &Float::with_val(prec, 0.0), &Float::with_val(prec, 1.0), prec, eps, iter);
        let minimizer = Float::with_val(prec, 3.0).recip_sqrt();
        let diff = (result.unwrap() - minimizer).to_f64();
        println!("diff = {diff:.6e}");
        assert!(diff.abs() < eps);
    }

    #[test]
    fn test_remez() {
        let prec = 1000;
        let exps = (2..=50).step_by(2).collect::<Vec<_>>();
        println!("exps = {exps:?}, len = {}", exps.len());
        let f = |x: &Float| {
            let y = x.clone().exp();
            x * (y.clone() + 1.0f64) / (y.clone() - 1.0f64) - 2.0
        };
        let coeffs = gen_coeff_remez(f, &exps, 0.0, 2.0f64.ln() / 2.0, prec);
        assert!(false);
    }
}
