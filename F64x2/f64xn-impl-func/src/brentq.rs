use rug::{ops::CompleteRound, Assign, Float};

// scipy defaults: xtol = 2e-12, rtol = 8.881784197001252e-16
#[allow(dead_code)]
#[allow(clippy::float_cmp)]
pub fn find_zeros(
    f: impl Fn(&Float) -> Float,
    xa: &Float,
    xb: &Float,
    prec: u32,
    xtol: f64,
    rtol: f64,
    iter: usize,
) -> Option<Float> {
    let mut xpre = xa.clone();
    let mut xcur = xb.clone();
    let mut xblk = Float::new(prec);
    let mut fpre = f(&xpre);
    let mut fcur = f(&xcur);
    let mut fblk = Float::new(prec);
    let mut spre = Float::new(prec);
    let mut scur = Float::new(prec);
    let mut sbis = Float::new(prec);
    let mut stry = Float::new(prec);

    if fpre.to_f64() * fcur.to_f64() > 0.0 {
        return None;
    }
    if fpre.is_zero() {
        return Some(xpre);
    }
    if fcur.is_zero() {
        return Some(xcur);
    }

    for _ in 0..iter {
        if fpre.to_f64() * fcur.to_f64() < 0.0 {
            xblk.assign(&xpre);
            fblk.assign(&fpre);
            spre.assign(&xcur - &xpre);
            scur.assign(&spre);
        }

        if fblk.clone().abs() < fcur.clone().abs() {
            xpre.assign(&xcur);
            xcur.assign(&xblk);
            xblk.assign(&xpre);

            fpre.assign(&fcur);
            fcur.assign(&fblk);
            fblk.assign(&fpre);
        }

        let delta = (xtol + rtol * xcur.clone().abs()) / 2.0f64;
        sbis.assign(&xblk - &xcur);
        sbis /= 2.0f64;
        // println!("iter {t}, f(cur) = {}");

        if fcur.is_zero() || sbis.to_f64().abs() < delta {
            // println!("finish at iter {t}");
            return Some(xcur);
        }

        let mut use_bis = true;
        if spre.clone().abs() > delta && fcur.clone().abs() < fpre.clone().abs() {
            stry.assign(if xpre == xblk {
                // interpolate
                -(xcur.clone() - xpre.clone()) / (fcur.clone() - fpre.clone()) * &fcur
            } else {
                // extrapolate
                let dpre = (&fpre - &fcur).complete(prec) / (&xpre - &xcur).complete(prec);
                let dblk = (&fblk - &fcur).complete(prec) / (&xblk - &xcur).complete(prec);
                -((&fblk * &dblk).complete(prec) - (&fpre * &dpre).complete(prec))
                    / ((&dblk * &dpre).complete(prec) * (&fblk - &fpre).complete(prec))
                    * &fcur
            });
            if 2.0 * stry.to_f64().abs()
                < spre.to_f64().abs().min(3.0 * sbis.to_f64().abs() - delta.to_f64())
            {
                use_bis = false;
            }
        }
        if use_bis {
            spre.assign(&sbis);
            scur.assign(&sbis);
        } else {
            spre.assign(&scur);
            scur.assign(&stry);
        }

        xpre.assign(&xcur);
        fpre.assign(&fcur);

        if scur.to_f64().abs() > delta {
            xcur += &scur;
        } else {
            xcur += if sbis.is_sign_positive() { delta } else { -delta };
        }

        fcur = f(&xcur);
    }
    Some(xcur)
}

// pub fn find_zeros<F>(f: F, xa: f64, xb: f64, xtol: f64, rtol: f64, iter: usize) -> Option<Float>
// where
//     F: Fn(f64) -> f64,
// {
//     let ratio = 0.5 + 1.25.sqrt();  // ϕ = (1 + \sqrt(5)) / 2
//     for _ in 0..iter {

//     }
// }

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_brentq() {
        let f = |x: &Float| x.clone() * x * x - 0.5f64;
        let prec = 1000u32;
        let eps = 1e-100;
        let result = find_zeros(
            f,
            &Float::with_val(prec, 0.0f64),
            &Float::with_val(prec, 1.0f64),
            prec,
            eps,
            eps,
            1000,
        );
        // let solution = 0.5f64.powf(1.0 / 3.0);
        let value = f(&result.unwrap()).to_f64();
        println!("value = {value:.6e}");
        assert!(value <= eps);
    }
}
