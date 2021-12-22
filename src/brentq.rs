use crate::traits::MyReal;

// scipy defaults: xtol = 2e-12, rtol = 8.881784197001252e-16

pub fn brentq<T: MyReal, F>(f: F, xa: T, xb: T, xtol: T, rtol: T, iter: usize) -> Option<T>
where
    F: Fn(T) -> T,
{
    let mut xpre = xa;
    let mut xcur = xb;
    let mut xblk = T::zero();
    let mut fpre = f(xpre);
    let mut fcur = f(xcur);
    let mut fblk = T::zero();
    let mut spre = T::zero();
    let mut scur = T::zero();

    if (fpre * fcur).is_sign_positive() {
        return None;
    }
    if fpre.is_zero() {
        return Some(xpre);
    }
    if fcur.is_zero() {
        return Some(xcur);
    }

    for _ in 0..iter {
        if (fpre * fcur).is_sign_negative() {
            xblk = xpre;
            fblk = fpre;
            spre = xcur - xpre;
            scur = spre;
        }

        if fblk.abs() < fcur.abs() {
            xpre = xcur;
            xcur = xblk;
            xblk = xpre;

            fpre = fcur;
            fcur = fblk;
            fblk = fpre;
        }

        let delta = (xtol + rtol * xcur.abs()) / 2.0;
        let sbis = (xblk - xcur) / 2.0;

        if fcur.is_zero() || sbis.abs() < delta {
            return Some(xcur);
        }
        if spre.abs() > delta && fcur.abs() < fpre.abs() {
            let stry = if xpre == xblk {
                -fcur * (xcur - xpre) / (fcur - fpre)
            } else {
                let dpre = (fpre - fcur) / (xpre - xcur);
                let dblk = (fblk - fcur) / (xblk - xcur);
                -fcur * (fblk * dblk - fpre * dpre) / (dblk * dpre * (fblk - fpre))
            };
            if stry.abs() * 2.0 < spre.abs().min(sbis.abs() * 3.0 - delta) {
                spre = scur;
                scur = stry;
            } else {
                spre = sbis;
                scur = sbis;
            }
        } else {
            spre = sbis;
            scur = sbis;
        }

        xpre = xcur;
        fpre = fcur;

        if scur.abs() > delta {
            xcur += scur;
        } else {
            xcur += if sbis.is_sign_positive() { delta } else { -delta };
        }

        fcur = f(xcur);
    }
    Some(xcur)
}

#[cfg(test)]
mod tests {
    use super::brentq;

    #[test]
    fn test_brentq() {
        let f = |x| x * x * x - 0.5;
        let result = brentq(f, 0.0, 1.0, 1e-8, 1e-8, 100);
        let solution = 0.5f64.powf(1.0 / 3.0);
        assert!((result.unwrap() - solution).abs() <= 1e-8);
    }
}
