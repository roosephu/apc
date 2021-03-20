// scipy defaults: xtol = 2e-12, rtol = 8.881784197001252e-16

pub fn brentq<F>(f: F, xa: f64, xb: f64, xtol: f64, rtol: f64, iter: i64) -> Option<f64>
where
    F: Fn(f64) -> f64,
{
    let mut xpre = xa;
    let mut xcur = xb;
    let mut xblk = 0.0;
    let mut fpre = f(xpre);
    let mut fcur = f(xcur);
    let mut fblk = 0.0;
    let mut spre = 0.0;
    let mut scur = 0.0;

    if fpre * fcur > 0.0 {
        return None;
    }
    if fpre == 0.0 {
        return Some(xpre);
    }
    if fcur == 0.0 {
        return Some(xcur);
    }

    for _ in 0..iter {
        if fpre * fcur < 0.0 {
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

        if fcur == 0.0 || sbis.abs() < delta {
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
            if 2.0 * stry.abs() < spre.abs().min(3.0 * sbis.abs() - delta) {
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
            xcur += if sbis > 0.0 { delta } else { -delta };
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
