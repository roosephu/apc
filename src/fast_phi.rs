use crate::cache_stat::CacheStat;

pub struct LittlePhiSum {
    radius: i64,
    eps: f64,
    t: i64,
    x: f64,
    λ: f64,
    drv: [f64; 3],
    sum: f64,
    sum0: i64,
    sum1: i64,
    sum2: i64,
    stat: CacheStat,
}

impl LittlePhiSum {
    pub fn stat(&self) -> &CacheStat { &self.stat }

    pub fn sum(self) -> f64 { self.sum }

    pub fn new(λ: f64, x: f64, eps: f64) -> Self {
        Self {
            radius: -1,
            eps,
            t: 0,
            x,
            λ,
            drv: [0.0, 0.0, 0.0],
            stat: CacheStat::new(),
            sum0: 0,
            sum1: 0,
            sum2: 0,
            sum: 0.0,
        }
    }

    pub fn once(&self, t: i64) -> f64 {
        rgsl::error::erfc((t as f64 / self.x).ln_1p() / self.λ / std::f64::consts::SQRT_2) / 2.0
    }

    fn fast(&mut self, t: i64) -> Option<()> {
        let δ = t - self.t;
        if δ < -self.radius || δ > self.radius {
            None
        } else {
            self.stat.hit();
            self.sum0 += 1;
            self.sum1 += δ;
            self.sum2 += δ * δ;
            Some(())
        }
    }

    pub fn flush(&mut self) {
        let s0 = self.sum0 as f64;
        let s1 = self.sum1 as f64;
        let s2 = self.sum2 as f64;
        self.sum += s0 * self.drv[0] + s1 * self.drv[1] + s2 * self.drv[2];

        self.sum0 = 0;
        self.sum1 = 0;
        self.sum2 = 0;
    }

    fn slow(&mut self, t: i64) {
        self.stat.miss();
        self.flush();

        let λ = self.λ;
        let x = self.x;
        let t_ = t as f64;

        let ρ = (t_ / x).ln_1p() / λ;
        let y0 = rgsl::error::erfc(ρ / std::f64::consts::SQRT_2) / 2.0;

        let coeff = (ρ * ρ / -2.0).exp() / (std::f64::consts::PI * 2.0).sqrt();
        let y1 = coeff / (x + t_).powi(1) * -(1.0 / λ);
        let y2 = coeff / (x + t_).powi(2) * (1.0 / λ + ρ / λ / λ);

        // See Readme.md for derivation.
        let max_y3 = (-2.0 / λ + 1.8 / λ / λ + 1.0 / λ / λ / λ)
            / (std::f64::consts::PI * 2.0).sqrt()
            / x.powi(3);
        let radius = (self.eps / max_y3 * 6.0).powf(1.0 / 3.0) as i64;
        let radius = std::cmp::min(radius, 1 << 20);

        self.radius = radius;
        self.t = t;
        self.drv = [y0, y1, y2 / 2.0];

        self.sum0 = 1;
    }

    pub fn add(&mut self, t: i64) { self.fast(t).unwrap_or_else(|| self.slow(t)) }
}
