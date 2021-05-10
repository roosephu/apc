use crate::cache_stat::CacheStat;

pub struct LittlePhiFn {
    radius: f64,
    eps: f64,
    t: f64,
    x: f64,
    λ: f64,
    drv: [f64; 3],
    stat: CacheStat,
}

impl LittlePhiFn {
    pub fn stat(&self) -> &CacheStat { &self.stat }

    pub fn new(λ: f64, x: f64, eps: f64) -> Self {
        Self { radius: -1e100, eps, t: 0.0, x, λ, drv: [0.0, 0.0, 0.0], stat: CacheStat::new() }
    }

    pub fn once(&self, t: f64) -> f64 {
        rgsl::error::erfc((t / self.x).ln_1p() / self.λ / std::f64::consts::SQRT_2) / 2.0
    }

    fn fast(&mut self, t: f64) -> Option<f64> {
        let δ = t - self.t;
        if δ < -self.radius || δ > self.radius {
            None
        } else {
            self.stat.hit();
            Some((self.drv[2] * δ + self.drv[1]) * δ + self.drv[0])
        }
    }

    fn slow(&mut self, t: f64) -> f64 {
        self.stat.miss();
        let λ = self.λ;
        let x = self.x;

        let ρ = (t / x).ln_1p() / λ;
        let y0 = rgsl::error::erfc(ρ / std::f64::consts::SQRT_2) / 2.0;

        let coeff = (ρ * ρ / -2.0).exp() / (std::f64::consts::PI * 2.0).sqrt();
        let y1 = coeff / (x + t).powi(1) * -(1.0 / λ);
        let y2 = coeff / (x + t).powi(2) * (1.0 / λ + ρ / λ / λ);

        // See Readme.md for derivation.
        let max_y3 = (-2.0 / λ + 1.8 / λ / λ + 1.0 / λ / λ / λ)
            / (std::f64::consts::PI * 2.0).sqrt()
            / x.powi(3);
        let radius = (self.eps / max_y3 * 6.0).powf(1.0 / 3.0);

        self.radius = radius;
        self.t = t;
        self.drv = [y0, y1, y2 / 2.0];

        y0
    }

    pub fn query(&mut self, ρ: f64) -> f64 { self.fast(ρ).unwrap_or_else(|| self.slow(ρ)) }
}
