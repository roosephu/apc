use crate::cache_stat::CacheStat;

pub struct PhiFn {
    radius: f64,
    eps: f64,
    ρ: f64,
    drv: [f64; 3],
    stat: CacheStat,
}

impl PhiFn {
    pub fn stat(&self) -> &CacheStat { &self.stat }

    pub fn new(eps: f64) -> Self {
        Self { radius: -1.0, eps, ρ: 0.0, drv: [0.0, 0.0, 0.0], stat: CacheStat::new() }
    }

    pub fn once(&self, ρ: f64) -> f64 { rgsl::error::erfc(ρ / std::f64::consts::SQRT_2) / 2.0 }

    fn fast(&mut self, ρ: f64) -> Option<f64> {
        let δ = ρ - self.ρ;
        if δ < -self.radius || δ > self.radius {
            None
        } else {
            self.stat.hit();
            Some((self.drv[2] * δ + self.drv[1]) * δ + self.drv[0])
        }
    }

    fn slow(&mut self, ρ: f64) -> f64 {
        self.stat.miss();
        let y0 = rgsl::error::erfc(ρ / std::f64::consts::SQRT_2) / 2.0;
        let y1 = -(ρ * ρ / -2.0).exp() / (std::f64::consts::PI * 2.0).sqrt();
        let y2 = y1 * -ρ / 2.0;

        // Lagrange remainder: error <= |Φ^(3)(x)| / 6 * radius^3
        // |Φ^(3)(x)| <= 0.4 (by plotting...)
        let radius = (self.eps * 15.0).powf(1.0 / 3.0);

        self.radius = radius;
        self.ρ = ρ;
        self.drv = [y0, y1, y2];

        y0
    }

    pub fn query(&mut self, ρ: f64) -> f64 { self.fast(ρ).unwrap_or_else(|| self.slow(ρ)) }
}
