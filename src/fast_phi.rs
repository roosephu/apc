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

pub struct LittlePhiFn {
    radius: f64,
    eps: f64,
    t: f64,
    λ: f64,
    drv: [f64; 3],
    stat: CacheStat,
}

impl LittlePhiFn {
    pub fn stat(&self) -> &CacheStat { &self.stat }

    pub fn new(λ: f64, eps: f64) -> Self {
        Self { radius: -1.0, eps, t: 0.0, λ, drv: [0.0, 0.0, 0.0], stat: CacheStat::new() }
    }

    pub fn once(&self, t: f64) -> f64 {
        rgsl::error::erfc(t.ln_1p() / self.λ / std::f64::consts::SQRT_2) / 2.0
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
        let ρ = t.ln_1p() / λ;
        let y0 = rgsl::error::erfc(ρ / std::f64::consts::SQRT_2) / 2.0;

        let coeff = (ρ * ρ / -2.0).exp() / (std::f64::consts::PI * 2.0).sqrt();
        let y1 = coeff / (1.0 + t).powi(1) * -(1.0 / λ);
        let y2 = coeff / (1.0 + t).powi(2) * (1.0 / λ + ρ / λ / λ);

        // let y3 =
        //     coeff / (1.0 + t).powi(3) * -(2.0 / λ + 3.0 * ρ / λ / λ - (1.0 - ρ * ρ) / λ / λ / λ);
        // here we compute max_t y3(t).
        // by plotting: the maximum is positive
        // exp(-ρ^2/2) ρ >= -0.63
        // exp(-ρ^2/2) ρ^2 >= 0
        let max_y3 =
            (-2.0 / λ + 1.8 / λ / λ + 1.0 / λ / λ / λ) / (std::f64::consts::PI * 2.0).sqrt();

        // Lagrange remainder: error <= |ϕ^(3)(t)| / 6 * radius^3
        let radius = (self.eps / max_y3 * 6.0).powf(1.0 / 3.0);
        // println!("radius = {}", radius);

        self.radius = radius;
        self.t = t;
        self.drv = [y0, y1, y2 / 2.0];

        y0
    }

    pub fn query(&mut self, ρ: f64) -> f64 { self.fast(ρ).unwrap_or_else(|| self.slow(ρ)) }
}
