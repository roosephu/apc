use std::marker::PhantomData;

pub trait Interpolation<T, U> {
    fn new(x: T, eps: f64) -> Self;
    fn query(&self, x: T) -> Option<U>;
}

pub struct AdaptiveInterp<T: Copy, U, I: Interpolation<T, U>> {
    interpolator: Option<I>,
    eps: f64,
    hit: usize,
    miss: usize,
    _marker: PhantomData<(T, U)>,
}

impl<T: Copy, U, I: Interpolation<T, U>> AdaptiveInterp<T, U, I> {
    pub fn new(eps: f64) -> Self {
        Self { interpolator: None, eps, _marker: PhantomData, hit: 0, miss: 0 }
    }

    #[inline(never)] // for flamegraph
    pub fn query(&mut self, x: T) -> U {
        if let Some(ref interpolator) = self.interpolator {
            if let Some(y) = interpolator.query(x) {
                self.hit += 1;
                return y;
            }
        }
        self.miss += 1;
        let interpolator = I::new(x, self.eps);
        let y = interpolator.query(x).unwrap();
        self.interpolator = Some(interpolator);
        y
    }

    pub fn stat(&self) -> (usize, usize) { (self.hit, self.miss) }
}

pub struct FastPhiInterp {
    radius: f64,
    eps: f64,
    x: f64,
    drv: [f64; 3],
}

impl Interpolation<f64, f64> for FastPhiInterp {
    fn new(x: f64, eps: f64) -> Self {
        let y0 = rgsl::error::erfc(x / std::f64::consts::SQRT_2) / 2.0;
        let y1 = -(x * x / -2.0).exp() / (std::f64::consts::PI * 2.0).sqrt();
        let y2 = y1 * -x / 2.0;

        // Lagrange remainder: error <= |Φ^(3)(x)| / 6 * radius^3
        // |Φ^(3)(x)| <= 0.4 (by plotting...)
        let radius = (eps * 15.0).powf(1.0 / 3.0);
        Self { radius, eps, x, drv: [y0, y1, y2] }
    }

    fn query(&self, x: f64) -> Option<f64> {
        let δ = x - self.x;
        if δ < -self.radius || δ > self.radius {
            None
        } else {
            let y = (self.drv[2] * δ + self.drv[1]) * δ + self.drv[0];

            #[cfg(feature = "verify_fast_phi")]
            {
                let gt = rgsl::error::erfc(x / std::f64::consts::SQRT_2) / 2.0;
                assert!((y - gt).abs() <= self.eps, "y = {}, value = {}", y, gt);
            }
            Some(y)
        }
    }
}

pub type FastPhi = AdaptiveInterp<f64, f64, FastPhiInterp>;
