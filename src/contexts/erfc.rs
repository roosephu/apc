use F64x2::f64x2;

pub trait Erfc {
    fn erfc(&self, eps: f64) -> Self;
}

impl Erfc for f64 {
    fn erfc(&self, _: f64) -> Self { rgsl::error::erfc(*self) }
}

impl Erfc for f64x2 {
    fn erfc(&self, eps: f64) -> Self { self.erfc_eps(eps) }
}
