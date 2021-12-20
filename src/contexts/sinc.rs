use F64x2::f64x2;
use num::{Zero, One};

pub trait Sinc {
    fn sinc(self) -> Self;
}

impl Sinc for f64 {
    fn sinc(self) -> Self {
        rgsl::Trigonometric::sinc(&self)
    }
}

impl Sinc for f64x2 {
    fn sinc(self) -> Self {
        if self.is_zero() {
            Self::one()
        } else {
            self.sin() / self
        }
    }
}
