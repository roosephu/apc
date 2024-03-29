use num::{One, Zero};
use F64x2::f64x2;

/// Normalized sinc function: sinc(x) = sin(\pi x) / (\pi x)
pub trait Sinc {
    fn sinc(&self) -> Self;
}

impl Sinc for f64 {
    fn sinc(&self) -> Self { rgsl::Trigonometric::sinc(self) }
}

impl Sinc for f64x2 {
    fn sinc(&self) -> Self {
        if self.is_zero() {
            Self::one()
        } else {
            // TODO:
            use crate::traits::MyReal;
            let t = *self * Self::PI();
            t.sin() / t
        }
    }
}
