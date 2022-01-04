use crate::{brentq::EvalPoint, traits::MyReal};

pub struct Illinois<T: MyReal> {
    pub xa: EvalPoint<T>,
    pub xb: EvalPoint<T>,
    pub side: i8,
}

impl<T: MyReal> Illinois<T> {
    pub fn new(xa: EvalPoint<T>, xb: EvalPoint<T>) -> Self { Self { xa, xb, side: 0 } }

    pub fn step(&mut self, f: impl FnOnce(T) -> T) {
        let r = (self.xa.x * self.xb.f - self.xb.x * self.xa.f) / (self.xb.f - self.xa.f);
        let c = EvalPoint::new(r, f);
        if c.has_same_sign(&self.xa) {
            self.xa = c;
            if self.side == -1 {
                self.xb.f = self.xb.f * 0.5;
            }
            self.side = -1;
        } else {
            self.xb = c;
            if self.side == 1 {
                self.xa.f = self.xa.f * 0.5;
            }
            self.side = 1;
        }
    }
}
