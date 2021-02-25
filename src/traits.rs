use rug::{Complex, Float};
use rug::Assign;

pub trait New {
    fn new(prec: u32) -> Self;
}

pub trait Alloc {
    fn alloc<O: Assign<Self> + New>(self, prec: u32) -> O where Self: Sized {
        let mut ret = O::new(prec);
        ret.assign(self);
        ret
    }

    fn to<O: Assign<Self>>(self, buf: &mut O) where Self: Sized {
        buf.assign(self)
    }
}

impl<T> Alloc for T {}

impl New for Float {
    fn new(prec: u32) -> Float {
        Float::new(prec)
    }
}

impl New for Complex {
    fn new(prec: u32) -> Complex {
        Complex::new(prec)
    }
}
