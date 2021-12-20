use std::cell::Cell;
use std::mem::MaybeUninit;

pub struct UninitCell<T>(Cell<MaybeUninit<T>>);

impl<T> UninitCell<T> {
    pub const fn uninit() -> Self { Self(Cell::new(MaybeUninit::uninit())) }

    #[inline]
    pub fn set(&self, x: T) { self.0.set(MaybeUninit::new(x)); }

    #[inline]
    pub fn get(&self) -> &'static mut T { unsafe { (&mut *self.0.as_ptr()).assume_init_mut() } }
}

impl<T> std::ops::Deref for UninitCell<T> {
    type Target = T;
    fn deref(&self) -> &'static Self::Target { self.get() }
}

unsafe impl<T> Sync for UninitCell<T> {}

#[macro_export]
macro_rules! make_mut_static {
    (static $N: ident : $T: ty = $e: expr; ) => {
        #[inline]
        fn $N() -> &'static mut $T {
            use std::sync::Once;
            use std::mem::MaybeUninit;
            static once: Once = Once::new();
            static mut value: MaybeUninit<$T> = MaybeUninit::uninit();
            once.call_once(|| unsafe { value = MaybeUninit::new($e) });
            unsafe { &mut *value.as_mut_ptr() }  // or unsafe { value.assume_init_mut() }
        }
    };
}

pub use make_mut_static;
