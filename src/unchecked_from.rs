pub trait UncheckedFrom<S> {
    fn unchecked_from(_: S) -> Self;
}

pub trait UncheckedInto<T> {
    fn unchecked_into(self) -> T;
}

pub trait UncheckedCast {
    fn unchecked_cast<T>(self) -> T
    where
        Self: UncheckedInto<T> + Sized,
    {
        self.unchecked_into()
    }
}

impl<S, T: UncheckedFrom<S>> UncheckedInto<T> for S {
    fn unchecked_into(self) -> T { T::unchecked_from(self) }
}

macro_rules! impl_primitive_uncheck_from {
    ($S:ty => $T:ty) => {
        impl UncheckedFrom<$S> for $T {
            fn unchecked_from(x: $S) -> $T {
                x as $T
            }
        }
    };
    ($S:ty => $($T:ty),*) => {
        $(impl_primitive_uncheck_from!($S => $T); )*
    };
}

impl_primitive_uncheck_from!(usize => usize, u8, u16, u32, u64, u128, isize, i8, i16, i32, i64, i128, f32, f64);
impl_primitive_uncheck_from!(u8 => usize, u8, u16, u32, u64, u128, isize, i8, i16, i32, i64, i128, f32, f64);
impl_primitive_uncheck_from!(u16 => usize, u8, u16, u32, u64, u128, isize, i8, i16, i32, i64, i128, f32, f64);
impl_primitive_uncheck_from!(u32=> usize, u8, u16, u32, u64, u128, isize, i8, i16, i32, i64, i128, f32, f64);
impl_primitive_uncheck_from!(u64 => usize, u8, u16, u32, u64, u128, isize, i8, i16, i32, i64, i128, f32, f64);
impl_primitive_uncheck_from!(u128 => usize, u8, u16, u32, u64, u128, isize, i8, i16, i32, i64, i128, f32, f64);
impl_primitive_uncheck_from!(isize => usize, u8, u16, u32, u64, u128, isize, i8, i16, i32, i64, i128, f32, f64);
impl_primitive_uncheck_from!(i8 => usize, u8, u16, u32, u64, u128, isize, i8, i16, i32, i64, i128, f32, f64);
impl_primitive_uncheck_from!(i16 => usize, u8, u16, u32, u64, u128, isize, i8, i16, i32, i64, i128, f32, f64);
impl_primitive_uncheck_from!(i32 => usize, u8, u16, u32, u64, u128, isize, i8, i16, i32, i64, i128, f32, f64);
impl_primitive_uncheck_from!(i64 => usize, u8, u16, u32, u64, u128, isize, i8, i16, i32, i64, i128, f32, f64);
impl_primitive_uncheck_from!(i128 => usize, u8, u16, u32, u64, u128, isize, i8, i16, i32, i64, i128, f32, f64);
impl_primitive_uncheck_from!(f32 => usize, u8, u16, u32, u64, u128, isize, i8, i16, i32, i64, i128, f32, f64);
impl_primitive_uncheck_from!(f64 => usize, u8, u16, u32, u64, u128, isize, i8, i16, i32, i64, i128, f32, f64);

macro_rules! impl_primitive_uncheck_cast {
    ($($T:ty),*) => {
        $(
            impl UncheckedCast for $T {}
        )*
    };
}

impl_primitive_uncheck_cast!(
    usize, u8, u16, u32, u64, u128, isize, i8, i16, i32, i64, i128, f32, f64
);
