use F64x2::f64x2;
use super::utils::{read_data, mpf_to_f64x2, impl_from_uninit_cell};

pub trait Factorial {
    fn factorial(n: usize) -> Self;
}
    // fn binom(n: usize, m: usize) -> Self {
    //     Self::factorial(n) / Self::factorial(m) / Self::factorial(n - m)
    // }

impl_from_uninit_cell!(Factorial, factorial, f64);
impl_from_uninit_cell!(Factorial, factorial, f64x2);


pub fn init() {
    let data = read_data("data/factorial.txt", 1000)
        .expect("can't load Bernoulli numbers from `data/factorial.txt`");

    TABLE_f64.set(data.iter().map(|x| x.to_f64()).collect());
    TABLE_f64x2.set(data.iter().map(mpf_to_f64x2).collect());
}
