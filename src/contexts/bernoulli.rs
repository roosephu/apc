use F64x2::f64x2;
use super::utils::{read_data, mpf_to_f64x2, impl_from_uninit_cell};

pub trait Bernoulli {
    fn bernoulli(n: usize) -> Self;
}

impl_from_uninit_cell!(Bernoulli, bernoulli, f64);
impl_from_uninit_cell!(Bernoulli, bernoulli, f64x2);


pub fn init() {
    let data = read_data("data/bernoulli.txt", 1000)
        .expect("can't load Bernoulli numbers from `data/bernoulli.txt`");

    TABLE_f64.set(data.iter().map(|x| x.to_f64()).collect());
    TABLE_f64x2.set(data.iter().map(mpf_to_f64x2).collect());
}
