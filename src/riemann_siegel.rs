// Riemann-Siegel formula in the critical line.

use crate::{power_series::PowerSeries, traits::MyFloat, unchecked_cast::UncheckedCast};

fn gen_rs_power_series<T: MyFloat, const N: usize>(z: T) -> PowerSeries<T, N> {
    let mut numer = PowerSeries::<T, N>::new(
        z * z * T::FRAC_PI_2() + 3.0f64.unchecked_cast::<T>() * T::FRAC_PI_8(),
    );
    numer.data[1] = T::PI() * z;
    numer.data[2] = T::FRAC_PI_2();
    numer.cos_();

    let mut denom = PowerSeries::<T, N>::new(z * T::PI());
    denom.data[1] = T::PI();
    denom.cos_();

    numer /= &denom;
    numer
}

pub struct RiemannSiegelZeta<T, const N: usize> {
    coeff: [[T; N]; N],
}

impl<T: MyFloat, const N: usize> RiemannSiegelZeta<T, N> {
    fn new() -> Self { Self { coeff: [[T::zero(); N]; N] } }

    // Riemann-Siegel Z function
    fn Z(t: T, eps: f64) -> T {
        let a = (t / T::TAU()).sqrt();

        todo!()
    }
}
