use num::{Complex};
use num_traits::{Float, FloatConst, NumAssignOps, NumOps, One, Zero, FromPrimitive};
use std::{fmt::{Debug, Display}, ops::{Add, Sub, Div, Mul, Neg}};

pub trait GenericFloat = Float + FloatConst + Zero + One + NumOps + NumAssignOps + FromPrimitive +
    Display + Debug + Default +
    Neg<Output=Self> +
    Add<Complex<Self>, Output=Complex<Self>> +
    Sub<Complex<Self>, Output=Complex<Self>> +
    Mul<Complex<Self>, Output=Complex<Self>> +
    Div<Complex<Self>, Output=Complex<Self>> +
    'static;
