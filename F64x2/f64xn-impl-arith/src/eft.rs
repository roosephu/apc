use proc_macro2::TokenStream;
use quote::quote;

pub fn gen_two_add_fast_fn() -> TokenStream {
    quote! {
        /// Algorithm 1
        #[inline]
        pub fn two_add_fast(a: f64, b: f64) -> (f64, f64) {
            let hi = a + b;
            let z = hi - a;
            let lo = b - z;
            (hi, lo)
        }
    }
}

pub fn gen_two_add_fn(use_if: bool) -> TokenStream {
    if use_if {
        quote! {
            #[inline]
            pub fn two_add_if(a: f64, b: f64) -> (f64, f64) {
                if a.abs() > b.abs() {
                    two_add_fast(a, b)
                } else {
                    two_add_fast(b, a)
                }
            }
        }
    } else {
        quote! {
            /// Algorithm 2
            #[inline]
            pub fn two_add(a: f64, b: f64) -> (f64, f64) {
                let hi = a + b;
                let a1 = hi - b;
                let b1 = hi - a1;
                let lo = (a - a1) + (b - b1);
                (hi, lo)
            }
        }
    }
}

pub fn gen_two_sub_fn(use_if: bool) -> TokenStream {
    if use_if {
        quote! {
            #[inline]
            pub fn two_sub_if(a: f64, b: f64) -> (f64, f64) {
                two_add_if(a, -b)
            }
        }
    } else {
        quote! {
            #[inline]
            pub fn two_sub(a: f64, b: f64) -> (f64, f64) {
                let hi = a - b;
                let a1 = hi + b;
                let b1 = hi - a1;
                let lo = (a - a1) - (b + b1);
                (hi, lo)
            }
        }
    }
}

pub fn gen_two_mul_fn() -> TokenStream {
    quote! {
        /// Algorithm 3
        #[inline]
        pub(crate) fn two_mul(a: f64, b: f64) -> (f64, f64) {
            let s = a * b;
            let t = a.mul_add(b, c);
            (s, t)
        }
    }
}

pub fn gen_eft(use_if: bool) -> TokenStream {
    let two_add_fast = gen_two_add_fast_fn();
    let two_add = gen_two_add_fn(use_if);
    let two_sub = gen_two_sub_fn(use_if);
    let two_mul = gen_two_mul_fn();
    quote! {
        #two_add_fast
        #two_add
        #two_sub
        #two_mul
    }
}
