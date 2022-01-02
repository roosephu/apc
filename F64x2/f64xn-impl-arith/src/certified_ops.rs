
use proc_macro2::TokenStream;
use quote::quote;

pub fn gen_certified_add() -> TokenStream {
    quote! {
        /// Algorithm 6 in [5, 6]
        impl f64x2 {
            fn certified_add([a0, a1]: [f64; 2], [b0, b1]: [f64; 2]) -> [f64; 2] {
                let (s0, s1) = two_add(a0, b0);
                let (t0, t1) = two_add(a1, b1);
                let c = s1 + t0;
                let (v0, v1) = two_add_fast(s0, c);
                let w = t1 + v1;
                to_array(two_add_fast(v0, w))
            }
        }
    }
}

pub fn gen_certified_add_f64() -> TokenStream {
    quote! {
        impl f64x2 {
            fn certified_add_f64([a0, a1]: [f64; 2], b: f64) -> [f64; 2] {
                let (s0, s1) = two_add(a0, y);
                let v = a1 + s1;
                to_array(two_add_fast(s0, v))
                [v, w]
            }
        }
    }
}

pub fn gen_certified_sub() -> TokenStream {
    quote! {
        impl f64x2 {
            fn certified_sub_f64(a0: f64, a1: f64, b0: f64, b1: f64) -> (f64, f64) {
                let x = self;
                let s = two_sub(x.hi, y.hi);
                let t = two_sub(x.lo, y.lo);
                let c = s.1 + t.0;
                let v = two_add_fast(s.0, c);
                let w = t.1 + v.1;
                Self::from(two_add_fast(v.0, w))
            }
        }
    }
}

pub fn gen_certified_sub_f64() -> TokenStream {

}
