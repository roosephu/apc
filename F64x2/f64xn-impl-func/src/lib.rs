#![allow(dead_code)]
mod brentq;
mod remez;

use proc_macro2::{TokenStream, Ident};
use quote::{format_ident, quote};
use rug::{float::Constant, Float};
use std::ops::SubAssign;
use syn::{parse_macro_input, LitInt};

fn gen_float_consts(n: usize) -> TokenStream {
    let prec = n as u32 * 64;
    let constants = [
        ("PI", Float::with_val(prec, Constant::Pi)),
        ("E", Float::with_val(prec, 1.0f64).exp()),
        ("LOG2_E", Float::with_val(prec, 1.0f64).exp().log2()),
        ("LN_10", Float::with_val(prec, 10.0f64).ln()),
        ("LN_2", Float::with_val(prec, 2.0f64).ln()),
        ("LOG10_2", Float::with_val(prec, 2.0f64).log10()),
        ("LOG10_E", Float::with_val(prec, 1.0f64).exp().log10()),
        ("LOG2_10", Float::with_val(prec, 10.0f64).log2()),
        ("FRAC_1_PI", Float::with_val(prec, Constant::Pi).recip()),
        ("FRAC_1_SQRT_2", Float::with_val(prec, 2.0f64).recip_sqrt()),
        ("FRAC_2_PI", Float::with_val(prec, Constant::Pi).recip() * 2i32),
        ("FRAC_2_SQRT_PI", Float::with_val(prec, Constant::Pi).recip_sqrt() * 2i32),
        ("FRAC_PI_2", Float::with_val(prec, Constant::Pi) / 2i32),
        ("FRAC_PI_3", Float::with_val(prec, Constant::Pi) / 3i32),
        ("FRAC_PI_4", Float::with_val(prec, Constant::Pi) / 4i32),
        ("FRAC_PI_6", Float::with_val(prec, Constant::Pi) / 6i32),
        ("FRAC_PI_8", Float::with_val(prec, Constant::Pi) / 8i32),
        ("SQRT_2", Float::with_val(prec, 2.0f64).sqrt()),
        ("TAU", Float::with_val(prec, Constant::Pi) * 2i32),
    ];
    let mut def = vec![];
    let mut implement = vec![];
    for (name, mut value) in constants {
        let mut data = vec![0.0; n];
        for i in 0..n {
            data[i] = value.to_f64();
            value.sub_assign(data[i]);
        }
        let name = format_ident!("{}", name);
        def.push(quote! { const #name: Self = Self { data: [#(#data),*] }; });
        implement.push(quote! { #[inline] fn #name() -> Self { Self::#name } });
    }

    let prec = 53 * n as u32;
    let epsilon = 0.5f64.powi(prec as i32);

    quote! {
        impl f64x::<#n> {
            const EPSILON: f64 = #epsilon;
            const PREC: u32 = #prec;
            #(#def)*
        }

        impl num::traits::FloatConst for f64x::<#n> {
            #(#implement)*
        }
    }
}

fn impl_f64xn_consts(n: usize) -> TokenStream {
    let float_const = gen_float_consts(n);
    quote! {
        #float_const
    }
}

fn gen_poly_eval(x: &Ident, n: usize, coeffs: &Ident) -> TokenStream {
    let y = format_ident!("y");
    let mut code = vec![];
    code.push(quote! { let mut #y = #x; });

    for i in (0..n).rev() {
        code.push(quote! { #y = #y * #x + #coeffs[#i]; });
    }

    quote! {
        {
            #(#code)*
            #y
        }
    }
}

#[proc_macro]
pub fn f64xn_impl_func(tokens: proc_macro::TokenStream) -> proc_macro::TokenStream {
    let input = parse_macro_input!(tokens as LitInt);
    let n = input.base10_parse::<usize>().unwrap();

    let expanded = impl_f64xn_consts(n);
    expanded.into()
}
