#![allow(dead_code)]

extern crate proc_macro;

use std::{borrow::Borrow, collections::VecDeque, ops::SubAssign};

use proc_macro2::TokenStream;
use quote::{format_ident, quote};
use syn::{parse_macro_input, Ident, LitInt};
use rug::{Float, float::Constant};

fn gen_sum(s: &Ident, idents: &[impl Borrow<Ident>]) -> TokenStream {
    let idents: Vec<_> = idents.iter().map(|x| x.borrow()).collect();
    quote! { let #s = #(#idents)*+; }
}

/// generates `(x, y) = two_sum(a, b)`
fn gen_two_add(x: &Ident, y: &Ident, a: &Ident, b: &Ident) -> TokenStream {
    quote! { let (#x, #y) = two_add(#a, #b); }
}

fn gen_two_sub(x: &Ident, y: &Ident, a: &Ident, b: &Ident) -> TokenStream {
    quote! { let (#x, #y) = two_sub(#a, #b); }
}

fn gen_two_mul(x: &Ident, y: &Ident, a: &Ident, b: &Ident) -> TokenStream {
    quote! { let (#x, #y) = two_mul(#a, #b); }
}

/// generates `(x, y) = quick_two_sum(a, b)`
fn gen_two_add_fast(x: &Ident, y: &Ident, a: &Ident, b: &Ident) -> TokenStream {
    quote! { let (#x, #y) = two_add_fast(#a, #b); }
}

fn vec_tokens(n: usize, prefix: &'static str) -> Vec<Ident> {
    (0..n).map(|x| format_ident!("{}{}", prefix, x)).collect()
}

fn gen_output(s: &[impl Borrow<Ident>], n: usize, sloppy: bool) -> TokenStream {
    let s: Vec<_> = s.iter().map(|x| x.borrow()).collect();
    if sloppy {
        assert!(s.len() == n);
        quote! { (#(#s),*) }
    } else {
        assert!(s.len() == n + 1);
        let prev = &s[0..n - 1];
        let si = &s[n - 1];
        let sn = &s[n];
        quote! { (#(#prev),*, #si + #sn) }
    }
}

fn gen_two_pass_renorm(n: usize, sloppy: bool) -> TokenStream {
    let a = vec_tokens(n + !sloppy as usize, "a");
    if n == 2 && sloppy {
        return quote! {
            #[inline]
            fn renormalize(a0: f64, a1: f64) -> Self {
                Self::from_tuple(two_add_fast(a0, a1))
            }
        };
    }

    let mut pass1 = vec![];
    for i in (0..n - sloppy as usize).rev() {
        pass1.push(gen_two_add_fast(&a[i], &a[i + 1], &a[i], &a[i + 1]));
    }

    let mut pass2 = vec![];
    for i in 1..n - 1 {
        pass2.push(gen_two_add_fast(&a[i], &a[i + 1], &a[i], &a[i + 1]));
    }

    let output = gen_output(&a, n, sloppy);
    quote! {
        #[inline]
        fn renormalize(#(#a: f64),*) -> Self {
            #(#pass1)*
            #(#pass2)*
            Self::new #output
        }
    }
}

fn gen_one_pass_renorm(n: usize, sloppy: bool) -> TokenStream {
    let a = vec_tokens(n + !sloppy as usize, "a");
    let s = vec_tokens(n + 1, "s");

    if n == 2 && sloppy {
        return quote! {
            #[inline]
            fn renormalize(a0: f64, a1: f64) -> Self {
                Self::from_tuple(quick_two_sum(a0, a1))
            }
        };
    }
    let mut pass = vec![];
    for i in 0..n {
        pass.push(gen_two_add_fast(&s[i], &s[i + 1], &a[i], &a[i + 1]));
    }
    let output = gen_output(&s, n, sloppy);

    quote! {
        #[inline]
        fn renormalize(#(#a: f64),*) -> Self {
            #(#pass)*
            Self::new #output
        }
    }
}

fn gen_merge_sum(
    inputs: &[(impl Borrow<Ident>, usize)],
    n_outputs: usize,
) -> (TokenStream, Vec<&Ident>) {
    assert!(inputs.len() >= n_outputs);

    let mut code = vec![];

    // a build has many *levels*
    let mut building = vec![VecDeque::new(); n_outputs];

    // put each ident to its corresponding level
    for &(ref x, l) in inputs {
        if l < n_outputs {
            building[l].push_back(x.borrow());
        }
    }

    // Goal: each level has at most 1 ident.
    for l in 0..n_outputs {
        while building[l].len() > 1 {
            let a = building[l].pop_front().unwrap();
            let b = building[l].pop_front().unwrap();

            if l == n_outputs - 1 {
                // last floor, we don't care about error term
                code.push(quote! { let #b = #a + #b; });
                building[l].push_back(b);
            } else {
                // the error term is then throwed to the next level
                code.push(gen_two_add(b, a, b, a));
                building[l].push_back(b);
                building[l + 1].push_back(a);
            }
        }
        assert!(building[l].len() == 1);
    }

    let outputs: Vec<_> = building.iter().map(|x| x[0]).collect();
    let code = quote! {
        {
            #(#code)*
            (#(#outputs,)*)
        }
    };
    (code, outputs)
}

fn destruct(src: &Ident, n: usize, prefix: &'static str) -> (TokenStream, Vec<Ident>) {
    let a = vec_tokens(n, prefix);
    let code = quote! { let (#(#a),*) = #src.to_tuple(); };
    (code, a)
}

fn gen_constructor(n: usize) -> TokenStream {
    let inputs: Vec<_> = (0..n).map(|x| format_ident!("a{}", x)).collect();
    let types = vec![format_ident!("f64"); n];
    let destruction: Vec<_> = (0..n).map(|x| quote! { self.data[#x] }).collect();

    quote! {
        #[inline]
        fn new(#(#inputs: f64),*) -> Self {
            Self { data: [#(#inputs),*] }
        }

        #[inline]
        fn from_tuple((#(#inputs),*): (#(#types),*)) -> Self {
            Self { data: [#(#inputs),*] }
        }

        #[inline]
        fn to_tuple(&self) -> (#(#types),*) {
            (#(#destruction),*)
        }
    }
}

fn gen_add_f64(multifloats: &Ident, n: usize, sloppy: bool) -> TokenStream {
    let b = format_ident!("b");
    let (destruct_a, a) = destruct(&format_ident!("self"), n, "a");

    let mut inputs: Vec<_> = a.iter().zip(0..n).collect();
    inputs.push((&b, 0));
    let (merge_sum, outputs) = gen_merge_sum(&inputs, n + !sloppy as usize);

    quote! {
        impl Add<f64> for #multifloats<#n> {
            type Output = Self;
            #[inline]
            fn add(self, #b: f64) -> Self {
                #destruct_a;
                let (#(#outputs,)*) = #merge_sum;
                Self::renormalize(#(#outputs,)*)
            }
        }
    }
}

fn gen_add(multifloats: &Ident, n: usize, sloppy: bool) -> TokenStream {
    let other = format_ident!("b");
    let (destruct_a, a) = destruct(&format_ident!("self"), n, "a");
    let (destruct_b, b) = destruct(&other, n, "b");

    let mut idents = vec![];
    let mut code = vec![];
    for i in 0..n {
        code.push(gen_two_add(&a[i], &b[i], &a[i], &b[i]));
        idents.push((&a[i], i));
        idents.push((&b[i], i + 1));
    }

    let (merge_sum, outputs) = gen_merge_sum(&idents, n + !sloppy as usize);

    quote! {
        impl Add for #multifloats<#n> {
            type Output = Self;
            #[inline]
            fn add(self, #other: Self) -> Self {
                #destruct_a;
                #destruct_b;
                #(#code)*
                let (#(#outputs,)*) = #merge_sum;
                Self::renormalize(#(#outputs,)*)
            }
        }
    }
}

fn gen_sub_f64(multifloats: &Ident, n: usize, sloppy: bool) -> TokenStream {
    let b = format_ident!("b");
    let (destruct_a, a) = destruct(&format_ident!("self"), n, "a");

    let mut code = vec![];
    code.push(gen_two_sub(&a[0], &b, &a[0], &b));
    let mut inputs: Vec<_> = a.iter().zip(0..n).collect();
    inputs.push((&b, 1));
    let (merge_sum, outputs) = gen_merge_sum(&inputs, n + !sloppy as usize);

    quote! {
        impl Sub<f64> for #multifloats<#n> {
            type Output = Self;
            #[inline]
            fn sub(self, #b: f64) -> Self {
                #destruct_a;
                #(#code)*
                let (#(#outputs,)*) = #merge_sum;
                Self::renormalize(#(#outputs,)*)
            }
        }
    }
}

fn gen_sub(multifloats: &Ident, n: usize, sloppy: bool) -> TokenStream {
    let other = format_ident!("b");
    let (destruct_a, a) = destruct(&format_ident!("self"), n, "a");
    let (destruct_b, b) = destruct(&other, n, "b");

    let mut idents = vec![];
    let mut code = vec![];
    for i in 0..n {
        code.push(gen_two_sub(&a[i], &b[i], &a[i], &b[i]));
        idents.push((&a[i], i));
        idents.push((&b[i], i + 1));
    }

    let (merge_sum, outputs) = gen_merge_sum(&idents, n + !sloppy as usize);

    quote! {
        impl Sub for #multifloats<#n> {
            type Output = Self;
            #[inline]
            fn sub(self, #other: Self) -> Self {
                #destruct_a;
                #destruct_b;
                #(#code)*
                let (#(#outputs,)*) = #merge_sum;
                Self::renormalize(#(#outputs,)*)
            }
        }
    }
}

fn gen_mul_f64(multifloats: &Ident, n: usize, sloppy: bool) -> TokenStream {
    let other = format_ident!("rhs");
    let (destruct_a, a) = destruct(&format_ident!("self"), n, "a");
    let b = vec_tokens(n, "b");

    let mut idents = vec![];
    let mut code = vec![];
    for i in 0..n {
        let ai = &a[i];
        let bi = &b[i];
        if i == n - 1 && sloppy {
            code.push(quote! { let #ai = #ai * #other; });
        } else {
            code.push(gen_two_mul(ai, bi, ai, &other));
            idents.push((ai, i));
            idents.push((bi, i + 1));
        }
    }

    let (merge_sum, outputs) = gen_merge_sum(&idents, n + !sloppy as usize);

    quote! {
        impl Mul<f64> for #multifloats<#n> {
            type Output = Self;
            #[inline]
            fn mul(self, #other: f64) -> Self {
                #destruct_a;
                #(#code);*
                let (#(#outputs,)*) = #merge_sum;
                Self::renormalize(#(#outputs,)*)
            }
        }
    }
}

fn gen_mul(multifloats: &Ident, n: usize, sloppy: bool) -> TokenStream {
    let other = format_ident!("rhs");
    let (destruct_a, a) = destruct(&format_ident!("self"), n, "a");
    let (destruct_b, b) = destruct(&other, n, "b");

    let mut idents = vec![];
    let mut code = vec![];
    for s in 0..=n - sloppy as usize {
        for i in 0..=s {
            let j = s - i;
            if i >= n || j >= n {
                continue;
            }
            let ai = &a[i];
            let bj = &b[j];
            if s == n - sloppy as usize {
                let x = format_ident!("s_{}", i);
                code.push(quote! { let #x = #ai * #bj; });
                idents.push((x, s));
            } else {
                let p = format_ident!("p_{}_{}", i, j);
                let q = format_ident!("q_{}_{}", i, j);
                code.push(gen_two_mul(&p, &q, ai, bj));
                idents.push((p, s));
                idents.push((q, s + 1));
            }
        }
    }
    let (merge_sum, outputs) = gen_merge_sum(&idents, n + !sloppy as usize);

    quote! {
        impl Mul for #multifloats<#n> {
            type Output = Self;
            #[inline]
            fn mul(self, #other: Self) -> Self {
                #destruct_a;
                #destruct_b;
                #(#code);*
                let (#(#outputs,)*) = #merge_sum;
                Self::renormalize(#(#outputs,)*)
            }
        }
    }
}

fn gen_div_f64(multifloats: &Ident, n: usize, sloppy: bool) -> TokenStream {
    let other = format_ident!("rhs");
    let m = n + !sloppy as usize;
    let q = vec_tokens(m, "q");
    let r = format_ident!("r");

    let mut code = vec![];
    code.push(quote! { let #r = self; });

    for i in 0..m {
        let qi = &q[i];
        code.push(quote! { let #qi = self.data[0] / #other; });
        if i < m - 1 {
            code.push(quote! { let #r = #r - #other * #qi; }); // TODO: #other * #qi is floating point operation!
        }
    }
    quote! {
        impl Div<f64> for #multifloats<#n> {
            type Output = Self;
            fn div(self, #other: f64) -> Self {
                self / Self::mp(#other)
                // #(#code)*;
                // Self::renormalize(#(#q),*)
            }
        }
    }
}

fn gen_div(multifloats: &Ident, n: usize, sloppy: bool) -> TokenStream {
    let other = format_ident!("rhs");
    let m = n + !sloppy as usize;
    let q = vec_tokens(m, "q");
    let r = format_ident!("r");

    let mut code = vec![];
    code.push(quote! { let #r = self; });

    for i in 0..m {
        let qi = &q[i];
        code.push(quote! { let #qi = r.data[0] / #other.data[0]; });
        if i < m - 1 {
            code.push(quote! { let #r = #r - #other * #qi; });
        }
    }
    quote! {
        impl Div for #multifloats<#n> {
            type Output = Self;
            fn div(self, #other: Self) -> Self {
                #(#code)*;
                Self::renormalize(#(#q),*)
            }
        }
    }
}

fn gen_float_consts(multifloats: &Ident, n: usize) -> TokenStream {
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

    quote! {
        impl #multifloats<#n> {
            #(#def)*
        }

        impl num::traits::FloatConst for #multifloats<#n> {
            #(#implement)*
        }
    }
}

#[proc_macro]
pub fn impl_f64x(tokens: proc_macro::TokenStream) -> proc_macro::TokenStream {
    let input = parse_macro_input!(tokens as LitInt);
    let n = input.base10_parse::<usize>().unwrap();
    let multifloats = format_ident!("f64x");
    let f64xn = format_ident!("f64x{}", n);
    let sloppy = false;

    let renormalize = gen_two_pass_renorm(n, sloppy);
    let constructor = gen_constructor(n);
    let add = gen_add(&multifloats, n, sloppy);
    let sub = gen_sub(&multifloats, n, sloppy);
    let mul = gen_mul(&multifloats, n, sloppy);
    let div = gen_div(&multifloats, n, sloppy);
    let add_f64 = gen_add_f64(&multifloats, n, sloppy);
    let sub_f64 = gen_sub_f64(&multifloats, n, sloppy);
    let mul_f64 = gen_mul_f64(&multifloats, n, sloppy);
    let div_f64 = gen_div_f64(&multifloats, n, sloppy);
    let float_consts = gen_float_consts(&multifloats, n);
    let expanded = quote! {
        impl #multifloats<#n> {
            #constructor
            #renormalize
        }
        #add
        #sub
        #mul
        #div
        #add_f64
        #sub_f64
        #mul_f64
        #div_f64
        #float_consts

        type #f64xn = #multifloats<#n>;
    };
    println!("{}", expanded);
    expanded.into()
}
