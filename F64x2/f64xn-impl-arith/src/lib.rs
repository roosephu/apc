#![allow(dead_code)]

// mod certified_ops;
// mod eft;

use std::{borrow::Borrow, collections::VecDeque};

use proc_macro2::TokenStream;
use quote::{format_ident, quote};
use syn::{parse_macro_input, Ident, LitInt};

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
        quote! { [#(#s),*] }
    } else {
        assert!(s.len() == n + 1);
        let prev = &s[0..n - 1];
        let si = &s[n - 1];
        let sn = &s[n];
        quote! { [#(#prev),*, #si + #sn] }
    }
}

fn gen_two_pass_renorm(n: usize, sloppy: bool) -> TokenStream {
    let a = vec_tokens(n + !sloppy as usize, "a");
    if n == 2 && sloppy {
        return quote! {
            #[inline]
            fn renormalize([a0, a1]: [f64; 2]) -> [f64; 2] {
                let (a0, a1) = two_add_fast(a0, a1);
                [a0, a1]
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
    let n_inputs = n + !sloppy as usize;
    quote! {
        #[inline]
        fn renormalize([#(#a),*]: [f64; #n_inputs]) -> [f64; #n] {
            #(#pass1)*
            #(#pass2)*
            #output
        }
    }
}

fn gen_one_pass_renorm(n: usize, sloppy: bool) -> TokenStream {
    let a = vec_tokens(n + !sloppy as usize, "a");
    let s = vec_tokens(n + 1, "s");

    if n == 2 && sloppy {
        return quote! {
            #[inline]
            fn renormalize([a0, a1]: [f64; 2]) -> [f64; 2] {
                let (a0, a1) = quick_two_sum(a0, a1);
                [a0, a1]
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
        fn renormalize([#(#a),*]: [f64; #n]) -> [f64; #n] {
            #(#pass)*
            #output
        }
    }
}

fn gen_merge_sum(inputs: &[(impl Borrow<Ident>, usize)], n_outputs: usize) -> TokenStream {
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
            [#(#outputs,)*]
        }
    };
    code
}

fn destruct(src: &Ident, n: usize, prefix: &'static str) -> (TokenStream, Vec<Ident>) {
    let a = vec_tokens(n, prefix);
    // let code = quote! { let (#(#a),*) = #src.to_tuple(); };
    let code = quote! { let [#(#a),*] = #src.data; };
    (code, a)
}

fn gen_constructor(n: usize) -> TokenStream {
    let inputs: Vec<_> = (0..n).map(|x| format_ident!("a{}", x)).collect();

    quote! {
        #[inline]
        fn new(#(#inputs: f64),*) -> Self {
            Self { data: [#(#inputs),*] }
        }
    }
}

fn gen_add_f64(n: usize, sloppy: bool) -> TokenStream {
    let a = vec_tokens(n, "a");
    let b = format_ident!("b");

    let mut inputs: Vec<_> = a.iter().zip(0..n).collect();
    inputs.push((&b, 0));
    let merge_sum = gen_merge_sum(&inputs, n + !sloppy as usize);

    quote! {
        #[inline]
        fn add_f64([#(#a,)*]: [f64; #n], #b: f64) -> [f64; #n] {
            let output = #merge_sum;
            Self::renormalize(output)
        }
    }
}

fn gen_add(n: usize, sloppy: bool) -> TokenStream {
    let a = vec_tokens(n, "a");
    let b = vec_tokens(n, "b");

    let mut idents = vec![];
    let mut code = vec![];
    for i in 0..n {
        code.push(gen_two_add(&a[i], &b[i], &a[i], &b[i]));
        idents.push((&a[i], i));
        idents.push((&b[i], i + 1));
    }

    let merge_sum = gen_merge_sum(&idents, n + !sloppy as usize);

    quote! {
        #[inline]
        fn add([#(#a,)*]: [f64; #n], [#(#b,)*]: [f64; #n]) -> [f64; #n] {
            #(#code)*
            let output = #merge_sum;
            Self::renormalize(output)
        }
    }
}

fn gen_sub_f64(n: usize, sloppy: bool) -> TokenStream {
    let a = vec_tokens(n, "a");
    let b = format_ident!("b");

    let mut code = vec![];
    code.push(gen_two_sub(&a[0], &b, &a[0], &b));
    let mut inputs: Vec<_> = a.iter().zip(0..n).collect();
    inputs.push((&b, 1));
    let merge_sum = gen_merge_sum(&inputs, n + !sloppy as usize);

    quote! {
        #[inline]
        fn sub_f64([#(#a,)*]: [f64; #n], #b: f64) -> [f64; #n] {
            #(#code)*
            let output = #merge_sum;
            Self::renormalize(output)
        }
    }
}

fn gen_sub(n: usize, sloppy: bool) -> TokenStream {
    let a = vec_tokens(n, "a");
    let b = vec_tokens(n, "b");

    let mut idents = vec![];
    let mut code = vec![];
    for i in 0..n {
        code.push(gen_two_sub(&a[i], &b[i], &a[i], &b[i]));
        idents.push((&a[i], i));
        idents.push((&b[i], i + 1));
    }

    let merge_sum = gen_merge_sum(&idents, n + !sloppy as usize);

    quote! {
        #[inline]
        fn sub([#(#a,)*]: [f64; #n], [#(#b,)*]: [f64; #n]) -> [f64; #n] {
            #(#code)*
            let output = #merge_sum;
            Self::renormalize(output)
        }
    }
}

fn gen_mul_f64(n: usize, sloppy: bool) -> TokenStream {
    let a = vec_tokens(n, "a");
    let b = format_ident!("b");
    let c = vec_tokens(n, "c");

    let mut idents = vec![];
    let mut code = vec![];
    for i in 0..n {
        let ai = &a[i];
        let ci = &c[i];
        if i == n - 1 && sloppy {
            code.push(quote! { let #ai = #ai * #b; });
        } else {
            code.push(gen_two_mul(ai, ci, ai, &b));
            idents.push((ai, i));
            idents.push((ci, i + 1));
        }
    }

    let merge_sum = gen_merge_sum(&idents, n + !sloppy as usize);

    quote! {
        #[inline]
        fn mul_f64([#(#a,)*]: [f64; #n], #b: f64) -> [f64; #n] {
            #(#code);*
            let output = #merge_sum;
            Self::renormalize(output)
        }
    }
}

fn gen_mul(n: usize, sloppy: bool) -> TokenStream {
    let a = vec_tokens(n, "a");
    let b = vec_tokens(n, "b");

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
    let merge_sum = gen_merge_sum(&idents, n + !sloppy as usize);

    quote! {
        #[inline]
        fn mul([#(#a,)*]: [f64; #n], [#(#b,)*]: [f64; #n]) -> [f64; #n] {
            #(#code);*
            let output = #merge_sum;
            Self::renormalize(output)
        }
    }
}

fn gen_div_f64(n: usize, _: bool) -> TokenStream {
    // let other = format_ident!("rhs");
    // let m = n + !sloppy as usize;
    // let q = vec_tokens(m, "q");
    // let r = format_ident!("r");

    // let mut code = vec![];
    // code.push(quote! { let #r = self; });

    // for i in 0..m {
    //     let qi = &q[i];
    //     code.push(quote! { let #qi = self.data[0] / #other; });
    //     if i < m - 1 {
    //         code.push(quote! { let #r = #r - #other * #qi; }); // TODO: #other * #qi is floating point operation!
    //     }
    // }

    let a = format_ident!("a");
    let b = format_ident!("b");
    quote! {
        #[inline]
        fn div_f64(#a: [f64; #n], #b: f64) -> [f64; #n] {
            Self::div(#a, Self::mp(#b).data)
            // #(#code)*;
            // Self::renormalize(#(#q),*)
        }
    }
}

fn gen_div(n: usize, sloppy: bool) -> TokenStream {
    let a = format_ident!("a");
    let b = format_ident!("b");
    let m = n + !sloppy as usize;
    let q = vec_tokens(m, "q");

    let mut code = vec![];

    for i in 0..m {
        let qi = &q[i];
        code.push(quote! { let #qi = a[0] / #b[0]; });
        if i < m - 1 {
            code.push(quote! { let #a = Self::sub(#a, Self::mul_f64(#b, #qi)); });
        }
    }
    quote! {
        fn div(#a: [f64; #n], #b: [f64; #n]) -> [f64; #n] {
            #(#code)*;
            Self::renormalize([#(#q,)*])
        }
    }
}

fn impl_f64xn_array_ops(n: usize, sloppy: bool) -> TokenStream {
    let renormalize = gen_two_pass_renorm(n, sloppy);
    let constructor = gen_constructor(n);
    let add = gen_add(n, sloppy);
    let sub = gen_sub(n, sloppy);
    let mul = gen_mul(n, sloppy);
    let div = gen_div(n, sloppy);
    let add_f64 = gen_add_f64(n, sloppy);
    let sub_f64 = gen_sub_f64(n, sloppy);
    let mul_f64 = gen_mul_f64(n, sloppy);
    let div_f64 = gen_div_f64(n, sloppy);
    quote! {
        impl f64x::<#n> {
            #constructor
            #renormalize
            #add
            #sub
            #mul
            #div
            #add_f64
            #sub_f64
            #mul_f64
            #div_f64
        }
    }
}

#[proc_macro]
pub fn f64xn_impl_arith(tokens: proc_macro::TokenStream) -> proc_macro::TokenStream {
    let input = parse_macro_input!(tokens as LitInt);
    let n = input.base10_parse::<usize>().unwrap();
    let sloppy = false;

    let expanded = impl_f64xn_array_ops(n, sloppy);
    // println!("{}", expanded);
    expanded.into()
}
