/*!
# Analytic Prime Counting

## Questions that nobody cares

### Is it fast? How fast is it?

No, it's slow, compared to [primecount](https://github.com/kimwalisch/primecount).

With that being said, it's not that slow as people might think. For example, with the help of LMFDB, it successfully calculates $\pi(10^{18}) = 24739954287740860$ in about 8.5min in Macbook Pro 13-inch (2017) with a single thread. As a comparison, `primecount` gives the same answer in about 100s in the same environment (`primecount 1e18 -t 1`), so the analytic method is ~5x slower.

### If a number is produced, is it guaranteed to be correct?

No. There are some possible reasons:

1. I added my own heuristics and traded correctness for speed. For example, I use heuristics to decide the interval to compute $\Delta$ (see technical details below).
2. Unlike [Platt], I don't use interval/ball arithmetic (MPFI/Arb) and don't do error analysis. I don't even use existing high-precision arithmetic library, so tracking the rounding error alone is nearly impossible.
3. Also for the logic code: it's not tested, so there might be bugs.

### Why Rust?

Because it's new to me either. What can the result of writing an analytic number theory algorithm by Rust, both of which I'm not familiar with, be?

Unfortunately, after writing these code, I still don't think I've learnt Rust. Some of the features confuse me. For example, can anyone pointing me out the correct way to enable `a + b` where `a` is a `T` and `b` is `Complex<T>` for a custom float type `T`? We can't `impl<T> Add<Complex<T>> for T`  due to the orphan rule, and `num::Complex` does this by using macros, which are not exported.

### Why [Platt] instead of [Büthe] or [FKBJ]?

Because they seems too difficult to me. Welcomed by the fancy Greek letters and subscripts, I clicked X immediately.

[Platt] seems to be quite straightforward, especially given that I've read [Galway] and understood the high-level idea. I appreciate a lot for how easy to follow [Galway] is. Without it, I wouldn't ever have the idea to implement such an algorithm, probably. (Perhaps because it's a thesis instead of a paper?)


## LMFDB

Using LMFDB is kind of cheating. The algorithm needs all non-trivial zeros of $\zeta(\frac{1}{2} + it)$ for $t \leq \softO(x^{1/2})$ to compute $\pi(x)$ in $\tilde O(x^{1/2})$ time. (To make life easier, I assume RH holds. ) However, a table of $\pi(k^2)$ for $k^2 \leq x$ also enables computing $\pi(x)$ in $\softO(x^{1/2})$ time, with a much more trivial algorithm. The only point is that the table of $\pi(k^2)$ might be very hard to get :( I don't know any algorithm which can output such a table in $N^{o(1)}$ time, while it seems that the non-trivial zeros of $\zeta(\frac{1}{2} + it)$ for $t \leq T$ can be found in $\softO(T)$ time.

Finding all non-trivial zeros of $\zeta$ seems to be another line of work. I tried to read [Platt2] but unfortunately didn't understand it. [Gourdon] also seems promising and I have some preliminary code for that.

## Float data type

I've implemented my own double-double data type. Although there are many existing libraries, like [TwoFloat](https://github.com/ajtribick/twofloat), I'd like to invent a new wheel for my own needs. With that being said, I appreciate their work and my implementations are heavily inspired by them. Kudos to them.

As mentioned [here](http://www.math.uni-bonn.de/people/jbuethe/topics/AnalyticPiX.html),

>Our implementation focuses rather on speed, using a 128 Bit fixed-point data type developed by J. Franke.

so I think double-double suffices for this project. I'm surprised to see that Jens Franke, who's Jan Büthe's advisor, developed the data type. Perhaps it's from his other projects.


## References

+ [[Platt](https://arxiv.org/abs/1203.5712)] : David Platt, “Computing $\pi(x)$ Analytically.”
+ [[Platt2](https://www.ams.org/journals/mcom/2017-86-307/S0025-5718-2017-03198-7/S0025-5718-2017-03198-7.pdf)] : David Platt, “Isolating some non-trivial zeros of zeta”
+ [[FKBJ](https://www.ams.org/journals/mcom/2017-86-308/S0025-5718-2017-03038-6/S0025-5718-2017-03038-6.pdf)] : Franke et al., “A Practical Analytic Method for Calculating $\pi (x)$.”
+ [[Büthe](https://arxiv.org/abs/1410.7008)] : Jan Büthe, “An Improved Analytic Method for Calculating $\pi(x)$.”
+ [[Büthe2](https://arxiv.org/abs/1410.7561)] : Jan Büthe, "A Brun-Titchmarsh inequality for weighted sums over prime numbers"
+ [[Galway](https://faculty.math.illinois.edu/~galway/PhD_Thesis/thesis-twoside.pdf)] : William Galway, "Analytic Computation of the Prime-Counting Function".
+ [[Gourdon](http://numbers.computation.free.fr/Constants/Miscellaneous/zetazeros1e13-1e24.pdf)] : Xavier Gourdon, "The $10^{13}$ first zeros of the Riemann Zeta function, and zeros computation at very large height"

I also briefly summarized these papers [here](https://www.roosephu.xyz/2021/03/08/analytic-pi/) and [here](https://www.roosephu.xyz/2021/03/08/analytic-pi-2/) (in Chinese).

*/

#![allow(non_snake_case)]
#![allow(dead_code)]
#![allow(confusable_idents)]
#![allow(uncommon_codepoints)]
#![allow(clippy::float_cmp)]
#![allow(clippy::needless_range_loop)]
#![allow(clippy::let_and_return)]
#![allow(clippy::approx_constant)]
#![feature(trait_alias)]
#![feature(destructuring_assignment)]
#![feature(type_alias_impl_trait)]

cfg_if::cfg_if! {
    if #[cfg(feature = "zeta")] {
        mod sum_trunc_dirichlet;
        mod bandwidth_interp;
        mod gamma;
        // mod riemann_siegel;
        // mod zeta;
        // pub use riemann_siegel::{RiemannSiegelZ, RiemannSiegelZeta};
        // pub use zeta::ZetaGalway;
    }
}

cfg_if::cfg_if! {
    if #[cfg(feature = "galway")] {
        mod galway;
        mod context;
        pub use galway::{Galway, GalwayHints};
        pub use context::Context;
    }
}

mod brentq;
mod cache_stat;
mod constants;
mod fast_phi;
mod lmfdb;
mod platt;
mod platt_integral;
mod power_series;
mod sieve;
mod traits;

pub use fast_phi::LittlePhiSum;
pub use platt::PlattBuilder;

cfg_if::cfg_if! {
    if #[cfg(feature = "doc")] {
        pub use crate::platt::Platt;
        pub use crate::platt_integral::{ExpansionIntegrator, HybridPrecIntegrator};
    }
}
