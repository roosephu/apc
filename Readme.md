# $\newcommand{\d}{\mathrm{d}}$About

This is a side project. I'm new to analytic number theory and can't understand most of the equations in the reference papers. I know little computer architecture so the program is also highly-unoptimized. Use at your own risk. 

The program hasn't been tested so it might produce incorrect numbers or even fail to produce reasonable results. 

### Is it fast? How fast is it? 

No, it's very slow, compared to [primecount](https://github.com/kimwalisch/primecount). It's even slow compared without LMFDB. 

With that being said, it's not that slow as people might think. For example, with the help of LMFDB, it successfully calculates $\pi(10^{15}) = 29844570422669$ in about 20s in my laptop (Macbook Pro 13-inch, 2017) with a single thread. As a comparison, `primecount` gives the same answer in about 1.7s in the same environment (`primecount 1e15 -t 1`), so it's 12x faster than my code. 

The algorithm can run in parallel easily but I don't do that. 

### Is using LMFDB cheating? 

Yes, kind of. The algorithm needs all non-trivial zeros of $\zeta(\frac{1}{2} + it)$ for $t \leq \tilde O(x^{1/2})$ to compute $\pi(x)$ in $\tilde O(x^{1/2})$ time. (I always assumes RH holds. ) However, a table of $\{\pi(k^2)\}_{k^2 \leq x}$ also enables computing $\pi(x)$ in $\tilde O(x^{1/2})$ time, with a more trivial algorithm. The only point is that the table of $\{\pi(k^2)\}_{k^2 \leq x}$ might be very hard to get :( I don't know any algorithm which can output such a table in $N^{o(1)}$ time, while it seems that the non-trivial zeros of $\zeta(\frac{1}{2} + it)$ for $t \leq T$ can be found in $\tilde O(T)$ time. 

### If a number is produced, is it guaranteed to be correct? 

No. There are some possible reasons:

1. I added my own heuristics and traded correctness for speed. For example, I use heuristics to decide the interval to compute $\Delta$ (see technical details below). 
2. Unlike [Platt], I don't use interval/ball arithmetic (MPFI/Arb) and don't do error analysis. I don't even use existing high-precision arithmetic library, so tracking the rounding error alone is nearly impossible. 
3. Also for the logic code: it's not tested, so there might be bugs. 

### Why Rust?

Because it's new to me either. What can the result of writing an analytic number theory algorithm by Rust, both of which I'm not familiar with, be? 

Unfortunately, after writing these code, I still don't think I've learnt Rust. Some of the features confuse me. For example, can anyone pointing me out the correct way to enable `a + b` where `a` is a `T` and `b` is `Complex<T>` for a custom float type `T`? We can't `impl<T> Add<Complex<T>> for T`  due to the orphan rule, and `num::Complex` does this by using macros, which are not exported. 

# Some technical details

## Multi-precision data type

I've implemented my own double-double data type. Although there are many existing libraries, like [TwoFloat](https://github.com/ajtribick/twofloat), I'd like to invent a new wheel for my own needs. With that being said, I appreciate their work and my implementations are heavily inspired by them. Kudos to them. 

## Heuristic to decide the interval to calculate $\Delta$

People might notice that my sieving interval seems to be shorter than others. That's true because I applied my heuristics, which is based on Cramér's random model. 

For brevity, we define $f(u) = \phi(u) - [u \leq x]$, i.e., $f(u) = \phi(u)$ for $u > x$ and $f(u) = \phi(u) - 1$ for $u < x$. 

Given a real $d > 0$,  define $U_L = [1, xe^{-d\lambda}], U_R = [xe^{d \lambda}, \infty)$ and $U = U_L \cup U_R$. We'd like to sieve out all primes in $(xe^{-d\lambda}, xe^{d \lambda})$ to compute $\Delta$ and estimate the truncation error $S = \sum_{p \in U}f(p)$. Please keep in mind that $\lambda = O(x^{-1/2})$ and $d = O(1)$, so $x e^{d\lambda} = x - O(x^{1/2})$ and $x e^{d\lambda} = x + O(x^{1/2})$. 

By Cramér's random model, we assume
$$
S = \sum_{p \in U} f(p) \approx  \sum_{u \in U} \frac{f(u)}{\ln u},
$$
One way to improve the approximation is to add an error bar, which is done by studying the behavior of a stochastic estimation $\hat S$: 
$$
\hat S = \sum_{u \in U} X_u f(u), \quad X_u \sim \text{Bernoulli}(1/\ln u),
$$
We use Bernstein's inequality to approximate $S = \mathbb{E} [\hat S] + t$ by solving $t$ for $\epsilon = \exp(-20)$ in
$$
\exp\left( -\frac{\frac{1}{2} t^2}{\text{Var}[\hat S] + \frac{1}{3} M t} \right) \leq \epsilon, \quad M := \max_{u \in U} |f(u)| = \Phi(d).
$$
Let $q = 1/\ln x$, so $q \approx 1/\ln u$ for reasonable $u \in U$, i.e., those with not-too-small $f(u)$. 

Next, we approximate $\mathbb{E}[\hat S]$:
$$
\mathbb{E}[\hat S] = \sum_{u \in U} \frac{f(u)}{\ln u} \approx \sum_{u \in U} q f(u) \approx q (I_r - I_l),
$$
where
$$
I_r = x \lambda \int_{d}^{\infty} \Phi(t) e^{t\lambda} \d t, \quad I_l = x \lambda \int_{d}^\infty \Phi(t) e^{-t\lambda} \d t.
$$
The last approximation is similar to [Section 3.2, Galway] and both $I_l$ and $I_r$ have a closed form solution. Different from [Galway], we are making use of the fact that $f(u) < 0$ for $u < x$ and $f(u) > 0$ for $u > x$ and the sum of them might cancel. 

The variance of $\hat S$ can be approximated by:
$$
\text{Var}[\hat S] = \sum_{u \in U} f(u)^2 \frac{1}{\ln u} (1 - 1/\ln u) \approx q(1-q) \sum_{u \in U} f(u)^2 \leq q(1-q) M \sum_{u \in U} |f(u)| \approx M q(1-q) (I_l + I_r).
$$
Finally, we minimize $d$ so that $S \leq 0.24$. In this way, we can get a shorter interval to sieve at a loss of guarantee.

In practice, when $x = 10^{11}$ and $\lambda = 3 \times 10^{-5}$, [Galway] predicts an interval of length 3e7 while this predicts an interval of length 2e7, while the difference between two sums is around 0.0112. 

## Precision loss when computing $\Delta$ using `f64`

Recall that
$$
\Delta = \sum_{L \leq p \leq R} f(p),
$$
so there are two possible precision loss: 

## Deciding Integral upper bound

​    */// we are integrating N(h) \hat_\phi(s), which is approximately x^σ exp(-λ^2 h^2 / 2) with σ = 0.5 or 1.*

​    */// |N(h) \hat_\phi(0.5 + ih)| \approx \log(h) x^0.5 exp(-λ^2 h^2 / 2)*

​    */// |\hat_\phi(1 + ih)| ≈ x^0.5 \exp(-λ^2 h^2 / 2)*

## Integrating $\hat\phi(s)$

Here we'd like to integrate $\hat\phi(s) = \frac{x^s}{s} \exp(\lambda^2 s^2 / 2)$. We basically follow the method in [Section 6, Platt] but with a small modification. 

The high level idea is that: We express $\hat\phi(s_0 + ih) = \hat\phi(s_0) f(z) \exp(z)$ for some complex number $z = wh$ and function $f(z)$, and find a polynomial $P(z)$ to approximate $f(z)$ locally, then we apply integration by parts repeatedly:
$$
\int P(z) \exp(z) \d z = \exp(z) \sum_{k=0}^\infty P^{(k)}(z) (-1)^k = \exp(z) Q(z), \quad Q(z):=\sum_{k=0}^\infty P^{(k)}(z)(-1)^k.
$$
[Platt] chooses the following form
$$
\begin{aligned}
\hat\phi(s_0 + ih) & = \frac{x^{s_0 + ih}}{s_0 + ih} \exp(\lambda^2 (s_0 + ih)^2 / 2) \\
 & = \hat\phi(s_0) \frac{x^{ih}}{1 + \frac{ih}{s0}} \exp(\lambda^2 (-h^2+2ihs_0)/2) \\
 & = \hat\phi(s_0) \frac{\exp(-\lambda^2 h^2/2)}{1 + \frac{i}{s_0}h} \exp(ih (\lambda^2 s_0 + \ln x)),
\end{aligned}
$$
while we choose 
$$
\hat\phi(s_0 + ih)  = \hat\phi(s_0) \frac{\exp(-\lambda^2 h^2/2 + ih \lambda^2 s_0)}{1 + \frac{i}{s_0}h} \exp(ih \ln x).
$$
The benefit of ours is that $ih\ln x$ has only imaginary part, so it's faster to compute the exponential function. This can speed it up a little bit (5% perhaps) when $x$ is not too large that LMFDB has enough zeros of $\zeta(s)$. 

