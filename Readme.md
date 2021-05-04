$\newcommand{\d}{\mathrm{d}}$

# Integral in [Platt]

Here we'd like to integral $\hat\phi(s) = \frac{x^s}{s} \exp(\lambda^2 s^2 / 2)$. 

We follow the method in [Section 6, Platt]. Let $s = \sigma + it$. 

## Method 1

We expand $\hat\phi(s)$ around $s_0$:
$$
\begin{aligned}
\hat\phi(s_0 + ih) & = \frac{x^{s_0 + ih}}{s_0 + ih} \exp(\lambda^2 (s_0 + ih)^2 / 2) \\
 & = \hat\phi(s_0) \frac{x^{ih}}{1 + \frac{ih}{s0}} \exp(\lambda^2 (-h^2+2ihs_0)/2) \\
 & = \hat\phi(s_0) \frac{\exp(-\lambda^2 h^2/2)}{1 + \frac{i}{s_0}h} \exp(ih (\lambda^2 s_0 + \ln x)),
\end{aligned}
$$
For brevity, we define $c = i(\lambda^2 s_0 + \ln x)$ and $a = -\lambda^2 / 2c^2$, $b = -\frac{i}{s_0 c}$, so it becomes: for $z = ch$,
$$
\hat\phi(s_0 + ih) = \hat\phi(s_0) \frac{\exp(a (ch)^2)}{1 - b(ch)} \exp(ch) = \hat\phi(s_0) \frac{\exp(az^2)}{1 - bz} \exp(z).
$$
Next, we find a polynomial $P(z)$ such that $P(z) \approx \frac{\exp(cz^2)}{1 - bz}$ for $|z| \leq \delta$. Simple truncated Taylor series should suffice. Now we are focusing on the integral:
$$
\int \hat\phi(s_0 + ih) \d h = \hat\phi(s_0) \int \frac{\exp(a(ch)^2)}{1 - b(ch)} \exp(ch) \d h \approx \frac{\hat\phi(s_0)}{c} \int P(z) \exp(z) \d z.
$$
By integration by parts, we have
$$
\int P(z) \exp(z) \d z = \exp(z) \sum_{k=0}^{\infty} P^{(k)}(z) (-1)^{k} = \exp(z) Q(z).
$$
Now what we do is:

1. compute $P(z)$ and $\delta$ such that $\left|P(z) - \frac{\exp(cz^2)}{1 - bz} \right| \leq \varepsilon$ for $|z| \leq \delta$. 
2. compute $Q(z) = \sum\limits_{k=0}^{\infty} P^{(k)}(z) (-1)^k$.

## Further optimization

The drawback in previous method is that $\exp(z)$ for complex is slow... Now we'd like to optimize it so that only $\exp(ih)$ for a real $h$ is computed. The optimization is simple: also expand $\exp(ih \lambda^2 s_0)$. 
$$
\hat\phi(s_0 + ih)  = \hat\phi \frac{\exp(-\lambda^2 h^2/2 + ih \lambda^2 s_0)}{1 + \frac{i}{s_0}h} \exp(ih \ln x).
$$
Note that $|\lambda^2 s_0|$ is typically very small as $\lambda = O(x^{-1/2})$ and $|s_0| = O(x^{1/2})$.

We simply set $c = i \ln x$, $a = -\lambda^2 / 2c^2$, $b = -\frac{i}{s_0 c}$ and $d = i h \lambda^2 s_0/c$, so 
$$
\hat\phi(s) = \hat\phi(s_0) \frac{\exp(a z^2 + dz)}{1 - bz} \exp(z),
$$
and everything else can continue.

