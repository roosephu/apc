# apc

This program counts the number of primes under $x$ using analytic method for $x \leq 10^{18}$ **assuming access to a table of non-trivial zeros of Riemann $\zeta(s)$ function**, based on the paper [Platt].

This is a side project. I'm new to analytic number theory and can't understand most of the equations in the reference papers. I know little computer architecture so the program is also highly-unoptimized. The program hasn't been tested so it might produce incorrect numbers or even fail to produce reasonable results. Use at your own risk.

[notes](https://apc.roosephu.xyz/apc)

# Build

1. Install Rust nightly, by first installing [Rustup](https://www.rust-lang.org/tools/install) and then installing Rust nightly.

   ```bash
   rustup toolchain install nightly
   ```

2. Download the files of zeros of $\zeta(s)$ from [LMFDB](https://beta.lmfdb.org/data/riemann-zeta-zeros/). To compute $\pi(10^{18})$, you might need all zeros up to $1.4 \times 10^8$.

   The directory should at least contain `md5.txt`. Also put data files (e.g.,`zeros_14.dat`) in this folder.

3. Build: `cargo +nightly build --bin apc --release`.

To build the doc:

```bash
RUSTDOCFLAGS="--html-in-header `pwd`/src/latex.html" cargo doc --features doc
```

# Usage

1. Run: `target/release/apc 10^15 --lambda-hint 10`. Time might vary, but probably it should finish under 1min. As a reference, it's ~8s in Macbook Pro 13-inch (2017).

   1. Add ` --zeta-zeros-path <your-LMFDB-path>` if you didn't put the LMFDB data under `data/zeros/`.

You may try different `--lambda-hint` value. Larger value sieves more and use fewer $\zeta(s)$ zeros, while smaller value does the opposite. The program also suggests `--lambda-hint` if you pass `--verbose` to the command line.

You may also want to try `target/release/apc 10^18 --lambda-hint 50`, which takes about 8.5min in Macbook Pro 13-inch (2017). Using more $\zeta$ zeros (and then reducing `--lambad-hint`) can reduce the time.

