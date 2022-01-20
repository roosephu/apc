# apc

This program counts the number of primes under $x$ using analytic method for $x \leq 10^{18}$ based on the paper [Platt].

This is a side project. I'm new to analytic number theory and can't understand most of the equations in the reference papers. I know little computer architecture so the program is also highly-unoptimized. The program hasn't been tested so it might produce incorrect numbers or even fail to produce reasonable results. Use at your own risk.

[notes](https://apc.roosephu.xyz/apc)

# Build

1. Install Rust nightly, by first installing [Rustup](https://www.rust-lang.org/tools/install) and then installing Rust nightly.

   ```bash
   rustup toolchain install nightly
   ```

2. Install [primesieve](https://github.com/kimwalisch/primesieve) which `apc` uses to sieve prime numbers.

3. Build: `cargo +nightly build --bin apc --release`.

To build the doc:

```bash
RUSTDOCFLAGS="--html-in-header `pwd`/src/latex.html" cargo doc --features doc --no-deps
```

# Usage

You may first try

```bash
target/release/apc pi "10^13" --lambda-hint 5000
```

which takes about 10s. There are two possible improvement to speed it up, by providing it a table of nontrivial roots of Riemann $\zeta$ function (see next subsection).

You may try different `--lambda-hint` value. Larger value sieves more and use fewer $\zeta$ zeros, vice versa. `apc` also suggests `--lambda-hint` if you export `RUST_LOG=info`.

Here are a few other settings that I've tried (if you have LMFDB under `./data/lmfdb`):

```bash
target/release/apc pi "10^15" --lambda-hint 10 --lmfdb ./data/lmfdb  # 10s
target/release/apc pi "10^18" --lambda-hint 50 --lmfdb ./data/lmfdb  # 8.5min
```

The time is measured in a Macbook Pro (either in 2017 version or M1 Pro version).

## Zeros of Riemann $\zeta$ function

As the algorithm uses the nontrivial roots of Riemann $\zeta$ function, which are expensive to compute. `apc` provides three methods to do so:
1. Computing them online. This is the default behavior.
2. Reading from [LMFDB](https://beta.lmfdb.org/data/riemann-zeta-zeros/) by `--lmfdb path/to/LMFDB`.
3. Reading from cached result by `--apcdb path/to/APCDB`.

In the last two methods, please make sure that the file `md5.txt` exists as `apc` uses it to index the database. However, you don't need the whole database. To generate the database by `apc`, use the following commands:

```bash
mkdir -p ./data/apcdb
target/release/apc zero "10^15" --apcdb ./data/apcdb
cd ./data/apcdb
md5sum *.dat > md5.txt
cd ../..
```

Please be aware that computing the table of nontrivial roots of Riemann $\zeta$ function can be quite slow (12k zeros per second).

