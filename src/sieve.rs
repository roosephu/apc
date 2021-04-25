use log::info;

#[inline(never)]
pub(crate) fn linear_sieve(n: u64) -> Vec<u64> {
    info!("sieve primes up to {} using linear sieve", n);
    let mut mark = bit_vec::BitVec::from_elem(n as usize + 1, false);
    let mut primes = vec![];

    for i in 2..=n {
        if !mark[i as usize] {
            primes.push(i);
        }
        for &p in &primes {
            let t = p * i;
            if t > n {
                break;
            }
            mark.set(t as usize, true);
            if i as u64 % p == 0 {
                break;
            }
        }
    }

    primes
}

/// sieve all primes number in [l, r), using primes in `primes`.
#[inline(never)]
pub(crate) fn sieve(primes: &[u64], l: u64, r: u64) -> Vec<u64> {
    let l = l | 1;
    let n = r - l;
    info!("sieve primes in [{}, {}), interval length = {}", l, r, n);

    // reserve a little bit more space to avoid complicated boundary computation
    let mut mark = bit_vec::BitVec::from_elem((n as usize >> 1) + 50, false);
    for &p in primes {
        if p <= 5 {
            continue;
        }
        // compute the first index
        // the `| 1` here finds the first odd one
        let mut offset = (std::cmp::max((l - 1) / p + 1, 2) | 1) * p - l;

        // multipliers of p
        // we want to save half of space, by only saving the status of odd numbers.
        let p2 = p + p;
        let p4 = p2 + p2;
        let p6 = p2 + p4;
        let p30 = p * 30;

        // enable wheel sieve for small p
        if n >= offset + 2 * p30 {
            assert!((offset + l) % 2 == 1 && (offset + l) % p == 0);
            // original sieve until the start of the wheel
            while (offset + l) % p30 != p {
                // assert!(offset <= n, "p = {}, n = {}, offset = {}, l = {}", p, n, offset, l);
                mark.set(offset as usize >> 1, true);
                offset += p2;
            }

            // the largest number x, such that x <= r - 30p and x % 30p = 1
            let max_x = r - (r - 1) % p30 - p30;
            let max_offset = max_x - l;

            // the 2/3/5 wheel 1 -> 7 -> 11 -> 13 -> 17 -> 19 -> 23 -> 29 -> 31
            while offset <= max_offset {
                // assert!((offset + l) % p30 == p, "p = {}, l = {}, offset = {}", p, l, offset);
                mark.set(offset as usize >> 1, true);
                offset += p6;
                mark.set(offset as usize >> 1, true);
                offset += p4;
                mark.set(offset as usize >> 1, true);
                offset += p2;
                mark.set(offset as usize >> 1, true);
                offset += p4;
                mark.set(offset as usize >> 1, true);
                offset += p2;
                mark.set(offset as usize >> 1, true);
                offset += p4;
                mark.set(offset as usize >> 1, true);
                offset += p6;
                mark.set(offset as usize >> 1, true);
                offset += p2;
            }
        }

        // the last ~30 numbers
        while offset < n {
            mark.set(offset as usize >> 1, true);
            offset += p2;
        }
    }
    let mut ret = vec![];
    for (idx, &s) in mark.storage().iter().enumerate() {
        let mut s = !s;
        while s != 0 {
            let w = s & (s - 1);
            let offset = (s ^ w).trailing_zeros();
            s = w;

            // note that we skipped all odd numbers;
            let x = l + ((idx as u64) << 6) + ((offset as u64) << 1);
            if x % 3 != 0 && x % 5 != 0 && x <= r {
                ret.push(x);
            }
        }
    }
    info!("found {} primes in the interval", ret.len());

    ret
}

#[inline(never)]
pub(crate) fn sieve2(primes: &[u64], l: u64, r: u64) -> Vec<u64> {
    let l = l | 1;
    let ub = (r - l) as usize;
    let mut mark = bit_vec::BitVec::from_elem(ub + 1, false);
    for &p in primes {
        if p == 2 {
            continue;
        }
        let mut x = std::cmp::max((l - 1) / p + 1, 2) * p;
        if x % 2 == 0 {
            x += p;
        }
        let mut y = (x - l) as usize;
        while y <= ub {
            mark.set(y, true);
            y += 2 * p as usize;
        }
    }
    let mut ret = vec![];
    for (idx, &s) in mark.storage().iter().enumerate() {
        let mut s = !s;
        s &= 0x55555555u32;
        while s != 0 {
            let w = s & (s - 1);
            let offset = (s ^ w).trailing_zeros();
            let x = l + ((idx as u64) << 5) + offset as u64;
            if x <= r {
                ret.push(x);
            }
            s = w;
        }
    }
    info!("found {} primes in [{}, {}]", ret.len(), l, r);

    ret
}
