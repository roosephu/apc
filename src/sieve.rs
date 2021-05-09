use num::ToPrimitive;
use primesieve_sys::*;

pub struct PrimesieveResult<'a> {
    pub primes: &'a [u64],
    raw_ptr: *mut libc::c_void,
}

impl<'a> PrimesieveResult<'a> {
    pub fn new(raw_ptr: *mut libc::c_void, size: libc::size_t) -> Self {
        let size = size.to_usize().unwrap();
        let primes = unsafe { std::slice::from_raw_parts(raw_ptr as *mut u64, size) };
        Self { primes, raw_ptr }
    }
}

impl<'a> Drop for PrimesieveResult<'a> {
    fn drop(&mut self) {
        unsafe {
            primesieve_sys::primesieve_free(self.raw_ptr);
        }
    }
}

// https://github.com/pthariensflame/primesieve.rs/blob/master/primesieve-rs/src/lib.rs#L422-L437
#[inline(never)]
pub fn sieve_primesieve<'a>(start: u64, stop: u64) -> PrimesieveResult<'a> {
    let mut size: libc::size_t = 0;
    unsafe {
        primesieve_set_num_threads(1);
        let raw_ptr = primesieve_generate_primes(start, stop, &mut size, INT64_PRIMES);
        PrimesieveResult::new(raw_ptr, size)
    }
}
