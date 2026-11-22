//! Enumeration and factorisation of `d` over a block of a contiguous interval.
//!
//! Segmented factorisation: sieve the block with every prime up to `sqrt(hi)`,
//! dividing out prime powers. Whatever survives is prime (a composite survivor
//! would have a factor <= sqrt(hi), already divided out). This keeps coverage
//! complete — unlike a smoothness-bounded DFS, which cannot reach `d` that is
//! prime near 1e18 (and `d` for the k = 33 solution *is* prime).

pub const MAX_FACTORS: usize = 16;
/// Marker in `nfac` for "more distinct primes than we can store".
pub const NFAC_OVERFLOW: u8 = 255;

#[derive(Clone, Copy)]
pub struct DEntry {
    pub d: u64,
    pub nfac: u8,
    pub p: [u64; MAX_FACTORS],
    pub e: [u8; MAX_FACTORS],
}

impl DEntry {
    fn new(d: u64) -> Self {
        DEntry {
            d,
            nfac: 0,
            p: [0; MAX_FACTORS],
            e: [0; MAX_FACTORS],
        }
    }

    #[inline]
    pub fn push(&mut self, p: u64, e: u32) {
        let i = self.nfac as usize;
        if i < MAX_FACTORS && e <= u8::MAX as u32 {
            self.p[i] = p;
            self.e[i] = e as u8;
            self.nfac += 1;
        } else {
            self.nfac = NFAC_OVERFLOW;
        }
    }
}

/// Odd-only prime bitmap.
pub struct PrimeBits {
    limit: u64,
    bits: Vec<u64>, // bit i set  <=>  2i + 3 is prime
}

impl PrimeBits {
    pub fn new(limit: u64) -> Self {
        let n_odd = if limit < 3 { 0 } else { (limit - 3) / 2 + 1 };
        let words = (n_odd as usize + 63) / 64;
        let mut bits = vec![u64::MAX; words];
        if n_odd as usize % 64 != 0 && words > 0 {
            let keep = n_odd as usize % 64;
            bits[words - 1] = (1u64 << keep) - 1;
        }
        let root = isqrt_u64(limit);
        let mut i: u64 = 0;
        while 2 * i + 3 <= root {
            let p = 2 * i + 3;
            if (bits[(i / 64) as usize] >> (i % 64)) & 1 == 1 {
                let mut m = p * p;
                while m <= limit {
                    let idx = (m - 3) / 2;
                    bits[(idx / 64) as usize] &= !(1u64 << (idx % 64));
                    m += 2 * p;
                }
            }
            i += 1;
        }
        PrimeBits { limit, bits }
    }

    pub fn limit(&self) -> u64 {
        self.limit
    }

    pub fn for_each(&self, mut f: impl FnMut(u64)) {
        f(2);
        for (wi, &w0) in self.bits.iter().enumerate() {
            let mut w = w0;
            while w != 0 {
                let b = w.trailing_zeros() as u64;
                let idx = wi as u64 * 64 + b;
                let n = 2 * idx + 3;
                if n <= self.limit {
                    f(n);
                }
                w &= w - 1;
            }
        }
    }
}

#[inline]
pub fn isqrt_u64(n: u64) -> u64 {
    if n < 2 {
        return n;
    }
    let mut x = (n as f64).sqrt() as u64;
    // Correct the floating-point estimate with exact checks.
    while x.saturating_mul(x) > n {
        x -= 1;
    }
    while (x + 1).saturating_mul(x + 1) <= n {
        x += 1;
    }
    x
}

/// Factorise every `d` in `[lo, hi)`. Returns one entry per `d`, in order.
pub fn factor_block(lo: u64, hi: u64, sieve: &PrimeBits) -> Vec<DEntry> {
    assert!(lo < hi);
    let n = (hi - lo) as usize;
    let mut rem: Vec<u64> = (lo..hi).collect();
    let mut entries: Vec<DEntry> = (lo..hi).map(DEntry::new).collect();

    sieve.for_each(|p| {
        let p2 = p * p;
        let mut m = if p2 >= lo { p2 } else { ((lo + p - 1) / p) * p };
        while m < hi {
            let idx = (m - lo) as usize;
            let mut cnt = 0u32;
            while rem[idx] % p == 0 {
                rem[idx] /= p;
                cnt += 1;
            }
            if cnt > 0 && entries[idx].nfac != NFAC_OVERFLOW {
                entries[idx].push(p, cnt);
            }
            m += p;
        }
    });

    for idx in 0..n {
        if rem[idx] > 1 && entries[idx].nfac != NFAC_OVERFLOW {
            entries[idx].push(rem[idx], 1);
        }
    }
    entries
}

#[cfg(test)]
mod tests {
    use super::*;

    fn naive_factor(n: u64) -> Vec<(u64, u32)> {
        let mut out = Vec::new();
        let mut m = n;
        let mut p = 2u64;
        while p * p <= m {
            if m % p == 0 {
                let mut c = 0;
                while m % p == 0 {
                    m /= p;
                    c += 1;
                }
                out.push((p, c));
            }
            p += 1;
        }
        if m > 1 {
            out.push((m, 1));
        }
        out
    }

    fn check_range(lo: u64, hi: u64) {
        let sieve = PrimeBits::new(isqrt_u64(hi - 1) + 1);
        let entries = factor_block(lo, hi, &sieve);
        assert_eq!(entries.len(), (hi - lo) as usize);
        for e in entries {
            let got: Vec<(u64, u32)> = (0..e.nfac as usize).map(|i| (e.p[i], e.e[i] as u32)).collect();
            assert_eq!(got, naive_factor(e.d), "wrong factorisation for {}", e.d);
            let prod: u64 = got.iter().map(|&(p, ex)| p.pow(ex)).product();
            assert_eq!(prod, e.d);
        }
    }

    #[test]
    fn small_ranges() {
        check_range(2, 500);
        check_range(1000, 1600);
        check_range(1_000_000, 1_010_000);
    }

    #[test]
    fn prime_is_reported_as_huge_factor() {
        // 87723532425289 is prime and is d for the k = 33 solution.
        let d = 87723532425289u64;
        let sieve = PrimeBits::new(isqrt_u64(d));
        let mut has_small_factor = false;
        sieve.for_each(|p| {
            if p < d && d % p == 0 {
                has_small_factor = true;
            }
        });
        assert!(!has_small_factor, "{} must have no prime factor <= sqrt", d);
    }

    #[test]
    fn isqrt_correct() {
        for n in [0u64, 1, 2, 3, 4, 8, 9, 10, 1_000_000, u64::MAX / 4] {
            let r = isqrt_u64(n);
            assert!(r * r <= n);
            assert!((r + 1).saturating_mul(r + 1) > n);
        }
    }
}
