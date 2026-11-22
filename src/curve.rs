//! One curve `(d, r)`: walk `z = r + step*q` over the requested window and run
//! the filter chain (DESIGN §4.4) followed by the exact test (§2.3).
//!
//! The sign of `s = x + y` matters for both the mod-243 wheel and the
//! quadratic-residue masks, and it is forced by `sign(s) = sign(k - z^3)`.
//! Since `|z|` is bounded away from `k^(1/3)` in every region we scan, the sign
//! is constant on each half of the window — but it is resolved exactly here
//! rather than assumed.

use crate::exact::{self, Triple};
use crate::filters::{Filters, WHEEL};
use crate::roots::invmod;

pub const MAX_PRIMES: usize = 15;

#[derive(Clone, Copy, Debug, Default)]
pub struct Stats {
    /// Size of the full `q` range (before the wheel).
    pub q_total: u64,
    /// Points the mod-243 wheel never visits.
    pub wheel_skipped: u64,
    pub parity_rejected: u64,
    pub prime_rejected: u64,
    pub exact_tests: u64,
    pub curves: u64,
    pub no_roots: u64,
    pub dense_skipped: u64,
    /// `d` values skipped because `3 | d` (provably empty when `k ≡ 3,6 mod 9`).
    pub div3_skipped: u64,
    pub hits: u64,
}

impl Stats {
    pub fn merge(&mut self, o: &Stats) {
        self.q_total += o.q_total;
        self.wheel_skipped += o.wheel_skipped;
        self.parity_rejected += o.parity_rejected;
        self.prime_rejected += o.prime_rejected;
        self.exact_tests += o.exact_tests;
        self.curves += o.curves;
        self.no_roots += o.no_roots;
        self.dense_skipped += o.dense_skipped;
        self.div3_skipped += o.div3_skipped;
        self.hits += o.hits;
    }

    /// Points actually tested (the right denominator for throughput).
    pub fn points(&self) -> u64 {
        self.q_total - self.wheel_skipped
    }

    /// Every point is accounted for exactly once.
    pub fn balanced(&self) -> bool {
        self.q_total == self.wheel_skipped + self.parity_rejected + self.prime_rejected + self.exact_tests
    }
}

#[inline]
pub fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

#[inline]
fn div_ceil(a: i128, b: i128) -> i128 {
    debug_assert!(b > 0);
    -((-a).div_euclid(b))
}

/// Largest `c >= 0` with `c^3 <= n`.
pub fn icbrt(n: u64) -> u64 {
    if n == 0 {
        return 0;
    }
    let mut r = (n as f64).cbrt() as u64;
    while (r + 1).saturating_mul(r + 1).saturating_mul(r + 1) <= n {
        r += 1;
    }
    while r.saturating_mul(r).saturating_mul(r) > n {
        r -= 1;
    }
    r
}

/// Scan both signs of `z` with `|z|` in `[zlo, zhi]`.
///
/// `step` is `d' = d / gcd(d, 4)`, the modulus of the root `r`.
#[allow(clippy::too_many_arguments)]
pub fn scan_all(
    f: &Filters,
    k: i64,
    d: u64,
    r: u64,
    step: u64,
    zlo: u128,
    zhi: u128,
    st: &mut Stats,
    hits: &mut Vec<Triple>,
) {
    debug_assert!(step > 0 && r < step && k > 0);
    debug_assert_eq!(
        crate::roots::mulmod(crate::roots::mulmod(r, r, step), r, step),
        (k as u64) % step,
        "r is not a cube root of k modulo d'"
    );
    if zhi < zlo {
        return;
    }
    let cbrt = icbrt(k as u64) as u128;

    // Negative z: z^3 < 0, so k - z^3 > 0 and s > 0.
    if zlo > 0 {
        scan_signed(f, k, d, r, step, -(zhi as i128), -(zlo as i128), false, st, hits);
    }
    // Positive z with z^3 <= k: s >= 0 (s = 0 has no solution and is rejected).
    let hi0 = zhi.min(cbrt);
    let lo0 = zlo.max(1);
    if hi0 >= lo0 {
        scan_signed(f, k, d, r, step, lo0 as i128, hi0 as i128, false, st, hits);
    }
    // Positive z with z^3 > k: s < 0.
    let lo1 = zlo.max(cbrt + 1);
    if zhi >= lo1 {
        scan_signed(f, k, d, r, step, lo1 as i128, zhi as i128, true, st, hits);
    }
}

#[allow(clippy::too_many_arguments)]
fn scan_signed(
    f: &Filters,
    k: i64,
    d: u64,
    r: u64,
    step: u64,
    za: i128,
    zb: i128,
    s_neg: bool,
    st: &mut Stats,
    hits: &mut Vec<Triple>,
) {
    if za > zb {
        return;
    }
    let step_i = step as i128;
    let ri = r as i128;
    let qlo = div_ceil(za - ri, step_i);
    let qhi = (zb - ri).div_euclid(step_i);
    if qlo > qhi {
        return;
    }

    let total = (qhi - qlo + 1) as u64;
    st.q_total += total;
    st.curves += 1;

    // Split the q range into the progressions allowed by the mod-243 wheel.
    let g = gcd(step, WHEEL);
    let m = WHEEL / g; // q step inside one residue class
    let dm = step / g;
    let inv = if m > 1 {
        invmod(dm % m, m).expect("gcd(d', 243/g) = 1")
    } else {
        0
    };

    let smod = if s_neg {
        ((WHEEL - d % WHEEL) % WHEEL) as usize
    } else {
        (d % WHEEL) as usize
    };
    let r243 = ri.rem_euclid(WHEEL as i128);

    let nprimes = f.prime_count().min(MAX_PRIMES);
    let mut s_mods = [0u64; MAX_PRIMES];
    let mut dz_mods = [0u64; MAX_PRIMES];
    for j in 0..nprimes {
        let p = f.prime(j);
        let dmod = d % p;
        s_mods[j] = if s_neg { (p - dmod) % p } else { dmod };
        dz_mods[j] = ((step as u128 * m as u128) % p as u128) as u64;
    }

    // k - s - z must be even, and parity(s) = parity(d).
    let par_want = (k as i128 - d as i128) & 1;

    let allowed = f.allowed243[smod].clone();
    let mut visited: u64 = 0;

    for a in allowed {
        let a = a as i128;
        let diff = (a - r243).rem_euclid(WHEEL as i128);
        if diff % g as i128 != 0 {
            continue;
        }
        let q_res = if m == 1 {
            0
        } else {
            (((diff / g as i128) as u64 % m) * inv) % m
        };
        let first = qlo + (q_res as i128 - qlo).rem_euclid(m as i128);

        let z0 = ri + step_i * first;
        let mut zmods = [0u64; MAX_PRIMES];
        for j in 0..nprimes {
            zmods[j] = z0.rem_euclid(f.prime(j) as i128) as u64;
        }

        let mut q = first;
        while q <= qhi {
            let z = ri + step_i * q;
            visited += 1;

            if (z & 1) == par_want {
                let mut ok = true;
                for j in 0..nprimes {
                    if !f.pass_prime(j, s_mods[j], zmods[j]) {
                        ok = false;
                        break;
                    }
                }
                if ok {
                    st.exact_tests += 1;
                    if let Some(t) = exact::solve(k, d, z) {
                        st.hits += 1;
                        hits.push(t);
                    }
                } else {
                    st.prime_rejected += 1;
                }
            } else {
                st.parity_rejected += 1;
            }

            for j in 0..nprimes {
                let p = f.prime(j);
                let v = zmods[j] + dz_mods[j];
                zmods[j] = if v >= p { v - p } else { v };
            }
            q += m as i128;
        }
    }

    st.wheel_skipped += total - visited;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::roots;

    /// End-to-end: with the true factorisation of d, a narrow window around the
    /// published z must be enough to find the solution.
    #[test]
    fn finds_known_solutions() {
        // (k, d, factorisation, z as published, i.e. signed)
        let cases: [(i64, u64, Vec<(u64, u32)>, i128); 4] = [
            (3, 108398887211, vec![(167, 1), (649095133, 1)], -472715493453327032),
            (33, 87723532425289, vec![(87723532425289, 1)], -2736111468807040),
            (
                42,
                102980666258459,
                vec![(11, 1), (43, 1), (215921, 1), (1008323, 1)],
                12602123297335631,
            ),
            (
                795,
                21083965616656,
                vec![(2, 4), (17, 1), (149, 1), (520232077, 1)],
                2337348783323923,
            ),
        ];
        for (k, d, mut factors, z) in cases {
            for f in factors.iter_mut() {
                if f.0 == 2 {
                    f.1 = f.1.saturating_sub(2);
                }
            }
            factors.retain(|x| x.1 > 0);
            let (modulus, roots_of_k) =
                roots::roots_from_factorization(k as u64, &factors).expect("roots must exist");
            assert_eq!(modulus, d / gcd(d, 4), "modulus must be d'");

            let filters = Filters::new(k);
            let mut st = Stats::default();
            let mut hits = Vec::new();
            let window = (modulus as u128) * 2;
            let az = z.unsigned_abs();
            let zlo = az - window;
            let zhi = az + window;
            for r in roots_of_k.iter() {
                scan_all(&filters, k, d, *r, modulus, zlo, zhi, &mut st, &mut hits);
            }
            assert!(st.balanced(), "point accounting is inconsistent for k={}", k);
            assert!(
                hits.iter().any(|t| t.z == z && exact::verify(k, t.x, t.y, t.z)),
                "known solution not found for k={} ({} hits)",
                k,
                hits.len()
            );
        }
    }

    #[test]
    fn wheel_step_handles_divisible_by_three() {
        // d' divisible by 3 exercises the lcm branch of the wheel. For k ≡ 6
        // (mod 9) that only happens for other k, because 3 ∤ d is forced (see
        // DESIGN §2.4); k = 2 with d = 3 covers the branch.
        // k = 8 is a cube mod 9, which is exactly when s ≡ 0 (mod 3) — hence
        // 3 | d — can be feasible at all.
        let f = Filters::new(8);
        let mut all = Stats::default();
        for d in [3u64, 9, 15, 45] {
            let mut facs: Vec<(u64, u32)> = Vec::new();
            let mut m = d;
            let mut p = 2u64;
            while p * p <= m {
                if m % p == 0 {
                    let mut c = 0;
                    while m % p == 0 {
                        m /= p;
                        c += 1;
                    }
                    facs.push((p, c));
                }
                p += 1;
            }
            if m > 1 {
                facs.push((m, 1));
            }
            let (modulus, roots_of_k) =
                roots::roots_from_factorization(8, &facs).expect("roots must exist");
            assert_eq!(modulus, d);
            let mut st = Stats::default();
            let mut hits = Vec::new();
            for r in roots_of_k {
                scan_all(&f, 8, d, r, modulus, 1_000, 4_000, &mut st, &mut hits);
            }
            assert!(st.balanced(), "d={}", d);
            for t in hits {
                assert!(exact::verify(8, t.x, t.y, t.z));
            }
            all.merge(&st);
        }
        // Not every root class of every d need be feasible; the union must be.
        assert!(all.points() > 0, "the 3 | d branch produced no points at all");
        assert!(all.balanced());
    }

    #[test]
    fn three_dividing_d_is_dead_for_k_equiv_3_or_6_mod_9() {
        // x^3+y^3 ≡ 0 (mod 9) whenever x+y ≡ 0 (mod 3), while k ≡ 3,6 (mod 9)
        // is not a cube residue mod 9. So every d divisible by 3 is empty.
        let f = Filters::new(114);
        for d in [3u64, 6, 9, 12, 243, 729] {
            let mut st = Stats::default();
            let mut hits = Vec::new();
            scan_all(&f, 114, d, 0, d, 1_000, 4_000, &mut st, &mut hits);
            assert!(st.balanced(), "d={}", d);
            assert_eq!(st.points(), 0, "d={} should be provably empty", d);
            assert!(hits.is_empty());
        }
    }

    #[test]
    fn icbrt_is_exact() {
        for n in [0u64, 1, 7, 8, 9, 26, 27, 28, 1_000_000, 1_000_000_000] {
            let r = icbrt(n);
            assert!(r * r * r <= n);
            assert!((r + 1).checked_pow(3).map_or(true, |v| v > n));
        }
    }
}
