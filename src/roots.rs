//! Cube roots of `k` modulo prime powers, and CRT assembly.
//!
//! The search fixes `d` and walks `z = r + d'q`, where `r^3 ≡ k (mod d')` and
//! `d' = d / gcd(d, 4)`. Getting *all* roots matters: missing one silently
//! loses a whole residue class of solutions.

use std::vec::Vec;

#[inline]
pub fn mulmod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}

pub fn powmod(base: u64, exp: u64, m: u64) -> u64 {
    if m == 1 {
        return 0;
    }
    let mut acc = 1u64 % m;
    let mut b = base % m;
    let mut e = exp;
    while e > 0 {
        if e & 1 == 1 {
            acc = mulmod(acc, b, m);
        }
        b = mulmod(b, b, m);
        e >>= 1;
    }
    acc
}

/// Modular inverse via the extended Euclidean algorithm; `None` if not coprime.
pub fn invmod(a: u64, m: u64) -> Option<u64> {
    if m == 1 {
        return Some(0);
    }
    let (mut old_r, mut r) = ((a % m) as i128, m as i128);
    let (mut old_s, mut s) = (1i128, 0i128);
    while r != 0 {
        let q = old_r / r;
        let tr = old_r - q * r;
        old_r = r;
        r = tr;
        let ts = old_s - q * s;
        old_s = s;
        s = ts;
    }
    if old_r != 1 {
        return None;
    }
    Some(old_s.rem_euclid(m as i128) as u64)
}

/// All `r mod p` with `r^3 ≡ k (mod p)`.
pub fn roots_mod_prime(k: u64, p: u64) -> Vec<u64> {
    debug_assert!(p >= 2);
    let a = k % p;
    match p {
        // On F_2 and F_3 the cube map is the identity (gcd(3, p-1) = 1).
        2 | 3 => return vec![a],
        _ => {}
    }
    if a == 0 {
        return vec![0];
    }
    if p % 3 == 2 {
        // The cube map is a bijection on F_p^*.
        let u = invmod(3, p - 1).expect("gcd(3, p-1) = 1");
        return vec![powmod(a, u, p)];
    }
    cube_roots_1mod3(a, p)
}

/// Adleman–Manders–Miller for `p ≡ 1 (mod 3)`.
fn cube_roots_1mod3(a: u64, p: u64) -> Vec<u64> {
    let n = p - 1;
    if powmod(a, n / 3, p) != 1 {
        return Vec::new(); // not a cubic residue
    }

    // n = 3^t * s with 3 ∤ s. s >= 2 because p is odd and p ≡ 1 (mod 3).
    let (t, s) = {
        let mut t = 0u32;
        let mut s = n;
        while s % 3 == 0 {
            s /= 3;
            t += 1;
        }
        (t, s)
    };
    debug_assert!(t >= 1 && s >= 2);
    debug_assert!(t <= 39, "3^t must fit in u64");

    // x0 = a^u with 3u ≡ 1 (mod s); then x0^3 = a * delta with delta = (a^s)^c.
    let u = invmod(3 % s, s).expect("gcd(3, s) = 1");
    let c = (3 * u - 1) / s; // 0 <= c <= 2
    let x0 = powmod(a, u, p);

    let order3t = 3u64.pow(t); // 3^t, fits u64 by the assertion above
    let as_ = powmod(a, s, p); // lives in the 3-Sylow subgroup
    let delta = powmod(as_, c, p);
    let delta_inv = powmod(delta, p - 2, p);

    // A non-cubic residue, and the generator of the 3-Sylow.
    let h = {
        let mut h = 2u64;
        while powmod(h, n / 3, p) == 1 {
            h += 1;
        }
        h
    };
    let eta = powmod(h, s, p);

    let y = cube_root_3sylow(delta_inv, eta, t, p);
    let r = mulmod(x0, y, p);
    debug_assert_eq!(mulmod(mulmod(r, r, p), r, p), a, "cube root self-check");

    // The three roots are r * omega^j.
    let omega = powmod(eta, order3t / 3, p);
    let mut out = Vec::with_capacity(3);
    let mut cur = r;
    for _ in 0..3 {
        if !out.contains(&cur) {
            out.push(cur);
        }
        cur = mulmod(cur, omega, p);
    }
    out
}

/// Cube root of an element of the 3-Sylow subgroup, by base-3 digit extraction.
fn cube_root_3sylow(delta: u64, eta: u64, t: u32, p: u64) -> u64 {
    if t == 0 {
        return 1;
    }
    let omega = powmod(eta, 3u64.pow(t - 1), p);
    let omega2 = mulmod(omega, omega, p);
    let eta_inv = invmod(eta, p).expect("eta is a unit");

    let mut acc = delta;
    let mut e: u128 = 0; // discrete log of delta to base eta
    let mut p3: u128 = 1;
    for j in 0..t {
        let exp = 3u64.pow(t - 1 - j);
        let v = powmod(acc, exp, p);
        let digit = if v == 1 {
            0u64
        } else if v == omega {
            1
        } else if v == omega2 {
            2
        } else {
            debug_assert!(false, "element is not in the 3-Sylow subgroup");
            0
        };
        e += digit as u128 * p3;
        p3 *= 3;
        if digit != 0 {
            let corr = powmod(eta_inv, digit * p3_u64(j), p);
            acc = mulmod(acc, corr, p);
        }
    }
    debug_assert_eq!(e % 3, 0, "delta was not a cube in the 3-Sylow");
    let f = (e / 3) as u64;
    let y = powmod(eta, f, p);
    debug_assert_eq!(mulmod(mulmod(y, y, p), y, p), delta, "3-Sylow root check");
    y
}

#[inline]
fn p3_u64(j: u32) -> u64 {
    3u64.pow(j)
}

/// All `r mod p^e` with `r^3 ≡ k (mod p^e)`.
pub fn roots_mod_prime_power(k: u64, p: u64, e: u32) -> Vec<u64> {
    debug_assert!(e >= 1);
    if k % p == 0 {
        // p | k. Write v = v_p(k) and z = p^w u; then v_p(z^3 - k) is
        //   3w          if 3w < v,
        //   v           if 3w > v,
        //   v + v_p(u^3 - k/p^v)  if 3w = v,
        // so for 3 ∤ v we need v >= e and w >= ceil((v+1)/3). For 3 | v we
        // enumerate directly (these moduli are tiny for the k we care about).
        let Some(pe) = p.checked_pow(e) else {
            return Vec::new();
        };
        let v = {
            let mut v = 0u32;
            let mut kk = k;
            while kk % p == 0 {
                kk /= p;
                v += 1;
            }
            v
        };
        if v % 3 != 0 {
            if e > v {
                return Vec::new();
            }
            let need = (v + 3) / 3; // ceil((v+1)/3)
            let Some(st) = p.checked_pow(need) else {
                return Vec::new();
            };
            if st >= pe {
                return vec![0];
            }
            return (0..pe).step_by(st as usize).collect();
        }
        if pe <= 1_000_000 {
            return (0..pe)
                .filter(|z| mulmod(mulmod(*z, *z, pe), *z, pe) == k % pe)
                .collect();
        }
        debug_assert!(false, "3 | v_p(k) with a large modulus is not implemented");
        return Vec::new();
    }

    let mut roots = roots_mod_prime(k % p, p);
    if roots.is_empty() {
        return roots;
    }
    let mut m: u64 = p;
    for _ in 1..e {
        let Some(mp) = m.checked_mul(p) else {
            return Vec::new();
        };
        let k_mp = k % mp;
        let mut next = Vec::with_capacity(roots.len());
        for &r in roots.iter() {
            let r3 = mulmod(mulmod(r % mp, r % mp, mp), r % mp, mp);
            let diff = (r3 + mp - k_mp) % mp;
            debug_assert_eq!(diff % m, 0, "Hensel precondition");
            let d3r2 = mulmod(3, mulmod(r % p, r % p, p), p);
            match invmod(d3r2, p) {
                Some(inv) => {
                    // Ordinary Hensel step: the derivative is a unit.
                    let q = diff / m;
                    let t = mulmod(p - q % p, inv, p);
                    next.push(r + t * m);
                }
                None => {
                    // p | 3r^2 — only for p = 3, where the derivative vanishes
                    // and each root branches three ways. Test the branches.
                    for t in 0..p {
                        let cand = r + t * m;
                        let c3 = mulmod(mulmod(cand % mp, cand % mp, mp), cand % mp, mp);
                        if c3 == k_mp {
                            next.push(cand);
                        }
                    }
                }
            }
        }
        roots = next;
        m = mp;
    }
    roots
}

/// CRT combine, for coprime moduli.
pub fn crt(r1: u64, m1: u64, r2: u64, m2: u64) -> u64 {
    let inv = invmod(m1 % m2, m2).expect("CRT moduli must be coprime");
    let diff = (r2 + m2 - r1 % m2) % m2;
    let t = (inv as u128 * diff as u128 % m2 as u128) as u64;
    let m = m1 as u128 * m2 as u128;
    ((r1 as u128 + m1 as u128 * t as u128) % m) as u64
}

/// Roots modulo `prod p^e`, given the factorisation.
///
/// Returns `Some((modulus, roots))` with `modulus = prod p^e`, or `None` when
/// `k` has no cube root modulo that product.
pub fn roots_from_factorization(k: u64, factors: &[(u64, u32)]) -> Option<(u64, Vec<u64>)> {
    let mut modulus: u64 = 1;
    let mut roots: Vec<u64> = vec![0];
    for &(p, e) in factors {
        let pe = p.checked_pow(e)?;
        let local = roots_mod_prime_power(k, p, e);
        if local.is_empty() {
            return None;
        }
        if modulus == 1 {
            modulus = pe;
            roots = local;
            continue;
        }
        let mut next = Vec::with_capacity(roots.len() * local.len());
        for &r1 in roots.iter() {
            for &r2 in local.iter() {
                next.push(crt(r1, modulus, r2, pe));
            }
        }
        modulus = modulus.checked_mul(pe)?;
        roots = next;
    }
    Some((modulus, roots))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn is_prime(n: u64) -> bool {
        if n < 2 {
            return false;
        }
        let mut d = 2u64;
        while d * d <= n {
            if n % d == 0 {
                return false;
            }
            d += 1;
        }
        true
    }

    #[test]
    fn roots_mod_prime_match_brute_force() {
        for p in 2u64..600 {
            if !is_prime(p) {
                continue;
            }
            for k in 0u64..200 {
                let mut expect: Vec<u64> = (0..p).filter(|z| (z * z % p) * z % p == k % p).collect();
                let mut got = roots_mod_prime(k, p);
                expect.sort_unstable();
                got.sort_unstable();
                assert_eq!(got, expect, "p={} k={}", p, k);
            }
        }
    }

    #[test]
    fn roots_mod_prime_large_prime() {
        // 87723532425289 | d for the k = 33 solution: a prime modulus near 1e14.
        let p = 87723532425289u64;
        assert!(is_prime(p));
        let roots = roots_mod_prime(33, p);
        assert!(!roots.is_empty());
        for r in roots.iter() {
            assert_eq!(mulmod(mulmod(*r, *r, p), *r, p), 33);
        }
        assert_eq!(roots.len(), 3, "p ≡ 1 mod 3 should give three roots");
    }

    #[test]
    fn prime_power_lifting() {
        // Check r^3 = k mod p^e for a spread of small primes.
        for &(k, p, e) in &[(114u64, 5u64, 4u32), (114, 7, 3), (114, 13, 2), (2, 7, 3), (795, 11, 2)] {
            let pe = p.pow(e);
            for r in roots_mod_prime_power(k, p, e) {
                assert_eq!(mulmod(mulmod(r, r, pe), r, pe), k % pe, "k={} p={} e={}", k, p, e);
            }
            // Completeness against brute force for one cheap case.
            if pe <= 2000 {
                let mut expect: Vec<u64> = (0..pe).filter(|z| (z * z % pe) * z % pe == k % pe).collect();
                let mut got = roots_mod_prime_power(k, p, e);
                expect.sort_unstable();
                got.sort_unstable();
                assert_eq!(got, expect, "k={} p={} e={}", k, p, e);
            }
        }
    }

    #[test]
    fn factorization_roots_are_consistent() {
        // The four known solutions' d values, with the roots they must contain.
        let cases: [(u64, Vec<(u64, u32)>, i128); 4] = [
            (3, vec![(167, 1), (649095133, 1)], -472715493453327032),
            (33, vec![(87723532425289, 1)], -2736111468807040),
            (42, vec![(11, 1), (43, 1), (215921, 1), (1008323, 1)], 12602123297335631),
            (795, vec![(2, 4), (17, 1), (149, 1), (520232077, 1)], 2337348783323923),
        ];
        for (k, factors, z) in cases {
            let (modulus, roots) =
                roots_from_factorization(k, &factors).expect("roots must exist for these d");
            // The congruence constrains the *signed* z: z^3 ≡ k (mod d').
            let zz = z.rem_euclid(modulus as i128) as u128;
            assert!(
                roots.iter().any(|&r| r as u128 == zz),
                "k={} z is not on any root class (modulus {})",
                k,
                modulus
            );
            for r in roots.iter() {
                assert_eq!(mulmod(mulmod(*r, *r, modulus), *r, modulus), k % modulus);
            }
        }
    }

    #[test]
    fn invmod_roundtrip() {
        for m in 2u64..500 {
            for a in 1u64..m {
                if let Some(i) = invmod(a, m) {
                    assert_eq!(mulmod(a, i, m), 1 % m, "a={} m={}", a, m);
                }
            }
        }
    }
}
