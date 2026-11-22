//! Modular filters: the mod-243 wheel and the small-prime quadratic-residue bitmaps.
//!
//! Both are *necessary* conditions for `v^2 = (4(k-z^3) - s^3)/(3s)` to be a
//! perfect square, so rejecting on them never loses a solution. Together they
//! remove about `1 - 1/5e5` of the candidates (measured on the upstream
//! campaign's counters), which is the entire reason the hot loop is cheap.

pub const WHEEL: u64 = 243;

/// Primes used for the quadratic-residue stage. All are distinct from 2, 3, 19
/// and 114's other prime factors are irrelevant: a vacuous test is skipped.
pub const QR_PRIMES: [u64; 15] = [5, 7, 11, 13, 17, 19, 23, 31, 37, 41, 43, 47, 53, 59, 61];

pub struct Filters {
    pub k: i64,
    /// `allowed243[s]` = admissible residues of `z` mod 243 for `s = x + y`.
    pub allowed243: Vec<Vec<u8>>,
    /// `masks[j][s mod p]` is a bitmap over `z mod p` for `QR_PRIMES[j]`.
    masks: Vec<Vec<u64>>,
}

impl Filters {
    pub fn new(k: i64) -> Self {
        let k243 = k.rem_euclid(WHEEL as i64) as u64;

        let mut cubes = [0u64; WHEEL as usize];
        for v in 0..WHEEL {
            cubes[v as usize] = v * v % WHEEL * v % WHEEL;
        }

        // pairs[s * 243 + r] = "some (x,y) has x+y = s and x^3+y^3 = r mod 243"
        let mut pairs = vec![false; (WHEEL * WHEEL) as usize];
        for x in 0..WHEEL {
            for y in 0..WHEEL {
                let s = (x + y) % WHEEL;
                let r = (cubes[x as usize] + cubes[y as usize]) % WHEEL;
                pairs[(s * WHEEL + r) as usize] = true;
            }
        }

        let mut allowed243 = vec![Vec::new(); WHEEL as usize];
        for s in 0..WHEEL as usize {
            let mut v = Vec::new();
            for z in 0..WHEEL as usize {
                let want = ((k243 + WHEEL - cubes[z]) % WHEEL) as usize;
                if pairs[s * WHEEL as usize + want] {
                    v.push(z as u8);
                }
            }
            allowed243[s] = v;
        }

        let mut masks = Vec::with_capacity(QR_PRIMES.len());
        for &p in QR_PRIMES.iter() {
            let pc = p as i64;
            let kp = k.rem_euclid(pc);
            let mut is_qr = vec![false; p as usize];
            for v in 0..pc {
                is_qr[((v * v) % pc) as usize] = true;
            }
            let mut m = vec![0u64; p as usize];
            for sm in 0..pc {
                let mut bits = 0u64;
                for zm in 0..pc {
                    let z3 = zm * zm % pc * zm % pc;
                    let inner =
                        (4 * ((kp - z3).rem_euclid(pc)) - (sm * sm % pc) * sm % pc).rem_euclid(pc);
                    // 3s(4(k-z^3) - s^3) must be a QR mod p; equivalent to the
                    // quotient being a square because (3s)^{-2} is itself a square.
                    let val = (3 * sm % pc) * inner % pc;
                    if is_qr[val as usize] {
                        bits |= 1u64 << zm;
                    }
                }
                m[sm as usize] = bits;
            }
            masks.push(m);
        }

        Filters {
            k,
            allowed243,
            masks,
        }
    }

    #[inline]
    pub fn prime_count(&self) -> usize {
        self.masks.len()
    }

    #[inline]
    pub fn prime(&self, j: usize) -> u64 {
        QR_PRIMES[j]
    }

    #[inline]
    pub fn pass_prime(&self, j: usize, s_mod: u64, z_mod: u64) -> bool {
        (self.masks[j][s_mod as usize] >> z_mod) & 1 == 1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn brute_allowed(k: i64, s: u64, z: u64) -> bool {
        for x in 0..WHEEL {
            let y = (s + WHEEL * 2 - x) % WHEEL; // x + y = s (mod 243)
            let r = ((x * x % WHEEL) * x + (y * y % WHEEL) * y) % WHEEL;
            let want = (k.rem_euclid(243) as u64 + WHEEL - (z * z % WHEEL) * z % WHEEL) % WHEEL;
            if r % WHEEL == want {
                return true;
            }
        }
        false
    }

    #[test]
    fn wheel_matches_brute_force() {
        let f = Filters::new(114);
        for s in 0..WHEEL as u64 {
            for z in 0..WHEEL as u64 {
                let in_table = f.allowed243[s as usize].contains(&(z as u8));
                assert_eq!(in_table, brute_allowed(114, s, z), "s={} z={}", s, z);
            }
        }
    }

    #[test]
    fn known_solutions_pass_every_filter() {
        let cases: [(i64, i128, i128, i128); 4] = [
            (3, 569936821221962380720, -569936821113563493509, -472715493453327032),
            (33, 8866128975287528, -8778405442862239, -2736111468807040),
            (42, -80538738812075974, 80435758145817515, 12602123297335631),
            (795, -14219049725358227, 14197965759741571, 2337348783323923),
        ];
        for &(k, x, y, z) in cases.iter() {
            let f = Filters::new(k);
            let s = x + y;
            assert!(
                f.allowed243[s.rem_euclid(243) as usize].contains(&(z.rem_euclid(243) as u8)),
                "wheel rejected a real solution for k={}",
                k
            );
            for j in 0..f.prime_count() {
                let p = f.prime(j) as i128;
                assert!(
                    f.pass_prime(j, s.rem_euclid(p) as u64, z.rem_euclid(p) as u64),
                    "prime {} rejected a real solution for k={}",
                    p,
                    k
                );
            }
        }
    }

    #[test]
    fn wheel_is_actually_selective() {
        let f = Filters::new(114);
        let total: usize = 243 * 243;
        let kept: usize = f.allowed243.iter().map(|v| v.len()).sum();
        assert!(
            (kept as f64) < 0.2 * total as f64,
            "wheel keeps {} of {} pairs — suspiciously weak",
            kept,
            total
        );
    }
}
