//! The exact acceptance test (DESIGN §2.3).
//!
//! Given `k`, `d = |x+y|` and a signed candidate `z`, decide whether a pair
//! `(x, y)` exists. Everything is exact integer arithmetic: a sign-magnitude
//! wrapper over `U256`, with the final identity re-checked from scratch.

use crate::u256::U256;
use std::cmp::Ordering;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Triple {
    pub x: i128,
    pub y: i128,
    pub z: i128,
}

impl std::fmt::Display for Triple {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}, {}, {})", self.x, self.y, self.z)
    }
}

/// Signed integer with a 256-bit magnitude: enough for every intermediate here.
#[derive(Clone, Copy, Debug)]
struct Sm {
    neg: bool,
    mag: U256,
}

impl Sm {
    const ZERO: Sm = Sm {
        neg: false,
        mag: U256::ZERO,
    };

    #[inline]
    fn new(neg: bool, mag: U256) -> Self {
        if mag.is_zero() {
            Sm::ZERO
        } else {
            Sm { neg, mag }
        }
    }

    fn from_i128(v: i128) -> Self {
        Sm::new(v < 0, U256::from_u128(v.unsigned_abs()))
    }

    #[inline]
    fn is_zero(self) -> bool {
        self.mag.is_zero()
    }

    #[inline]
    fn mul_u64(self, m: u64) -> Self {
        Sm::new(self.neg, self.mag.mul_u64(m))
    }

    fn add(self, o: Sm) -> Sm {
        if self.neg == o.neg {
            Sm::new(self.neg, self.mag.add(o.mag))
        } else if self.mag.cmp(&o.mag) == Ordering::Less {
            Sm::new(o.neg, o.mag.sub(self.mag))
        } else {
            Sm::new(self.neg, self.mag.sub(o.mag))
        }
    }

    #[inline]
    fn negate(self) -> Sm {
        Sm::new(!self.neg, self.mag)
    }

    #[inline]
    fn sub(self, o: Sm) -> Sm {
        self.add(o.negate())
    }

    #[inline]
    fn mul(self, o: Sm) -> Sm {
        Sm::new(self.neg ^ o.neg, self.mag.mul(o.mag))
    }

    fn cube_of_i128(v: i128) -> Sm {
        let s = Sm::from_i128(v);
        s.mul(s).mul(s)
    }
}

/// Solve for `(x, y)` given `k`, `d = |x + y|` and a signed `z`.
///
/// The sign of `s = x + y` is not free: `s(s^2 + 3v^2)/4 = k - z^3` forces
/// `sign(s) = sign(k - z^3)`, so it is derived here rather than passed in.
pub fn solve(k: i64, d: u64, z: i128) -> Option<Triple> {
    debug_assert!(d > 0, "d = 0 cannot yield a solution");

    let z_neg = z < 0;
    let az = z.unsigned_abs();
    let z3 = {
        let zz = U256::from_u128(az);
        Sm::new(z_neg, zz.mul(zz).mul(zz))
    };

    let t = Sm::from_i128(k as i128).sub(z3); // t = k - z^3
    if t.is_zero() {
        return None; // would force s = 0
    }

    let s = Sm::new(t.neg, U256::from_u64(d));
    let n = t.mul_u64(4).sub(s.mul(s).mul(s)); // N = 4t - s^3
    let den = s.mul_u64(3);

    let q_mag = if n.is_zero() {
        U256::ZERO
    } else {
        if n.neg != den.neg {
            return None; // v^2 < 0
        }
        let divisor = den.mag.to_u128()?;
        let (q, rem) = n.mag.div_rem_u128(divisor);
        if rem != 0 {
            return None; // 3s does not divide N
        }
        q
    };

    // v^2 = q_mag must be a perfect square.
    let v = q_mag.isqrt();
    if U256::from_u128(v).mul(U256::from_u128(v)) != q_mag {
        return None;
    }

    // x, y integers require s ≡ v (mod 2); parity of s equals parity of d.
    if (d & 1) != ((v & 1) as u64) {
        return None;
    }

    let vi = i128::try_from(v).ok()?;
    let si = if s.neg { -(d as i128) } else { d as i128 };
    let x = (si + vi) / 2;
    let y = (si - vi) / 2;

    verify(k, x, y, z).then_some(Triple { x, y, z })
}

/// Independent re-check of the identity, in exact arithmetic.
pub fn verify(k: i64, x: i128, y: i128, z: i128) -> bool {
    let sum = Sm::cube_of_i128(x)
        .add(Sm::cube_of_i128(y))
        .add(Sm::cube_of_i128(z));
    let kk = Sm::from_i128(k as i128);
    sum.neg == kk.neg && sum.mag == kk.mag
}

#[cfg(test)]
mod tests {
    use super::*;

    // The four published solutions used throughout the design document.
    const CASES: [(i64, i128, i128, i128); 4] = [
        (3, 569936821221962380720, -569936821113563493509, -472715493453327032),
        (33, 8866128975287528, -8778405442862239, -2736111468807040),
        (42, -80538738812075974, 80435758145817515, 12602123297335631),
        (795, -14219049725358227, 14197965759741571, 2337348783323923),
    ];

    #[test]
    fn known_solutions_are_accepted() {
        for &(k, x, y, z) in CASES.iter() {
            assert!(verify(k, x, y, z), "verify failed for k={}", k);
            let d = (x + y).unsigned_abs() as u64;
            let got = solve(k, d, z).unwrap_or_else(|| panic!("solve refused k={}", k));
            assert!(verify(k, got.x, got.y, got.z));
            let mut v = [got.x, got.y];
            v.sort();
            let mut w = [x, y];
            w.sort();
            assert_eq!(v, w, "wrong pair for k={}", k);
        }
    }

    #[test]
    fn neighbours_are_rejected() {
        for &(k, x, y, z) in CASES.iter() {
            let d = (x + y).unsigned_abs() as u64;
            for dz in -2i128..=2 {
                if dz == 0 {
                    continue;
                }
                if let Some(t) = solve(k, d, z + dz) {
                    panic!("false positive for k={}: {} (dz={})", k, t, dz);
                }
            }
            // A perturbed d must not produce the same triple either.
            for dd in [1u64, 7] {
                if let Some(t) = solve(k, d + dd, z) {
                    assert!(!verify(k, t.x, t.y, t.z), "unverified hit leaked out");
                }
            }
        }
    }

    #[test]
    fn zero_z_and_small_cases() {
        // 2 = 1^3 + 1^3 + 0^3
        let t = solve(2, 2, 0).expect("(1,1,0)");
        assert_eq!((t.x, t.y), (1, 1));
        // 43 = 9^3 - 7^3 - 7^3, the alpha-lemma extremal example
        let t = solve(43, 2, -7).expect("(9,-7,-7)");
        assert_eq!(t.x + t.y, 2);
    }
}
