//! Minimal fixed-width unsigned 256-bit arithmetic.
//!
//! Only what the exact acceptance test needs: multiply, divide by a `u128`,
//! compare, and an integer square root. Everything here is on the cold path
//! (roughly one call per 5e5 quotient points), so clarity beats micro-tuning.

use std::cmp::Ordering;

#[derive(Clone, Copy, PartialEq, Eq, Default)]
pub struct U256(pub [u64; 4]); // little-endian limbs

impl std::fmt::Debug for U256 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "U256({:016x}_{:016x}_{:016x}_{:016x})",
            self.0[3], self.0[2], self.0[1], self.0[0]
        )
    }
}

impl U256 {
    pub const ZERO: U256 = U256([0, 0, 0, 0]);
    pub const ONE: U256 = U256([1, 0, 0, 0]);

    #[inline]
    pub const fn from_u64(v: u64) -> Self {
        U256([v, 0, 0, 0])
    }

    #[inline]
    pub const fn from_u128(v: u128) -> Self {
        U256([v as u64, (v >> 64) as u64, 0, 0])
    }

    #[inline]
    pub fn is_zero(self) -> bool {
        self.0 == [0, 0, 0, 0]
    }

    #[inline]
    pub fn to_u128(self) -> Option<u128> {
        if (self.0[2] | self.0[3]) != 0 {
            None
        } else {
            Some(self.0[0] as u128 | ((self.0[1] as u128) << 64))
        }
    }

    #[inline]
    pub fn bit_len(self) -> u32 {
        for i in (0..4).rev() {
            if self.0[i] != 0 {
                return 64 * (i as u32) + (64 - self.0[i].leading_zeros());
            }
        }
        0
    }

    #[inline]
    fn bit(self, i: u32) -> u64 {
        (self.0[(i >> 6) as usize] >> (i & 63)) & 1
    }

    #[inline]
    fn set_bit(&mut self, i: u32) {
        self.0[(i >> 6) as usize] |= 1u64 << (i & 63);
    }

    #[inline]
    pub fn add(self, o: Self) -> Self {
        let mut r = [0u64; 4];
        let mut carry = 0u128;
        for i in 0..4 {
            let s = self.0[i] as u128 + o.0[i] as u128 + carry;
            r[i] = s as u64;
            carry = s >> 64;
        }
        debug_assert_eq!(carry, 0, "U256 addition overflow");
        U256(r)
    }

    /// Panics in debug builds if `self < o`.
    #[inline]
    pub fn sub(self, o: Self) -> Self {
        let mut r = [0u64; 4];
        let mut borrow = 0i128;
        for i in 0..4 {
            let s = self.0[i] as i128 - o.0[i] as i128 - borrow;
            if s < 0 {
                r[i] = (s + (1i128 << 64)) as u64;
                borrow = 1;
            } else {
                r[i] = s as u64;
                borrow = 0;
            }
        }
        debug_assert_eq!(borrow, 0, "U256 subtraction underflow");
        U256(r)
    }

    #[inline]
    pub fn mul_u64(self, m: u64) -> Self {
        let mut r = [0u64; 4];
        let mut carry = 0u128;
        for i in 0..4 {
            let s = self.0[i] as u128 * m as u128 + carry;
            r[i] = s as u64;
            carry = s >> 64;
        }
        debug_assert_eq!(carry, 0, "U256 scalar multiplication overflow");
        U256(r)
    }

    /// Truncated product; a debug build asserts that no bits were lost.
    #[inline]
    pub fn mul(self, o: Self) -> Self {
        let mut wide = [0u64; 8];
        let mut lost = false;
        for i in 0..4 {
            if self.0[i] == 0 {
                continue;
            }
            let mut carry = 0u128;
            for j in 0..4 {
                let cur = wide[i + j] as u128 + (self.0[i] as u128) * (o.0[j] as u128) + carry;
                wide[i + j] = cur as u64;
                carry = cur >> 64;
            }
            let mut k = i + 4;
            while carry != 0 && k < 8 {
                let cur = wide[k] as u128 + carry;
                wide[k] = cur as u64;
                carry = cur >> 64;
                k += 1;
            }
            if carry != 0 {
                lost = true;
            }
        }
        debug_assert_eq!(&wide[4..], &[0, 0, 0, 0], "U256 multiplication overflow");
        debug_assert!(!lost, "U256 multiplication overflow");
        U256([wide[0], wide[1], wide[2], wide[3]])
    }

    #[inline]
    pub fn shl(self, bits: u32) -> Self {
        debug_assert!(bits < 256);
        let mut r = [0u64; 4];
        let limb = (bits / 64) as usize;
        let rem = bits % 64;
        for i in limb..4 {
            let src = i - limb;
            let mut v = self.0[src] << rem;
            if rem > 0 && src > 0 {
                v |= self.0[src - 1] >> (64 - rem);
            }
            r[i] = v;
        }
        U256(r)
    }

    #[inline]
    pub fn shr(self, bits: u32) -> Self {
        debug_assert!(bits < 256);
        let mut r = [0u64; 4];
        let limb = (bits / 64) as usize;
        let rem = bits % 64;
        for i in 0..(4 - limb) {
            let src = i + limb;
            let mut v = self.0[src] >> rem;
            if rem > 0 && src + 1 < 4 {
                v |= self.0[src + 1] << (64 - rem);
            }
            r[i] = v;
        }
        U256(r)
    }

    /// Binary long division. The running remainder is kept in a `u128`, which is
    /// safe because every divisor used here is `3|s| < 2^63`.
    pub fn div_rem_u128(self, d: u128) -> (Self, u128) {
        assert!(d > 0, "division by zero");
        debug_assert!(d < (1u128 << 127), "divisor too wide for this routine");
        let mut q = U256::ZERO;
        let mut rem: u128 = 0;
        for i in (0..256u32).rev() {
            rem = (rem << 1) | self.bit(i) as u128;
            if rem >= d {
                rem -= d;
                q.set_bit(i);
            }
        }
        (q, rem)
    }

    /// Integer square root by the restoring (digit-by-digit) method.
    /// The root must fit in a `u128`, which holds for every value we feed it.
    pub fn isqrt(self) -> u128 {
        let mut rem = U256::ZERO;
        let mut root: u128 = 0;
        for i in (0..128u32).rev() {
            rem = rem.shl(2).add(U256::from_u64(self.shr(2 * i).0[0] & 3));
            debug_assert!(root < (1u128 << 120), "U256::isqrt root too wide");
            let trial = U256::from_u128((root << 2) | 1);
            if rem.cmp(&trial) != Ordering::Less {
                rem = rem.sub(trial);
                root = (root << 1) | 1;
            } else {
                root <<= 1;
            }
        }
        root
    }

    #[inline]
    pub fn lt(self, o: Self) -> bool {
        self.cmp(&o) == Ordering::Less
    }
}

impl Ord for U256 {
    fn cmp(&self, o: &Self) -> Ordering {
        for i in (0..4).rev() {
            match self.0[i].cmp(&o.0[i]) {
                Ordering::Equal => {}
                other => return other,
            }
        }
        Ordering::Equal
    }
}

impl PartialOrd for U256 {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip_u128() {
        let v = 123456789012345678901234567890u128;
        let u = U256::from_u128(v);
        assert_eq!(u.to_u128(), Some(v));
        assert_eq!(u.bit_len(), 128 - v.leading_zeros());
    }

    #[test]
    fn mul_matches_u128() {
        for &(a, b) in &[
            (1u128, 1u128),
            (1_000_000_007, 1_000_000_009),
            (u32::MAX as u128, u32::MAX as u128),
            (1u128 << 60, 1u128 << 60),
        ] {
            let p = U256::from_u128(a).mul(U256::from_u128(b));
            assert_eq!(p.to_u128(), Some(a * b), "{} * {}", a, b);
        }
    }

    #[test]
    fn div_rem_matches_u128() {
        for &(n, d) in &[
            (100u128, 7u128),
            (1u128 << 100, 3),
            (123456789012345678901234567890, 987654321),
            (5, 5),
        ] {
            let (q, r) = U256::from_u128(n).div_rem_u128(d);
            assert_eq!(q.to_u128(), Some(n / d), "{} / {}", n, d);
            assert_eq!(r, n % d, "{} % {}", n, d);
        }
    }

    #[test]
    fn isqrt_small() {
        for n in 0u128..2000 {
            let r = U256::from_u128(n).isqrt();
            assert!(r * r <= n && (r + 1) * (r + 1) > n, "isqrt({}) = {}", n, r);
        }
    }

    #[test]
    fn isqrt_large() {
        // (10^23)^2 and neighbours
        let base: u128 = 100_000_000_000_000_000_000_000;
        let sq = U256::from_u128(base).mul(U256::from_u128(base));
        assert_eq!(sq.isqrt(), base);
        let sq1 = sq.add(U256::from_u64(1));
        assert_eq!(sq1.isqrt(), base);
        let sqm1 = sq.sub(U256::ONE);
        assert_eq!(sqm1.isqrt(), base - 1);
    }

    #[test]
    fn shifts() {
        let v = U256::from_u64(0b1011);
        assert_eq!(v.shl(4).to_u128(), Some(0b1011_0000));
        assert_eq!(v.shl(70).to_u128(), Some(0b1011u128 << 70));
        assert_eq!(v.shl(200).to_u128(), None);
        assert_eq!(v.shl(70).shr(70), v);
    }
}
