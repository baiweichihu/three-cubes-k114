//! End-to-end checks against published solutions and against brute force.

use k114::curve::{self, Stats};
use k114::exact::{self, Triple};
use k114::filters::Filters;
use k114::roots;

fn gcd_u64(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

fn scan_d(
    k: i64,
    d: u64,
    factors: &[(u64, u32)],
    rlo: f64,
    rhi: f64,
    hits: &mut Vec<Triple>,
    st: &mut Stats,
) {
    let mut facs: Vec<(u64, u32)> = factors.to_vec();
    for f in facs.iter_mut() {
        if f.0 == 2 {
            f.1 = f.1.saturating_sub(2);
        }
    }
    facs.retain(|x| x.1 > 0);
    let Some((modulus, rs)) = roots::roots_from_factorization(k as u64, &facs) else {
        return; // no cube root of k modulo d' — nothing to scan
    };
    assert_eq!(modulus, d / gcd_u64(d, 4), "d' mismatch for d={}", d);
    let f = Filters::new(k);
    // The alpha lemma needs |z| > sqrt(k).
    let z_min = ((k as f64).sqrt() as u128) + 1;
    let zlo = (((d as f64) * rlo).floor() as u128).max(z_min);
    let zhi = ((d as f64) * rhi).ceil() as u128;
    for r in rs {
        curve::scan_all(&f, k, d, r, modulus, zlo, zhi, st, hits);
    }
}

#[test]
fn published_solutions_are_recovered() {
    // (k, d, factorisation of d, |z|)
    let cases: [(i64, u64, Vec<(u64, u32)>, u128); 4] = [
        (3, 108_398_887_211, vec![(167, 1), (649_095_133, 1)], 472_715_493_453_327_032),
        (33, 87_723_532_425_289, vec![(87_723_532_425_289, 1)], 2_736_111_468_807_040),
        (
            42,
            102_980_666_258_459,
            vec![(11, 1), (43, 1), (215_921, 1), (1_008_323, 1)],
            12_602_123_297_335_631,
        ),
        (
            795,
            21_083_965_616_656,
            vec![(2, 4), (17, 1), (149, 1), (520_232_077, 1)],
            2_337_348_783_323_923,
        ),
    ];

    for (k, d, factors, az) in cases {
        let mut hits = Vec::new();
        let mut st = Stats::default();
        let r = az as f64 / d as f64;
        scan_d(k, d, &factors, r * 0.999, r * 1.001, &mut hits, &mut st);
        assert!(st.balanced(), "accounting broken for k={}", k);
        assert!(
            hits.iter().any(|t| t.z.unsigned_abs() == az && exact::verify(k, t.x, t.y, t.z)),
            "did not recover the published solution for k={}",
            k
        );
    }
}

/// Every brute-force solution inside the scanned region must be found.
#[test]
fn brute_force_agrees_on_small_k() {
    let rlo = 3.4f64;
    let rhi = 64.0f64;
    let dmax = 240u64;
    let boxn = 90i128;

    for k in [2i64, 6, 10, 43, 71, 100, 128] {
        let ki = k as i128;
        // Brute force: solutions with |z| minimal and R in range.
        let mut expected: Vec<(i128, i128, i128)> = Vec::new();
        for x in -boxn..=boxn {
            for y in -boxn..=boxn {
                for z in -boxn..=boxn {
                    if x * x * x + y * y * y + z * z * z != ki {
                        continue;
                    }
                    if z == 0 || z * z <= ki {
                        continue;
                    }
                    if z.abs() > x.abs() || z.abs() > y.abs() {
                        continue;
                    }
                    let d = (x + y).unsigned_abs();
                    if d == 0 || d > dmax as u128 {
                        continue;
                    }
                    let r = z.abs() as f64 / d as f64;
                    if r < rlo || r > rhi {
                        continue;
                    }
                    expected.push((x, y, z));
                }
            }
        }
        if expected.is_empty() {
            continue;
        }
        // Scan the whole d range and every z the region allows.
        let mut hits: Vec<Triple> = Vec::new();
        let mut st = Stats::default();
        for d in 1..=dmax {
            // Brute-force factorisation is fine at this size.
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
            scan_d(k, d, &facs, rlo, rhi, &mut hits, &mut st);
        }
        for (x, y, z) in expected {
            // x and y are interchangeable: the parameterisation returns them in
            // the order forced by (s + v)/2, so compare the unordered pair.
            let mut want = [x, y];
            want.sort();
            assert!(
                hits.iter().any(|t| {
                    let mut got = [t.x, t.y];
                    got.sort();
                    got == want && t.z == z
                }),
                "missed solution ({}, {}, {}) for k={}",
                x,
                y,
                z,
                k
            );
        }
        for t in hits.iter() {
            assert!(exact::verify(k, t.x, t.y, t.z));
        }
    }
}

/// The exact test must never accept a non-solution, at any size.
#[test]
fn no_false_positives_in_a_small_window() {
    let k = 114i64;
    let d = 1_000_003u64; // prime, to exercise the generic root path
    let facs = vec![(1_000_003u64, 1u32)];
    let mut hits = Vec::new();
    let mut st = Stats::default();
    scan_d(k, d, &facs, 3.4, 3.5, &mut hits, &mut st);
    assert!(st.balanced());
    for t in hits {
        assert!(exact::verify(k, t.x, t.y, t.z), "false positive: {}", t);
    }
}
