//! CLI: scan an annulus `d in [dlo, dhi)`, `R = z/d in [rlo, rhi]` for
//! `x^3 + y^3 + z^3 = k`, writing hits and a region manifest.

use std::fs::File;
use std::io::Write;
use std::time::Instant;

use k114::curve::{self, Stats};
use k114::dgen::{self, DEntry, PrimeBits, MAX_FACTORS, NFAC_OVERFLOW};
use k114::exact::Triple;
use k114::filters::Filters;
use k114::roots;

struct Cfg {
    k: i64,
    dlo: u64,
    dhi: u64,
    rlo: f64,
    rhi: f64,
    threads: usize,
    block: u64,
    out: String,
}

fn usage() {
    eprintln!(
        "usage: k114 --dlo <u64|1.2e18> --dhi <...> [--k 114] [--rlo 3.4] [--rhi 64]\n\
         \x20            [--threads N] [--block N] [--out DIR]\n\n\
         Region: d in [dlo, dhi), and |z| in [rlo*d, rhi*d] with R = |z|/d.\n\
         The lower edge rlo must be below 1/(cuberoot(2)-1) = 3.8473 only for\n\
         degenerate ties; 3.4 is a safe default."
    );
}

fn parse_u64_loose(s: &str) -> u64 {
    if let Ok(v) = s.parse::<u64>() {
        return v;
    }
    let f: f64 = s.parse().unwrap_or_else(|_| panic!("bad integer: {}", s));
    assert!(f >= 0.0 && f < 1.9e19, "value out of range: {}", s);
    f as u64
}

fn gcd_u64(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

fn process(f: &Filters, k: i64, rlo: f64, rhi: f64, entries: &[DEntry]) -> (Stats, Vec<Triple>) {
    let mut st = Stats::default();
    let mut hits = Vec::new();
    let mut factors: Vec<(u64, u32)> = Vec::with_capacity(MAX_FACTORS);

    for e in entries {
        if e.nfac == NFAC_OVERFLOW {
            st.dense_skipped += 1;
            continue;
        }
        // For k ≡ 3,6 (mod 9): x^3+y^3 ≡ 0 (mod 9) whenever x+y ≡ 0 (mod 3),
        // so no solution can have 3 | d. Skip a third of the enumeration.
        if k.rem_euclid(3) == 0 && e.d % 3 == 0 {
            st.div3_skipped += 1;
            continue;
        }
        factors.clear();
        for i in 0..e.nfac as usize {
            factors.push((e.p[i], e.e[i] as u32));
        }
        // The congruence only needs d' = d / gcd(d, 4); the 2-adic exponent drops.
        let g = gcd_u64(e.d, 4);
        for fac in factors.iter_mut() {
            if fac.0 == 2 {
                fac.1 = fac.1.saturating_sub(2);
            }
        }
        factors.retain(|x| x.1 > 0);

        let Some((modulus, roots_of_k)) = roots::roots_from_factorization(k as u64, &factors) else {
            st.no_roots += 1;
            continue;
        };
        debug_assert_eq!(modulus, e.d / g, "modulus must be d'");

        let d = e.d;
        let zlo_f = (d as f64) * rlo;
        let zhi_f = (d as f64) * rhi;
        // The alpha lemma needs |z| > sqrt(k); below that the equation is a
        // different (Thue-type) problem.
        let z_min = dgen::isqrt_u64(k.unsigned_abs()) as u128 + 1;
        let zlo = if zlo_f < z_min as f64 {
            z_min
        } else {
            zlo_f.floor() as u128
        };
        let zhi = zhi_f.ceil() as u128;
        if zhi <= zlo {
            continue;
        }
        for &r in roots_of_k.iter() {
            curve::scan_all(f, k, d, r, modulus, zlo, zhi, &mut st, &mut hits);
        }
    }
    (st, hits)
}

fn main() {
    let mut cfg = Cfg {
        k: 114,
        dlo: 0,
        dhi: 0,
        rlo: 3.4,
        rhi: 64.0,
        threads: std::thread::available_parallelism()
            .map(|n| n.get().min(8))
            .unwrap_or(4),
        block: 1 << 14,
        out: String::new(),
    };

    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < args.len() {
        let key = args[i].as_str();
        let val = args.get(i + 1).cloned().unwrap_or_default();
        match key {
            "--k" => cfg.k = val.parse().expect("--k"),
            "--dlo" => cfg.dlo = parse_u64_loose(&val),
            "--dhi" => cfg.dhi = parse_u64_loose(&val),
            "--rlo" => cfg.rlo = val.parse().expect("--rlo"),
            "--rhi" => cfg.rhi = val.parse().expect("--rhi"),
            "--threads" => cfg.threads = val.parse().expect("--threads"),
            "--block" => cfg.block = val.parse().expect("--block"),
            "--out" => cfg.out = val,
            "-h" | "--help" => {
                usage();
                return;
            }
            other => {
                eprintln!("unknown argument: {}", other);
                usage();
                std::process::exit(2);
            }
        }
        i += 2;
    }

    if cfg.dlo == 0 || cfg.dhi <= cfg.dlo {
        eprintln!("need 0 < dlo < dhi");
        usage();
        std::process::exit(2);
    }
    if !(1.0..=1000.0).contains(&cfg.rlo) || cfg.rhi <= cfg.rlo {
        eprintln!("need 1 <= rlo < rhi");
        std::process::exit(2);
    }

    let sieve_limit = dgen::isqrt_u64(cfg.dhi - 1) + 1;
    eprintln!(
        "sieve: primes up to {} ({:.1} MB bitmap)",
        sieve_limit,
        (sieve_limit / 16) as f64 / 1e6
    );
    let t0 = Instant::now();
    let sieve = PrimeBits::new(sieve_limit);
    eprintln!("sieve built in {:.1}s", t0.elapsed().as_secs_f64());

    let filters = Filters::new(cfg.k);
    let out_dir = if cfg.out.is_empty() {
        format!("runs/k{}", cfg.k)
    } else {
        cfg.out.clone()
    };
    std::fs::create_dir_all(&out_dir).expect("create output directory");
    let mut hits_file = File::create(format!("{}/hits.jsonl", out_dir)).expect("create hits.jsonl");
    let mut progress_file =
        File::create(format!("{}/progress.log", out_dir)).expect("create progress.log");

    eprintln!(
        "region: k={} d in [{}, {}) R in [{}, {}] threads={} block={}",
        cfg.k, cfg.dlo, cfg.dhi, cfg.rlo, cfg.rhi, cfg.threads, cfg.block
    );

    let start = Instant::now();
    let mut total = Stats::default();
    let mut blocks_done = 0u64;
    let mut block_lo = cfg.dlo;

    while block_lo < cfg.dhi {
        let block_hi = (block_lo + cfg.block).min(cfg.dhi);
        let entries = dgen::factor_block(block_lo, block_hi, &sieve);

        let nthreads = cfg.threads.max(1);
        let chunk = (entries.len() + nthreads - 1) / nthreads;
        let results: Vec<(Stats, Vec<Triple>)> = std::thread::scope(|s| {
            let mut handles = Vec::new();
            for part in entries.chunks(chunk.max(1)) {
                let f = &filters;
                let k = cfg.k;
                let rlo = cfg.rlo;
                let rhi = cfg.rhi;
                handles.push(s.spawn(move || process(f, k, rlo, rhi, part)));
            }
            handles
                .into_iter()
                .map(|h| h.join().expect("worker thread panicked"))
                .collect()
        });

        for (st, hs) in results {
            total.merge(&st);
            for t in hs.iter() {
                let d = (t.x + t.y).unsigned_abs();
                writeln!(
                    hits_file,
                    "{{\"k\":{},\"x\":\"{}\",\"y\":\"{}\",\"z\":\"{}\",\"d\":{},\"verified\":true}}",
                    cfg.k, t.x, t.y, t.z, d
                )
                .expect("write hit");
                eprintln!("HIT k={} (x,y,z)=({}, {}, {})", cfg.k, t.x, t.y, t.z);
            }
            if !hs.is_empty() {
                hits_file.flush().ok();
            }
        }

        blocks_done += 1;
        block_lo = block_hi;

        let el = start.elapsed().as_secs_f64();
        let pts = total.points();
        writeln!(
            progress_file,
            "d<={} blocks={} points={} exact={} hits={} secs={:.1}",
            block_lo, blocks_done, pts, total.exact_tests, total.hits, el
        )
        .ok();
        if blocks_done % 8 == 0 || block_lo >= cfg.dhi {
            eprintln!(
                "d<={:.4e} blocks={} points={:.3e} ({:.2e}/s) exact={} curves={} hits={} {:.0}s",
                block_lo as f64,
                blocks_done,
                pts as f64,
                pts as f64 / el.max(1e-9),
                total.exact_tests,
                total.curves,
                total.hits,
                el
            );
        }
    }

    if !total.balanced() {
        eprintln!("ERROR: point accounting is inconsistent: {:?}", total);
        std::process::exit(1);
    }

    let el = start.elapsed().as_secs_f64();
    let manifest = format!(
        "{{\n  \"k\": {},\n  \"d_lo\": {},\n  \"d_hi\": {},\n  \"r_lo\": {},\n  \"r_hi\": {},\n  \
         \"engine\": \"k114 {}\",\n  \"threads\": {},\n  \"block\": {},\n  \"points\": {},\n  \
         \"curves\": {},\n  \"exact_tests\": {},\n  \"hits\": {},\n  \"no_roots\": {},\n  \
         \"dense_skipped\": {},\n  \"seconds\": {:.1},\n  \"points_per_second\": {:.3e},\n  \
         \"completed\": {}\n}}\n",
        cfg.k,
        cfg.dlo,
        cfg.dhi,
        cfg.rlo,
        cfg.rhi,
        env!("CARGO_PKG_VERSION"),
        cfg.threads,
        cfg.block,
        total.points(),
        total.curves,
        total.exact_tests,
        total.hits,
        total.no_roots,
        total.dense_skipped,
        el,
        total.points() as f64 / el.max(1e-9),
        if blocks_done > 0 && block_lo >= cfg.dhi { "true" } else { "false" }
    );
    std::fs::write(format!("{}/region.json", out_dir), manifest).expect("write region.json");

    eprintln!(
        "done: points={:.3e} in {:.1}s ({:.2e}/s), exact={}, hits={}, output={}",
        total.points() as f64,
        el,
        total.points() as f64 / el.max(1e-9),
        total.exact_tests,
        total.hits,
        out_dir
    );
}
