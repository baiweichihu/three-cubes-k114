# DESIGN — a local solver for $x^3+y^3+z^3=114$

> Developer-facing document: mathematics, algorithm selection, coverage survey, architecture, and verification plan.
> Status: **design frozen, implementation started**. Date: 2026-11-23.

---

## 0. One paragraph

Fold $x,y$ away with $d=x+y$. The equation becomes two-dimensional, $v^2=(4(k-z^3)-d^3)/(3d)$, where $v=x-y$ must make the right side a perfect square. The congruence $z^3\equiv k \pmod d$ collapses the admissible $z$ into arithmetic progressions; a mod-243 residue wheel, a parity test, and 15 small-prime quadratic-residue bitmaps reject candidates at a measured $1:5\times10^5$ rate; only the survivors get a 256-bit exact square-root test. Pure integer arithmetic, no network, no crowd-sourcing.

---

## 1. Goals / non-goals

### Goals

1. A correct, self-certifying, resumable $k=114$ kernel; speed target ≥ 2× the reference C implementation `zcubes`.
2. Generalise over $k$; regress against the known solutions for $k=3,33,42,795$.
3. Be able to state "*this $(d,z)$ rectangle has been scanned completely*" and leave a re-checkable manifest.
4. Serve as genuine Rust performance engineering: fixed-width big integers, bitmap sieving, cache blocking, SIMD.

### Non-goals

- No crowd-sourcing, no issue submission, no "banks", no network.
- No shared code or data with `../k=114` (Math Gambling console). That project *samples* the space and submits receipts; this one *scans contiguous regions* locally.
- **"Finding 114" is not the success criterion.** See §3.6: a laptop buys roughly a $10^{-5}$ hit probability per year. The deliverable is the engine.
- No asymptotic improvement is attempted — the literature is at $O(N^{1+o(1)})$ and there is nothing better (§2.7).

---

## 2. Mathematics

### 2.1 The equation and the local obstruction

$$x^3+y^3+z^3=k,\qquad k=114$$

- Cubes are $\{0,\pm1\}\bmod 9$, so $k\equiv\pm4\pmod 9$ is impossible. $114=12\cdot9+6\equiv6$: **no local obstruction**.
- Open for $k\le1000$: 114, 390, 627, 633, 732, 921, 975 (all $\equiv3,6\bmod9$). 114 is the smallest.
- Known large solutions (all verified by direct computation):

| $k$ | solution | $\max\lvert x\rvert$ |
|---|---|---|
| 33 | $(8866128975287528,\,-8778405442862239,\,-2736111468807040)$ | $8.9\times10^{15}$ |
| 42 | $(-80538738812075974,\,80435758145817515,\,12602123297335631)$ | $8.1\times10^{16}$ |
| 3 | $(569936821221962380720,\,-569936821113563493509,\,-472715493453327032)$ | $5.7\times10^{20}$ |
| 795 | $(-14219049725358227,\,14197965759741571,\,2337348783323923)$ | $1.4\times10^{16}$ |

### 2.2 Folding two variables away

Let $s=x+y$ and $v=x-y$. Then $x^3+y^3=s(s^2+3v^2)/4$, so

$$4(k-z^3)=s^3+3sv^2\quad\Longrightarrow\quad \boxed{\,v^2=\frac{4(k-z^3)-s^3}{3s}\,},\qquad x,y=\frac{s\pm v}{2}$$

Notation used throughout:

- $d:=|s|=|x+y|>0$, with $s=\pm d$. The sign of $s$ is forced: $s(s^2+3v^2)/4=k-z^3$ and $s^2+3v^2>0$, so $\operatorname{sign}(s)=\operatorname{sign}(k-z^3)$.
- The three variables may be relabelled freely; **we normalise $z$ to be the smallest in absolute value** ($|z|\le\min(|x|,|y|)$). This is what makes the cost proportional to $\min$ rather than $\max$.

### 2.3 The exact acceptance test

Given $(s,z)$:

1. $t:=k-z^3$, require $3s \mid 4t-s^3$;
2. $v^2:=(4t-s^3)/(3s)$, require $v^2\ge0$;
3. $v:=\lfloor\sqrt{v^2}\rfloor$, require $v\cdot v=v^2$;
4. require $s\equiv v\pmod 2$ (else $x,y\notin\mathbb Z$);
5. output $(x,y,z)=((s+v)/2,(s-v)/2,z)$, then **independently re-verify $x^3+y^3+z^3=k$ in exact arithmetic**.

All integer arithmetic; no floating point anywhere in the verdict path.

### 2.4 What the divisibility condition means (this is where the sieve enters)

Requiring $3d\mid 4t-d^3$ splits into two conditions:

- **(C1)** $d\mid 4t$. With $d':=d/\gcd(d,4)$ this is
  $$z^3\equiv k\pmod{d'}$$
- **(C2)** with $A:=4t/s\in\mathbb Z$ (sign of $s$ included), require
  $$A\equiv s^2 \pmod 3$$

(C2) is one cheap mod-3 test. In the "enumerate $d$" framework this *is* the whole of the literature's "3-adic lifting" — we never factor $k-z^3$, so nothing deeper is needed.

Consequences of (C1). Two of them are strong enough to remove whole thirds of the search space, and both were discovered while implementing the engine (the tests in §6 pin them down):

- **$3\nmid d$ whenever $k\equiv3,6\pmod 9$.** If $3\mid s=x+y$ then $y\equiv-x\pmod3$, so $x^2-xy+y^2\equiv3x^2\equiv0\pmod3$ and hence $x^3+y^3=s(x^2-xy+y^2)\equiv0\pmod9$. Then $z^3\equiv k\pmod 9$ with $k\equiv3$ or $6$, but cubes mod 9 are $0,\pm1$ — contradiction. So a third of the $d$ range is empty by pure arithmetic (the engine skips it and counts the skips).
- **For $k\equiv6\pmod 9$: $s\equiv1\pmod3$ and $z\equiv2\pmod3$.** Since $x+y+z\equiv k\equiv0\pmod3$ and $s\not\equiv0$, only $s\equiv1,2$ remain; if $s\equiv2$ then $z\equiv1$ and $z^3\equiv1\pmod9$, forcing $x^3+y^3\equiv5\pmod9$, which is impossible. The mod-243 wheel encodes exactly this: `allowed[s]` is empty for every $s\not\equiv1\pmod3$ when $k=114$.
- $9\mid d'$ is likewise unsatisfiable, hence $9\nmid d'$ — subsumed by $3\nmid d$ above for the $k$ we care about.
- For primes $p\equiv1\pmod3$ dividing $d'$, $k$ must be a cubic residue mod $p$; if not, skip. For $p\equiv2\pmod3$ the cube map is a bijection, so exactly one root always exists.
- $p=3$: the cube map is the identity on $\mathbb F_3$; when $3\mid k$ the root is $z\equiv0\pmod3$.
- $p=2$: the cube map is the identity on $\mathbb F_2$.
- If $p\mid k$, write $v=v_p(k)$ and $z=p^wu$. Then $v_p(z^3-k)$ equals $3w$ (if $3w<v$), $v$ (if $3w>v$), or $v+v_p(u^3-k/p^v)$. For $3\nmid v$ that requires $v\ge e$ and $w\ge\lceil(v+1)/3\rceil$; for $3\mid v$ the engine enumerates directly.

### 2.5 Curves and quotient points

Fix $d$ and a root $r$ of $z^3\equiv k\pmod{d'}$. All admissible $z$ form the arithmetic progression

$$z=r+d\,q,\qquad q\in\mathbb Z$$

- a **curve** is a pair $(d,r)$ (with a sign for $s$);
- a **quotient point** is one value of $q$.

One curve costs one factorisation + one cube-root computation and then walks an arithmetic progression; the $q$ range comes from the target $z$ window. Cost per $d$ is therefore $\approx (R_{\max}-R_{\min})$ points, where $R:=z/d$.

### 2.6 Height relation, the $\alpha$ lemma, and hit density

For $|z|^3\gg k$ we have $v^2\approx-4z^3/(3s)$, hence

$$\boxed{\,H:=\max(|x|,|y|)\approx\sqrt{\frac{z^3}{3d}}=d\,\frac{R^{3/2}}{\sqrt3}\,}$$

**Lemma (Booker–Sutherland).** If $|z|>\sqrt k$ then $0< d< \alpha|z|$ with

$$\alpha=\sqrt[3]{2}-1=0.2599210498948732,\qquad 1/\alpha=3.8473221018630723$$

i.e. **$R=z/d>3.847332$ always**. Verification of the constant: the extremal configuration has the two smallest coordinates equal, $|y|=|z|=M$, whence $x^3=2M^3$ and $d=|x+y|=(\sqrt[3]2-1)M=\alpha|z|$; e.g. $43=9^3-7^3-7^3$ gives $d=2$, $z=7$, $R=3.5$ — inside the bound, and excluded only because the ordering is not strict. Ties are Thue-type degenerate cases; a scan with $R\ge3.4$ covers them anyway at negligible cost.

Validation of the $H$ formula on the four known solutions:

| $k$ | $d=\lvert x+y\rvert$ | $z$ (the minimum) | $d\,R^{3/2}/\sqrt3$ | actual $\max\lvert x\rvert$ | error | $R=z/d$ |
|---|---|---|---|---|---|---|
| 33 | 87 723 532 425 289 | $-2.7361\times10^{15}$ | $8.822\times10^{15}$ | $8.866\times10^{15}$ | 0.5% | 31.19 |
| 42 | 102 980 666 258 459 | $1.2602\times10^{16}$ | $8.0487\times10^{16}$ | $8.0539\times10^{16}$ | 0.06% | 122.4 |
| 3 | 108 398 887 211 | $-4.7272\times10^{17}$ | $5.699368\times10^{20}$ | $5.699368\times10^{20}$ | $<10^{-6}$ | $4.36\times10^{6}$ |
| 795 | 21 083 965 616 656 | $2.3373\times10^{15}$ | $1.42085\times10^{16}$ | $1.42191\times10^{16}$ | 0.07% | 110.9 |

**Hit density.** For fixed $d$, $v^2$ is an essentially random integer of size $\approx 4|z|^3/(3d)$, so it is a perfect square with probability $\approx1/(2v)$:

$$\boxed{\,P\approx\frac{\sqrt{3d}}{4\,z^{3/2}}\,}\qquad\text{density per point}=\frac{P}{d}=\frac{\sqrt3}{4\sqrt d\,z^{3/2}}$$

Consequence: **large $d$ and small $z$ are worth more per point** (subject to $d\le\alpha z$). "Bigger is harder" is false — what grows is the volume to be covered, not the sparsity per point.

**Total expectation** (Heath-Brown heuristic; constants from Booker–Sutherland Table 1):

$$\mathbb E[\#\text{solutions},\ \text{height}\le B]=\rho(k)\ln B,\qquad \rho(114)=0.058459$$

So $e(B)=1$ only at $B=e^{1/\rho}\approx2.685\times10^{7}$, and:

| $B$ | $10^{5}$ | $10^{10}$ | $10^{15}$ | $10^{17}$ |
|---|---|---|---|---|
| $\mathbb E$ | 0.673 | 1.346 | 2.019 | 2.29 |

Every further decade adds $0.058$ — i.e. **5.9% per order of magnitude**. `k=114` has survived to $10^{17}$ with 0 found; Poisson gives a 13% probability of that. Heath-Brown's conjecture is stated to hold for any reasonable height function (min, $|xyz|^{1/3}$, min of pairwise sums), so switching height does not change the asymptotics.

### 2.7 Why this parameterisation and not another

| Route | Mechanism | Strengths | Weaknesses | Verdict |
|---|---|---|---|---|
| **A. Enumerate $d$ directly** (chosen) | enumerate $d$ → compute $\sqrt[3]{k}\bmod d$ → walk $q$ | **Contiguous** coverage; "this rectangle is done" is provable | needs factorisation of $d$ + cube roots | **adopt** |
| **B. Mordell norm form** (the `mg114-offset-*` scheme in `../k=114`) | small $(b,c)$, $d=N(a,b,c)/\ell$ with $N=a^3+kb^3+k^2c^3-3kabc$, root comes for free algebraically | extremely cheap per point, jumps straight to $|z|\sim10^{22}$ | only *samples*: the reachable $d$ are the image of $(b,c,t)$, not an interval | keep as a cross-check, never the main line |
| **C. Elkies lattice reduction** (Elsenhans–Jahnel 2009) | LLL short vectors | effective when $|z|$ is known to be small | worse for large ranges, superseded after 2009 | reject |

Published complexity ladder (Sutherland's lecture notes): naive $O(N^3)$ → naive-improved $O(N^{2+o(1)})$ → **(A) in $O(N^{1+o(1)})$**. Nothing better is published.

---

## 3. Survey of published search ranges, and the region we take

### 3.1 Every published brute-force range found

All rows are *brute-force style* searches of the same equation (analytic/lattice papers without a range are excluded).

| Search | Year | Height definition | Range | Method / notes |
|---|---|---|---|---|
| Beck–Pine–Tarrant–Jensen (BPTJ) | 2007 | $\max$ | bound not restated in later surveys | historical survey; 3 values of $k<100$ and 27 of $k<1000$ still open |
| Elsenhans–Jahnel | 2009 | $\max$ | $\le10^{14}$ | Elkies' lattice reduction (LLL); 14 values $<1000$ open |
| Huisman | 2016 | $\max$ | $\le10^{15}$ | C + sieve; 966 new solutions, 15 254 total; first solution for $k=74$; 13 values $<1000$ open |
| Booker | 2019 | $\min$ | $\le10^{16}$ | introduced the $(d,z)$ parameterisation; found 33 and 795; 3145 core-days |
| Booker–Sutherland, round 1 | Sep 2019 | $\min$ | $z\le10^{17}$, $d\le\alpha z$ | Charity Engine, ~500k PCs; **exhaustive** for the 11 then-open $k$; found 42, 165, 906 |
| Booker–Sutherland, round 2 | 2020 | $\min$ | $z\le10^{18}$, $d\le\alpha z/9$ | $k=3$ only; found the $5.7\times10^{20}$ solution |
| Booker–Sutherland, round 3 | 2020–21 | $\min$, **partial** | $z\le10^{19}$, $d\le z/54=1.8519\times10^{17}$ | *not* a full $\min\le10^{19}$ search; found 579 (Jan 2021); **includes $k=114$** |
| Booker–Sutherland, $k=3$ verification | 2020 | $\min$ | $z\le4.7272\times10^{17}$ | Google Compute Engine, 155 579 nodes, $<4$ h, ≈120 core-years; 5 faulty runs, redone |
| Booker–Sutherland, Remark 3.8 | 2021 | $d$ | $d\le100$ unconditional; $d\le20\,000$ under GRH | 3-isogeny descent (Bremner) — an *exclusion* by $d$, not a scan |
| Booker–Sutherland, $k=33$ benchmark | 2021 | $\min$ | $z\le2^{53}$, $d\le2^{47}$ ($R=64$) | 107 core-days; vs Booker's 3145 → 25× faster |
| Math Gambling `k=114` campaign | 2025–26 | $(d,z)$ sampling | $d\in(1.8519\times10^{17},\,1.4815\times10^{18}]$; $z/d\in[3.85,\,4096]$, bands $R\le64,256,4096$ | volunteer campaign, **sparse random sampling** (~$7.7\times10^5$ verified tasks out of a $10^{12}$-scale task space); targets 114 + 6 other open $k$ |
| Status trackers (TheoremDB P2650, ProblemHive, Epoch AI) | 2026 | — | record: exhaustive frontier is $10^{17}$ in the smallest coordinate; a later *partial* run reached $|z|=10^{19}$; the $\max\le10^{20}$ box is **not** certified | independently corroborates the rows above |

Two things worth flagging, because they are easy to get wrong:

1. **`min ≤ 1e19` was never exhausted.** Round 3 capped $d$ at $z/54$, while the geometry allows $d\le\alpha z$ — 14× more $d$. Everything with $1.8519\times10^{17}<d\le\alpha z$ at $z\le10^{19}$ was untouched by round 3.
2. **The Math Gambling console inherits that exact boundary**: its `D0 = 10**19 // 54 = 185185185185185185` *is* round 3's $d_{\max}$. Its $d$ shells are $D_0\cdot2^{j}$, $j=0,1,2$, so it covers $d\in(1.8519\times10^{17},1.4815\times10^{18}]$ — and nothing above.

### 3.2 Coverage map in the $(d,z)$ plane

With $z$ normalised to the minimum and $R=z/d>1/\alpha=3.8473$:

| # | Covered set | Source |
|---|---|---|
| C1 | $\{z\le10^{17}\}$ for all $d$ | round 1 (exhaustive) |
| C2 | $\{z\le10^{19}\}\cap\{d\le1.8519\times10^{17}\}$ | round 3 (partial) |
| C3 | $\{d\le100\}$ (GRH: $d\le20\,000$) | Remark 3.8 |
| C4 | $\{1.8519\times10^{17}<d\le1.4815\times10^{18}\}\cap\{3.85\le R\le4096\}$ | Math Gambling (sampled, so treat as *claimed*) |
| C5 | $\{\max\le10^{15}\}$ | Huisman / EJ (subsumed by C1) |

Derived quantities: round 3's effective height floor is $H\ge d\sqrt{(1/\alpha)^3/3}=4.356\,d$; with $d>D_0=1.8519\times10^{17}$ that is $H\gtrsim8.1\times10^{17}$. Math Gambling's ceiling is $z\le4096d\le6.07\times10^{21}$, i.e. $H\lesssim2.2\times10^{23}$.

### 3.3 Virgin regions

| Region | Definition | $H$ range | Status |
|---|---|---|---|
| **V1 (chosen)** | $d>1.4815\times10^{18}$, $3.4\le R\le64$ | $\ge5.4\times10^{18}$ | disjoint from C1–C5 |
| V2 | $d\le1.8519\times10^{17}$, $z>10^{19}$ | $\ge4.2\times10^{19}$ | disjoint, but lower density |
| V3 | $z>4096d$ inside C4's $d$ shells | $\ge2.2\times10^{23}$ | disjoint, very expensive |
| V4 | $1.4815\times10^{18}<d\le\alpha z$, $R>64$ | larger | disjoint, thin bands |

### 3.4 Why V1

**Disjointness proof.** For $d>1.4815\times10^{18}$:

- $d>1.8519\times10^{17}$ → outside C2 and C4;
- $d>20\,000$ → outside C3 (even under GRH);
- $R\le64$ and $R\ge3.4$ ⟹ $z\ge3.4d>5.0\times10^{18}>10^{17}$ → outside C1;
- $H=dR^{3/2}/\sqrt3\ge5.4\times10^{18}>10^{15}$ → outside C5.

**Density advantage.** At its cheapest corner ($R=3.4$, $d=1.4815\times10^{18}$): $P=\sqrt{3d}/(4z^{3/2})\approx4.7\times10^{-20}$, versus $5.9\times10^{-21}$ for V2's cheapest corner ($d=1.8519\times10^{17}$, $z=10^{19}$) — about 8× better per point. V1 also continues the *same axis* (d) that the Math Gambling campaign stops on, which makes the "unscanned" claim crisp and easy to state.

**Aspect ratio.** Booker–Sutherland's cost model (their Remark 2.2) is runtime $\propto\rho_{ap}d_{\max}$, and their tuning (Table 3, $k=33$ benchmark) finds the optimum near $R=z_{\max}/d_{\max}=64$ once $d_{\max}\ge2^{50}$. So $R\in[3.4,64]$ is the cost-optimal band; Math Gambling's $R\le64$ band agrees. Cost $\approx(R_{\max}-R_{\min})\cdot\Delta d\approx60.6\,\Delta d$ points.

### 3.5 Concrete starting block

$$k=114,\quad d\in(1.4815\times10^{18},\,1.4815\times10^{18}+\Delta d],\quad R=\frac zd\in[3.4,\,64]$$

with $\Delta d$ set from a measured throughput (see §3.6), starting at the smallest executable block and growing.

### 3.6 What a laptop can actually buy — measured

Expected successes in the region scale like $\mathbb E\approx0.042\,\ln\!\big((d_0+\Delta d)/d_0\big)$ (obtained by integrating the point density and calibrating against $\rho(114)$).

The cost driver turned out to be **per-$d$ setup** (factorisation + cube roots + wheel split), not the $q$ loop. First measurement, `k114` v0.1.0, `--threads 8`, $d\in[10^{12},1.0021\times10^{12})$ (2.1e9 values, 743 s):

| Quantity | Measured |
|---|---|
| $d$ processed per second | $2.8\times10^{6}$ |
| effective points (post-wheel) per second | $6.8\times10^{6}$ |
| curves scanned | $1.8\times10^{9}$ (~0.87 per $d$) |
| $d$ with no cube root mod $d'$ | 45% |
| exact tests | $2.0\times10^{5}$ (1 in $2.5\times10^4$ points) |
| hits | 0 (expected: none in this range) |

Same measurement repeated at the real target scale, $d\approx1.4815\times10^{18}$ (§7.3), gives $3.8\times10^{5}$ $d$/s — 7× slower, because exact segmented factorisation has to visit every prime up to $1.2\times10^{9}$ once per block.

So one year of wall-clock buys

| where | $\Delta d$ per year | $\mathbb E$ |
|---|---|---|
| $d\sim10^{12}$ (fast regime) | $8.8\times10^{13}$ | $2.5\times10^{-6}$ (1 in 4e5) |
| $d\sim1.5\times10^{18}$ (target, as built) | $1.2\times10^{13}$ | $3.4\times10^{-7}$ (1 in 3e6) |
| target, after the M5 enumeration rewrite (assume 7×) | $8.4\times10^{13}$ | $2.4\times10^{-6}$ |

A 10× optimisation would make that $2.5\times10^{-5}$; even a 1000× one leaves it a lottery. Read this as: **region choice buys honesty and non-duplication, not odds.** The odds are set by the heuristic and are brutal for any single machine. The engine, the verification discipline, and the coverage bookkeeping are the actual products.

> Cost model, corrected by measurement: the "~60 q points per $d$" of the design phase is the *pre-wheel* figure. The mod-243 wheel leaves only $|\text{allowed}(s)|/243$ of them, about $1/20$ for $k=114$, so ~2.4 effective points per $d$ — which is why setup cost dominates and why root caching is the first optimisation to try.

> Optional cheap win: $\rho(k)$ varies between the seven open $k$ (Booker–Sutherland Table 1). Picking the largest $\rho$ costs nothing architecturally and multiplies the odds. Keep $k$ configurable.

---

## 4. Algorithm

### 4.1 Pipeline

```
outer: d-block (Δd), in order across region V1
   ├─ 1. enumerate d in the block with full factorisation
   ├─ 2. cube roots r ≡ ∛k (mod d')            (multiple roots possible)
   ├─ 3. mod-3 part of (C2)                    (absorbed into the 243 wheel)
   ├─ 4. q expansion: z window → q range, split into APs of step lcm(d,243)/d
   └─ 5. filter chain (hot path)
          ├─ mod-243 wheel
          ├─ parity: k − s − z ≡ 0 (mod 2)
          ├─ 15 small-prime quadratic-residue bitmaps
          └─ (optional) composite wheel
        6. exact test (cold path, ≈1:5e5) → u256 isqrt → re-verify → record
```

### 4.2 Step 1: enumerating $d$ with factorisation

For a *contiguous* interval the only correct way is **segmented factorisation**: sieve the block with all primes $\le\sqrt{d_{\max}}$, dividing out each prime power; whatever remains $>1$ is a prime factor.

- Base sieve: odd-only bitset up to $\sqrt{d_{\max}}$ ($1.22\times10^{9}$ for $d\sim1.5\times10^{18}$ → 76 MB, ~10 s once).
- Per block of $\Delta d$: for each prime $p$ (streamed from the base bitset), step through its multiples in the block and divide out $p^{e}$; record $(p,e)$.
- Cost: $\Delta d\cdot\log\log d$ division steps plus one "first multiple" per prime per block. Measured estimate: ~50–200 ns per $d$ — comparable to the $q$ loop, and the reason the reference implementation exposes `pmin/pmax` (restricting the largest prime factor makes the root computation cheap and table-able; **hypothesis to confirm by reading `zcubes`'s `cbrts.h`**).

Alternative rejected for the main line: DFS over prime factorisations with a smoothness bound (that is effectively what `pmin/pmax` does) — it cannot cover $d$ that are prime near $10^{18}$, and note $k=33$'s own $d=87\,723\,532\,425\,289$ **is prime**.

### 4.3 Step 2: $\sqrt[3]{k}\bmod d'$

- $p\equiv2\pmod3$: bijective cube map, unique root $r\equiv k^{u}$, $3u\equiv1\pmod{p-1}$ — one exp.
- $p\equiv1\pmod3$: write $p-1=3^{t}s$, $3\nmid s$. Take $u$ with $3u\equiv1\pmod s$; then $x_0^3=a\cdot\delta$ with $\delta=(a^s)^{-c}$, $c=(3u-1)/s$, and $\delta$ lies in the 3-Sylow. Correct with $y\in\langle\eta\rangle$, $\eta=h^{s}$ for a non-cubic residue $h$, using base-3 digit extraction of $\log_\eta\delta$, then $y=\eta^{e/3}$. The three roots are $x_0y\omega^{j}$, $\omega=\eta^{3^{t-1}}$. (This is Adleman–Manders–Miller; the digit extraction costs $O(t)$ exponentiations, and $t=v_3(p-1)$ is 1 or 2 for most primes.)
- Hensel lift to $p^{e}$: $r\mapsto r+t p^{j}$ with $t=-((r^3-k)/p^{j})(3r^2)^{-1}\bmod p$.
- $p=3$: the derivative $3r^2$ vanishes mod 3, so ordinary Hensel fails; each root branches three ways and the engine tests the three candidates. (Unreachable for $k=114$ because $3\nmid d$ — but `k=8` exercises it in the tests.)
- $p\mid k$: see the case analysis in §2.4. For $k=114$ only $e\le v_p(k)=1$ survives, so the root class is a single residue.
- CRT to combine prime powers; the engine asserts $r^3\equiv k\pmod{d'}$ at every call site (the same defensive check the upstream kernel keeps).

### 4.4 Steps 4–5: the wheel and the bitmaps

Mod-243 wheel. Precompute `allowed[s mod 243] = {z mod 243 : feasible}` from the 243×243 table of $(x+y, x^3+y^3)$ residues. Then for one curve:

- $z=r+dq$ and $z\equiv a\pmod{243}$. With $g=\gcd(d,243)$, $d'=d/g$, $m=243/g$: solvable only if $a\equiv r\pmod g$, and then $q\equiv((a-r)/g)\,(d')^{-1}\pmod m$. So the $q$ range splits into $|\text{allowed}|$ arithmetic progressions of step $m=243/g$ — **not always 243**, which is the classic bug when $3\mid d$.

Quadratic-residue bitmaps. The value $v^2=(4(k-z^3)-s^3)/(3s)$ must be a square; since $(3s)^{-2}$ is itself a square this is equivalent to

$$3s\,\bigl(4(k-z^3)-s^3\bigr)\ \text{being a QR mod } p$$

for each $p\in\{5,7,11,13,17,19,23,31,37,41,43,47,53,59,61\}$. Precompute `mask[p][s mod p]` as a `u64` bitmask indexed by `z mod p`; the test is one shift and one AND. `s mod p` is computed once per curve; `z mod p` is advanced **incrementally** (`zmod += dmod·step; if zmod >= p { zmod -= p }`) — never a division in the hot loop.

Measured rejection: $1.37\times10^{10}$ quotient points yielded only $28\,183$ exact tests ≈ $1:4.85\times10^5$ (upstream campaign counters), so the hot loop really is adds, table lookups, and bit tests.

Caveat: if $p\mid s$ then $3s\equiv0$ and the QR test degenerates (0 counts as a square); such primes must be skipped for that curve.

### 4.5 Step 6: exact test and word sizes

| quantity | bound ($z\sim6\times10^{21}$) | width |
|---|---|---|
| $d,\ r$ | $\sim1.5\times10^{18}$ (larger if the region moves up) | **u64** |
| $z$ | $>1.8\times10^{19}$ leaves u64; here $\sim10^{20}$ | **u128** |
| $z^3$, $4(k-z^3)-d^3$, quotient | $\sim10^{66}$ | **u256** ($[u64;4]$) |
| $v$, $x$, $y$ | $\sim10^{23}$ | **u128 / i128** |
| re-verification $x^3+y^3+z^3$ | $\sim10^{69}$ | u256 (sign–magnitude) |

u128 alone is insufficient ($3.4\times10^{38}$ caps $z$ at $7\times10^{12}$), hence a hand-rolled u256 with `mul`, `div_rem_u128` (binary long division, remainder in u128 — safe because $3|s|<2^{63}$), and a restoring square root (128 iterations; the root fits u128).

**Forbidden**: `num-bigint`/GMP in the hot path (heap allocation), and any `f64` in the verdict (a cube of $10^{22}$ has only 15 significant digits).

### 4.6 Complexity

- Time: $O(\#\text{points})$ with $\#\text{points}\approx(R_{\max}-R_{\min})\Delta d$, a few ns per point in the hot loop, amortised curve setup (factorisation + cube roots) of order µs.
- Space: base-sieve bitset (76 MB at $d\sim1.5\times10^{18}$), block arrays (~50 MB), 243 wheel (<1 MB), prime masks (15×61×8 B).
- The crucial ratio is **curve setup / points per curve**. With $R$ up to 64 the ratio is favourable; with a narrow $R$ band the setup can dominate. Measure before choosing $\Delta d$ vs $R$ extent.

---

## 5. Architecture

### 5.1 Layout

```
local-computing/
├── Cargo.toml            no external dependencies
├── console.bat           one-click launcher for the web console
├── src/
│   ├── main.rs           CLI, region iteration, resume, JSON events, manifest
│   ├── u256.rs           U256: mul / div_rem_u128 / isqrt / cmp / shifts
│   ├── exact.rs          §2.3 acceptance test (sign–magnitude over U256)
│   ├── filters.rs        mod-243 wheel + 15 prime QR bitmaps
│   ├── roots.rs          cube roots mod p^e (bijections / AMM / Hensel / CRT)
│   ├── dgen.rs           segmented factorisation of a d block
│   └── curve.rs          one (d,r): AP walk + filter chain + exact test
├── app/
│   ├── server.py         local HTTP console (stdlib only), spawns the engine
│   └── static/index.html single page: start/stop, statistics, (z,d) plot
├── tests/
│   └── known.rs          known-solution regression
└── docs/DESIGN.md, docs/screenshot.png
```

There is no `store.rs`: persistence is three plain files in `<out>/` (`hits.jsonl`,
`progress.log`, `region.json`) plus the console's `console.json`, written by `main.rs` and
`app/server.py` respectively. No SQLite, no dependencies.

### 5.2 Data flow

```
Region{k,dlo,dhi,Rlo,Rhi}
  └─> d blocks ──> dgen (factorisation) ──> [(d, roots)]
                                                   │ std::thread::scope over d
                              curve::scan(d, r, s, q-range)
                                                   ↓ hits (rare)
                          hits.jsonl / progress.log / region.json
                                                   │ JSON event per block on stdout
                                     app/server.py ──> /api/state ──> index.html
```

### 5.3 Parallelism

- Granularity: **d-block** (never the $q$ loop), so the wheel and masks are shared and there is no synchronisation in the hot path.
- Per-thread hit vectors and counters, merged at block boundaries.
- Threads: the console default is 6, leaving two of the eight physical cores free, and the engine is launched below normal priority. Benchmarks used 8 threads, but that is the "machine is free" setting, not the daily one; 16 threads is never optimal here because of the power ceiling.

### 5.4 Persistence and self-certification

`<out>/`:

1. `progress.log` — appended, one line per completed d-block with cumulative counts. **This is the resumption state**: `--resume` takes the maximum `d<=` value and starts there.
2. `hits.jsonl` — appended, `{k, x, y, z, d, verified:true}`, flushed right after a hit.
3. `region.json` — the manifest of the *invocation*: $k$, `d_lo`/`d_hi` actually scanned, `resumed_from`, `frontier`, $R$ band, engine version, threads, block size, points, curves, exact tests, hits, provably-empty counts, seconds, throughput, `completed`.
4. `engine.stderr.log` — the engine's human log, when launched from the console.
5. `console.json` — the console's cumulative view of the *region* (scanned intervals, totals, hits), so the picture survives a console restart.

Division of responsibility: `region.json` certifies one invocation; `console.json` certifies the region. Neither is written `completed` unless every $d$ in the range was actually processed, with `dense_skipped = 0` (a `d` whose factorisation overflowed the 16-slot array is *not* scanned, and would have to be declared).

All three files are opened in **append** mode, so a resumed run never destroys the earlier record.

### 5.5 CLI

```
k114 --k 114 --dlo 1.4815e18 --dhi 1.4816e18 \
     --rlo 3.4 --rhi 64 --threads 2 --out runs/20261123-a [--resume] [--emit-json]
```

$z$ bounds are expressed as the ratio $R=z/d$, because $R>3.847$ is intrinsic and ratios are directly comparable with published rectangles.

### 5.6 Web console (`app/`)

The console is deliberately thin, and the engine stays dumb: it knows nothing about
pausing, cancelling or aggregate state. Standard library only, one page on `127.0.0.1`.

| Piece | Responsibility |
|---|---|
| `app/server.py` | `ThreadingHTTPServer`; spawns the engine with `--resume --emit-json` (below normal priority), one reader thread parses its stdout, serves `/` and `/api/{state,start,stop}`, persists `console.json` |
| `app/static/index.html` | single page: Start/Stop, statistics tiles, coverage card, canvas plot; polls `/api/state` every second |

**Event protocol** (one JSON object per line on stdout, only with `--emit-json`):

```
{"event":"start","k":114,"dlo":...,"dhi":...,"rlo":...,"rhi":...,"threads":...,"block":...,"resumed_from":...}
{"event":"block","d_lo":...,"d_hi":...,"points":...,"cum_points":...,"curves":...,"exact":...,"cum_exact":...,"hits":...,"no_roots":...,"div3":...,"dense":...,"secs":...}
{"event":"hit","x":"...","y":"...","z":"...","d":...}
{"event":"done","frontier":...,"points":...,"hits":...,"completed":true|false}
```

Human progress still goes to stderr, so a console-launched run leaves a readable log.

**Why the console owns resumption.** `Stop` is a process termination on Windows, so up to
one block is lost; that is acceptable because a block is 65536 values of $d$ — about 0.2 s
at the target scale. On `Start`, the engine is invoked with `--resume`, which re-reads
`progress.log`. The console therefore only has to aggregate, never to checkpoint.

**The plot.** Axes are $x=z$, $y=d$, both $\log_{10}$; the fitted view is the configured
region, so the picture shows progress against the plan. Each scanned interval
$[d_0,d_1]$ is drawn as the quadrilateral between $z=R_{\min}d$ and $z=R_{\max}d$, which
is the exact set of $(z,d)$ pairs proved not to be solutions. The implied height window
$H=d\,R^{3/2}/\sqrt3$ is reported separately in the coverage card, because the plot's
coordinates cannot show it directly.

Two presentation details that matter at these magnitudes: when the $d$ interval is narrow
relative to $d$ (the normal case — $10^{-4}$ relative width at $10^{18}$), the axis and
ranges are printed as `base + offset`, otherwise every label rounds to the same number;
and the wheel/drag zoom exists because a single auto-fit cannot resolve both the $z$
decades and the $d$ sliver.

---

## 6. Correctness

### 6.1 Known-solution regression (must pass before any long run)

For $k=3,33,42,795$:

1. §2.3 returns the known triple when fed the true $(d,s,z)$;
2. the $H$ formula's error stays below 1% (§2.6 table);
3. given $d$, the known $z$ satisfies $z\equiv r\pmod{d'}$ for one of the computed roots;
4. **end-to-end**: with $d$'s true factorisation, scanning a narrow window around the known $z$ must actually *find* it (no feeding the answer).

Factorisations (already verified, used as test fixtures):

| $k$ | $d$ | factorisation |
|---|---|---|
| 3 | 108 398 887 211 | $167\cdot649\,095\,133$ |
| 33 | 87 723 532 425 289 | prime |
| 42 | 102 980 666 258 459 | $11\cdot43\cdot215\,921\cdot1\,008\,323$ |
| 795 | 21 083 965 616 656 | $2^4\cdot17\cdot149\cdot520\,232\,077$ |

These cover the awkward cases: prime $d$ (hardest cube root), $2^4$ (2-adic and $d'=d/16$), many small factors (CRT fan-out).

### 6.2 Cross-validation

Same $(k,d_{lo},d_{hi},R_{lo},R_{hi})$: the hit *sets* of this engine and `zcubes` must be **equal**, not merely similar. That is the only effective way to catch a sieve that silently drops points.

### 6.3 Assertions kept in release

- $r^3\equiv k\pmod{d'}$ for every computed root;
- $3s\mid 4t-s^3$ on entering the exact test;
- `isqrt` bracket: $v^2\le n<(v+1)^2$;
- final re-verification $x^3+y^3+z^3=k$ for every hit (cost negligible);
- counter balance: `points == rejected_wheel + rejected_parity + rejected_prime + exact_tests`, in the spirit of the upstream kernel's accounting asserts.

---

### 6.4 Implemented coverage (v0.1.0)

27 tests, all green (`cargo test --release`):

| Test | What it pins down |
|---|---|
| `exact::known_solutions_are_accepted` | the four published solutions pass the acceptance test, and the recovered $(x,y)$ is the right *pair* |
| `exact::neighbours_are_rejected` | $\pm2$ perturbations of $z$, and perturbed $d$, never produce an unverified hit |
| `exact::zero_z_and_small_cases` | $2=1^3+1^3+0^3$ and the $\alpha$-extremal $43=9^3-7^3-7^3$ |
| `u256::*` | `mul`, `div_rem_u128` and `isqrt` against `u128` arithmetic; overflow brackets |
| `roots::roots_mod_prime_match_brute_force` | all primes $<600$, all $k<200$, exhaustively |
| `roots::roots_mod_prime_large_prime` | $p=87\,723\,532\,425\,289$ (a prime $d$!), three roots verified |
| `roots::prime_power_lifting` | $r^3\equiv k\pmod{p^e}$ plus completeness against brute force |
| `roots::factorization_roots_are_consistent` | each known solution's signed $z$ lies on a computed root class |
| `filters::wheel_matches_brute_force` | the whole $243\times243$ wheel against direct enumeration |
| `filters::known_solutions_pass_every_filter` | no filter may reject a real solution (all four $k$) |
| `curve::finds_known_solutions` | end-to-end recovery of all four, with balanced accounting |
| `curve::three_dividing_d_is_dead_for_k_equiv_3_or_6_mod_9` | the §2.4 lemma |
| `curve::wheel_step_handles_divisible_by_three` | the $\gcd(d',243)\ne1$ branch (via $k=8$) |
| `curve::icbrt_is_exact` | the window split at $k^{1/3}$ |
| `dgen::small_ranges` | segmented factorisation equals naive factorisation on 60k values |
| `known::brute_force_agrees_on_small_k` | for $k\in\{2,6,10,43,71,100,128\}$ every brute-force solution with $d\le240$, $R\in[3.4,64]$ is found |
| `known::no_false_positives_in_a_small_window` | $k=114$, $d=1\,000\,003$ prime: hits verify or there are none |

## 7. Performance

### 7.1 Build

```toml
[profile.release]
opt-level = 3
lto = "fat"
codegen-units = 1
```

`RUSTFLAGS=-C target-cpu=native`. Target machine: Ryzen 7 7735HS (Zen3+, 8C/16T, 16 MB L3, AVX2, no AVX-512).

### 7.2 Benchmarking

1. Measure the reference: `zcubes 8 57 1 1e9 1e9 1e10` (README claims ~10 s, 15 solutions) — derive its points/s on this machine.
2. Run this engine over the same rectangle; compare points/s *and* hit sets.
3. Microbenchmark separately: curve setup (factorisation + roots) vs the $q$ loop.

### 7.3 Targets, and the measured baseline

Baseline of v0.1.0 (see §3.6): **$2.8\times10^{6}$ $d$ per second**, $6.8\times10^{6}$ effective points per second, 8 threads, $d\sim10^{12}$.

Same engine at the *actual* target scale, $d\approx1.4815\times10^{18}$: **$3.8\times10^{5}$ $d$ per second** (4.2e8 values in 1095 s, 2.4 effective points per $d$, 1 in $2.5\times10^{4}$ reaching the exact test). That is **7× slower** than at $d\sim10^{12}$, and the cause is structural: exact segmented factorisation needs every prime up to $\sqrt{d_{\max}}=1.2\times10^{9}$, i.e. $6.2\times10^{7}$ primes, each of which pays a "first multiple in this block" computation *for every block*. With 65536-value blocks that is 6.2e7 × 65536⁻¹ ≈ 1000 prime-visits per $d$ — pure overhead.

So the first real optimisation is not the $q$ loop at all, but replacing this enumeration with the reference implementation's scheme: **enumerate $d=p\cdot m$ by its largest prime factor $p$ in a band $[p_{\min},p_{\max}]$**, so no prime list up to $\sqrt{d}$ is ever needed. That is almost certainly why `zcubes` exposes `pmin/pmax` (open question 2). Consequence for coverage honesty: a prime-band enumeration is *not* contiguous in $d$; the manifest must then state the band instead of an interval.

Profiling by arithmetic (the engine has no per-phase timing yet) points at per-$d$ setup:

- factorisation is amortised by the block sieve, and the prime-streaming dominates at large $d$ (primes up to $1.2\times10^{9}$ for $d\sim1.5\times10^{18}$);
- cube roots cost one `powmod` per $p\equiv2\pmod3$ and several per $p\equiv1\pmod3$ (AMM), i.e. ~1–3 µs per $d$;
- the $q$ loop itself only sees ~2.4 points per $d$ after the wheel.

Ordered optimisation list, to be validated by measurement:

1. **Root cache** keyed by $(p,e)$: a prime $p$ recurs $\Delta d/p$ times inside a block, so the hit rate is high for all but the largest primes.
2. **Prime-band restriction** (`pmin`/`pmax`, as `zcubes` exposes): bounds the root cost and makes caching near-perfect, at the price of a stated coverage restriction.
3. Composite wheel ($3^5\cdot5\cdot7$) to cut the per-curve wheel setup.
4. Redundant-scans removal: 45% of $d$ have no root at all — cheap to detect earlier by sieving residues mod a few small primes.

Targets: correctness first (equal hit sets), then ≥2× `zcubes` on the same rectangle; 3–5× is plausible only after items 1–2.

---

## 8. Risks and traps

| Trap | Note |
|---|---|
| u128 overflow | $z>7\times10^{12}$ already exceeds u128 for cubes; fix the width ladder up front |
| Rust `%` on negatives | truncates toward zero; use `rem_euclid` for residues, or keep magnitudes unsigned with an explicit sign |
| Sign of $s$ | forced by $\operatorname{sign}(k-z^3)$; never loop over both signs |
| $3\mid d$ in the wheel | AP step is $243/\gcd(d,243)$, not 243 |
| $p\mid s$ in the QR masks | $3s\equiv0$ makes the test vacuous; skip that prime for the curve |
| $9\mid d'$ | provably no roots; skip (saves ~1/9 of the $d$) |
| Big integers in the hot loop | any allocation in the $q$ loop destroys throughput by orders of magnitude |
| Floating point | not in the verdict, not in the re-verification, not in logging large integers |
| Laptop power ceiling | measure with 8 threads, not 16 |
| Bit-sieve bandwidth | long $q$ runs become memory-bound; block for L2/L3 |
| Calling sampling "scanning" | a region must be a stated $(d,z)$ rectangle with no skipped $d$, or the claim is void |
| Re-deriving $\rho$ from first principles | the naive integral is off by a factor ≈8.5 (singular series); always use the tabulated $\rho$ |

---

## 9. Milestones

| # | Deliverable | Acceptance | Status |
|---|---|---|---|
| M0 | this document | — | **done** |
| M1 | `u256.rs` + `exact.rs` | unit tests pass; all four known solutions accepted | **done** |
| M2 | `filters.rs` + `curve.rs` (no bit-sieve yet) | end-to-end recovery of all four; balanced accounting | **done** (bit-sieve deferred to M5) |
| M3 | `roots.rs` + `dgen.rs` (contiguous $d$ enumeration) | arbitrary $d$ intervals, no smoothness restriction | **done** |
| M4 | parallelism + progress + manifest | kill and restart loses nothing | **done** |
| M4.5 | `--resume` reading `progress.log`; append-only outputs | an interrupted run continues at the last completed block | **done** |
| M4.6 | web console: Start/Stop, statistics, $(z,d)$ coverage plot | page drives a real run; state survives a console restart | **done** |
| M5 | root cache, prime bands, composite wheel, bit-sieve | ≥2× the v0.1.0 baseline, benchmarked | **next** |
| M6 | cross-check hit sets against `zcubes` on one small rectangle | sets are equal, not merely similar | pending |
| M7 | long run, sized by measured throughput | `region.json` says `completed: true` | on demand (the console makes this a per-session decision) |

---

## 10. Open questions

1. `dgen`: segmented factorisation (chosen) vs DFS with a smoothness bound — measure the per-$d$ overhead at $d\sim1.5\times10^{18}$.
2. Why exactly does `zcubes` expose `pmin/pmax`? Performance (cheap roots, tables) or a correctness reason? Read `cbrts.h` / `admissible.c`.
3. Is a composite wheel ($3^5\cdot5\cdot7=8505$, ≈1 KB) worth it for mid-length $q$ runs?
4. AVX2 8-lane residue accumulation: only after profiling says the $q$ loop dominates (setup may dominate instead).
5. Optimal $R_{\max}$ at our $d_{\max}$: B–S measured 64 for $d_{\max}\ge2^{50}$; re-derive for $d\sim1.5\times10^{18}$.
6. Root caching: for a window of width $\Delta d$, a prime $p$ recurs $\Delta d/p$ times; cache roots for the large primes and measure. **This is the measured bottleneck (§7.3) and should be done first.**
7. Implement route B (Mordell norm form) as a second engine for complementary sampling? Nice-to-have.
8. Should the default target be the open $k$ with the largest $\rho$ instead of 114?
9. Resume semantics: the current `progress.log` records completed blocks; a real `--resume` flag that skips them is not implemented yet (re-running a sub-window is the workaround).

---

## 11. References

- A. R. Booker, A. V. Sutherland, *On a question of Mordell*, PNAS 118(11):e2022377118, 2021 — [arXiv:2007.01209](https://arxiv.org/abs/2007.01209). §2 (parameterisation, height), §3.8 ($d$-exclusions), §4 (cost model, optimal $R$), §5.2 (all rounds), Table 1 ($\rho(k)$).
- A. V. Sutherland, lecture notes *Sums of three cubes* — complexity ladder $O(N^3)\to O(N^{2+o(1)})\to O(N^{1+o(1)})$.
- Reference implementation (C, MIT, primesieve + GMP): <https://github.com/AndrewVSutherland/SumsOfThreeCubes>
- S. G. Huisman, *Newer sums of three cubes*, 2016 — [arXiv:1604.07746](https://arxiv.org/abs/1604.07746)
- A.-S. Elsenhans, J. Jahnel, *New sums of three cubes*, 2009 (Elkies lattice reduction)
- A. R. Booker, *Cracking the problem with 33*, 2019 — [arXiv:1903.04284](https://arxiv.org/abs/1903.04284)
- Upstream kernel used only to reverse-engineer the sampling scheme and to borrow measured rejection counters: `../k=114/math-gambling/tools/search_core.py` (GPL-2.0, not reused in code)

### Explicitly excluded sources

- Zenodo "Universal Cubic Pattern" papers claiming to settle 33/42/114/390/…: no verifiable content.
- Items claiming a solution of 114 that actually solve $114^3=18^3+60^3+108^3$ — that is "three cubes summing to a cube", a different problem.
