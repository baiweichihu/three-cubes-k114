//! Local search kernel for `x^3 + y^3 + z^3 = k`, default `k = 114`.
//!
//! See `docs/DESIGN.md` for the mathematics and the architecture.

pub mod curve;
pub mod dgen;
pub mod exact;
pub mod filters;
pub mod roots;
pub mod u256;
