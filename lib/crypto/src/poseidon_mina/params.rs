//! This module contains a trait with poseidon hash parameters.
//!
//! Consumer of this trait should implement the parameters for the specific
//! poseidon hash instance.
//! Or use the existing instances in the [`crate::poseidon_mina::instance`]
//! module.

use crate::field::prime::PrimeField;

/// Poseidon hash parameters.
pub trait PoseidonMinaParams<F: PrimeField> {
    const SPONGE_CAPACITY: usize = 1;
    const SPONGE_WIDTH: usize = 3;
    const SPONGE_RATE: usize = 2;
    const PERM_ROUNDS_FULL: usize;
    const PERM_ROUNDS_PARTIAL: usize;
    const PERM_HALF_ROUNDS_FULL: usize;
    const PERM_SBOX: u32;
    const PERM_FULL_MDS: bool;
    const PERM_INITIAL_ARK: bool;

    /// MDS (Maximum Distance Separable) matrix used in the Poseidon
    /// permutation.
    const MAT_INTERNAL_DIAG_M_1: &'static [&'static [F]];

    /// The round constants used in the full and partial rounds of the Poseidon
    /// permutation.
    const ROUND_CONSTANTS: &'static [&'static [F]];
}
