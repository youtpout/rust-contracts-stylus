//! Base on mina kimchi implementation
//! https://github.com/o1-labs/proof-systems/blob/94982ca896de874d6439e211e4b7d8fb91caa355/poseidon/src/pasta/fp_kimchi.rs#L13
//! Tweek by youtpout to be compatible with stylus
//! Also an other implementation in rust but not optimized https://github.com/youtpout/poseidon_hash_mina/blob/14be42efe898efe3ac806d71223a3f7a237ef5fc/src/lib.rs#L23
pub mod instance;
pub mod params;

use alloc::{boxed::Box, vec, vec::Vec};

use crate::{
    field::prime::PrimeField, poseidon_mina::params::PoseidonMinaParams,
};

#[derive(Clone, Debug)]
pub enum SpongeState {
    Absorbed(usize),
    Squeezed(usize),
}

#[derive(Clone, Default, Debug)]
pub struct ArithmeticSpongeParams<F: PrimeField> {
    pub round_constants: &'static [&'static [F]],
    pub mds: &'static [&'static [F]],
}

/// Poseidon mina sponge that can absorb any number of `F` field elements and be
/// squeezed to a finite number of `F` field elements.
///
/// ## Security Notice
///
/// This is a low-level primitive that does not implement padding or domain
/// separation. Users must ensure proper input formatting and security practices
/// for their specific cryptographic protocols.
#[derive(Clone, Debug)]
pub struct PoseidonMina<P: PoseidonMinaParams<F>, F: PrimeField> {
    phantom: core::marker::PhantomData<P>,
    rate: usize,
    state: Box<[F]>,
    pub sponge_state: SpongeState,
    params: ArithmeticSpongeParams<F>,
}

impl<P: PoseidonMinaParams<F>, F: PrimeField> Default for PoseidonMina<P, F> {
    fn default() -> Self {
        Self::new()
    }
}

impl<P: PoseidonMinaParams<F>, F: PrimeField> PoseidonMina<P, F> {
    /// Create a new Poseidon sponge.
    #[must_use]
    #[inline]
    pub fn new() -> Self {
        let params = ArithmeticSpongeParams {
            round_constants: &P::ROUND_CONSTANTS,
            mds: &P::MAT_INTERNAL_DIAG_M_1,
        };
        Self {
            phantom: core::marker::PhantomData,
            rate: P::SPONGE_RATE,
            state: vec![F::zero(); P::SPONGE_WIDTH].into_boxed_slice(),
            sponge_state: SpongeState::Absorbed(0),
            params,
        }
    }

    /// Size of poseidon sponge's state.
    #[must_use]
    pub const fn state_size() -> usize {
        P::SPONGE_WIDTH
    }

    #[must_use]
    pub const fn rate_size() -> usize {
        P::SPONGE_RATE
    }

    /// Start index of partial rounds.
    ///
    /// This represents the point where the algorithm transitions from full
    /// rounds to partial rounds in the Poseidon permutation.
    #[must_use]
    const fn partial_round_start() -> usize {
        P::PERM_HALF_ROUNDS_FULL
    }

    /// End index of partial rounds (noninclusive).
    ///
    /// This represents the point where the algorithm transitions from partial
    /// rounds back to full rounds in the Poseidon permutation.
    #[must_use]
    const fn partial_round_end() -> usize {
        Self::partial_round_start() + P::PERM_ROUNDS_PARTIAL
    }

    /// Total number of rounds.
    ///
    /// This is the sum of full rounds and partial rounds in the Poseidon
    /// permutation.
    #[must_use]
    const fn rounds() -> usize {
        P::PERM_ROUNDS_FULL + P::PERM_ROUNDS_PARTIAL
    }

    pub fn full_round<SC: PoseidonMinaParams<F>>(&mut self, r: usize) {
        full_round::<F, SC>(&self.params, &mut self.state, r);
    }

    pub fn poseidon_block_cipher<SC: PoseidonMinaParams<F>>(&mut self) {
        poseidon_block_cipher::<F, SC>(&self.params, &mut self.state);
    }

    /// Absorb a single element into the sponge.
    ///
    /// Transitions from [`Mode::Absorbing`] to [`Mode::Squeezing`] mode are
    /// unidirectional.
    ///
    /// # Panics
    ///
    /// May panic if absorbing while squeezing.
    #[inline]
    pub fn absorb(&mut self, elem: &[F]) {
        for x in elem.iter() {
            match self.sponge_state {
                SpongeState::Absorbed(n) => {
                    if n == self.rate {
                        self.poseidon_block_cipher::<P>();
                        self.sponge_state = SpongeState::Absorbed(1);
                        self.state[0].add_assign(x);
                    } else {
                        self.sponge_state = SpongeState::Absorbed(n + 1);
                        self.state[n].add_assign(x);
                    }
                }
                SpongeState::Squeezed(_n) => {
                    self.state[0].add_assign(x);
                    self.sponge_state = SpongeState::Absorbed(1);
                }
            }
        }
    }

    /// Squeeze a single element from the sponge.
    ///
    /// When invoked from [`Mode::Absorbing`] mode, this function triggers a
    /// permutation and transitions to [`Mode::Squeezing`] mode.
    #[inline]
    pub fn squeeze(&mut self) -> F {
        match self.sponge_state {
            SpongeState::Squeezed(n) => {
                if n == self.rate {
                    self.poseidon_block_cipher::<P>();
                    self.sponge_state = SpongeState::Squeezed(1);
                    self.state[0]
                } else {
                    self.sponge_state = SpongeState::Squeezed(n + 1);
                    self.state[n]
                }
            }
            SpongeState::Absorbed(_n) => {
                self.poseidon_block_cipher::<P>();
                self.sponge_state = SpongeState::Squeezed(1);
                self.state[0]
            }
        }
    }
}

fn apply_mds_matrix<F: PrimeField, SC: PoseidonMinaParams<F>>(
    params: &ArithmeticSpongeParams<F>,
    state: &Box<[F]>,
) -> Box<[F]> {
    if SC::PERM_FULL_MDS {
        let v: Vec<F> = params
            .mds
            .iter()
            .map(|m| {
                state
                    .iter()
                    .zip(m.iter())
                    .fold(F::zero(), |x, (s, &m)| m * s + x)
            })
            .collect();
        return v.into_boxed_slice();
    } else {
        vec![state[0] + state[2], state[0] + state[1], state[1] + state[2]]
            .into_boxed_slice()
    }
}

/// Apply a full round of the permutation.
/// A full round is composed of the following steps:
/// - Apply the S-box to each element of the state.
/// - Apply the MDS matrix to the state.
/// - Add the round constants to the state.
///
/// The function has side-effect and the parameter state is modified.
pub fn full_round<F: PrimeField, SC: PoseidonMinaParams<F>>(
    params: &ArithmeticSpongeParams<F>,
    state: &mut Box<[F]>,
    r: usize,
) {
    for state_i in state.iter_mut() {
        *state_i = sbox::<F, SC>(*state_i);
    }
    *state = apply_mds_matrix::<F, SC>(params, state);
    for (i, x) in params.round_constants[r].iter().enumerate() {
        state[i].add_assign(x);
    }
}

pub fn half_rounds<F: PrimeField, SC: PoseidonMinaParams<F>>(
    params: &ArithmeticSpongeParams<F>,
    state: &mut Box<[F]>,
) {
    for r in 0..SC::PERM_HALF_ROUNDS_FULL {
        for (i, x) in params.round_constants[r].iter().enumerate() {
            state[i].add_assign(x);
        }
        for state_i in state.iter_mut() {
            *state_i = sbox::<F, SC>(*state_i);
        }
        let res = apply_mds_matrix::<F, SC>(params, state);
        for (i, state_i) in state.iter_mut().enumerate() {
            *state_i = res[i]
        }
    }

    for r in 0..SC::PERM_ROUNDS_PARTIAL {
        for (i, x) in params.round_constants[SC::PERM_HALF_ROUNDS_FULL + r]
            .iter()
            .enumerate()
        {
            state[i].add_assign(x);
        }
        state[0] = sbox::<F, SC>(state[0]);
        let res = apply_mds_matrix::<F, SC>(params, state);
        res.iter().enumerate().for_each(|(i, x)| {
            state[i] = *x;
        });
    }

    for r in 0..SC::PERM_HALF_ROUNDS_FULL {
        for (i, x) in params.round_constants
            [SC::PERM_HALF_ROUNDS_FULL + SC::PERM_ROUNDS_PARTIAL + r]
            .iter()
            .enumerate()
        {
            state[i].add_assign(x);
        }
        for state_i in state.iter_mut() {
            *state_i = sbox::<F, SC>(*state_i);
        }
        let res = apply_mds_matrix::<F, SC>(params, state);
        res.iter().enumerate().for_each(|(i, x)| {
            state[i] = *x;
        });
    }
}

pub fn poseidon_block_cipher<F: PrimeField, SC: PoseidonMinaParams<F>>(
    params: &ArithmeticSpongeParams<F>,
    state: &mut Box<[F]>,
) {
    if SC::PERM_HALF_ROUNDS_FULL == 0 {
        if SC::PERM_INITIAL_ARK {
            for (i, x) in params.round_constants[0].iter().enumerate() {
                state[i].add_assign(x);
            }
            for r in 0..SC::PERM_ROUNDS_FULL {
                full_round::<F, SC>(params, state, r + 1);
            }
        } else {
            for r in 0..SC::PERM_ROUNDS_FULL {
                full_round::<F, SC>(params, state, r);
            }
        }
    } else {
        half_rounds::<F, SC>(params, state);
    }
}

pub fn sbox<F: PrimeField, SC: PoseidonMinaParams<F>>(mut x: F) -> F {
    if SC::PERM_SBOX == 7 {
        // This is much faster than using the generic `pow`. Hard-code to get
        // the ~50% speed-up that it gives to hashing.
        let mut square = x;
        square.square_in_place();
        x *= square;
        square.square_in_place();
        x *= square;
        x
    } else {
        x.pow(SC::PERM_SBOX)
    }
}
