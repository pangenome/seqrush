//! Fixed hashing seeds for stable hashing.
//!
//! Set to [`None`] to disable stable hashing.
//!
//! See [`rhai::config::hashing::set_hashing_seed`][set_hashing_seed].
//!
//! # Example
//!
//! ```rust
//! // Set the hashing seed to [1, 2, 3, 4]
//! rhai::config::hashing::set_hashing_seed(Some([1, 2, 3, 4])).unwrap();
//! ```
//! Alternatively, set this at compile time via the `RHAI_HASHING_SEED` environment variable.
//!
//! # Example
//!
//! ```sh
//! env RHAI_HASHING_SEED ="[236,800,954,213]"
//! ```

use super::hashing_env;
use crate::OnceCell;

static HASHING_SEED: OnceCell<Option<[u64; 4]>> = OnceCell::new();

#[allow(deprecated)]
pub use crate::api::deprecated::config::hashing::{get_ahash_seed, set_ahash_seed};

/// Set the hashing seed. This is used to hash functions etc.
///
/// This is a static global value and affects every Rhai instance.
/// This should not be used _unless_ you know you need it.
///
/// Set the hashing seed to all zeros effectively disables stable hashing.
///
/// # Warning
///
/// * You can only call this function **ONCE** for the entire duration of program execution.
/// * You **MUST** call this before performing **ANY** Rhai operation (e.g. creating an [`Engine`][crate::Engine]).
///
/// # Error
///
/// Returns an error containing the existing hashing seed if already set.
///
/// # Example
///
/// ```rust
/// # use rhai::Engine;
/// // Set the hashing seed to [1, 2, 3, 4]
/// rhai::config::hashing::set_hashing_seed(Some([1, 2, 3, 4])).unwrap();
///
/// // Use Rhai AFTER setting the hashing seed
/// let engine = Engine::new();
/// ```
#[inline(always)]
pub fn set_hashing_seed(new_seed: Option<[u64; 4]>) -> Result<(), Option<[u64; 4]>> {
    #[cfg(feature = "std")]
    return HASHING_SEED.set(new_seed);

    #[cfg(not(feature = "std"))]
    return HASHING_SEED.set(new_seed.into()).map_err(|err| *err);
}

/// Get the current hashing Seed.
///
/// If the seed is not yet defined, the `RHAI_HASHING_SEED` environment variable (if any) is used.
///
/// Otherwise, the hashing seed is randomized to protect against DOS attacks.
///
/// See [`rhai::config::hashing::set_hashing_seed`][set_hashing_seed] for more.
#[inline]
#[must_use]
pub fn get_hashing_seed() -> &'static Option<[u64; 4]> {
    HASHING_SEED.get().unwrap_or(&hashing_env::HASHING_SEED)
}
