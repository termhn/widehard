//! # `ultraviolet`
//!
//! This is a crate to computer-graphics and games-related linear and geometric algebra, but *fast*, both in terms
//! of productivity and in terms of runtime performance.
//!
//! In terms of productivity, ultraviolet uses no generics and is designed to be as straightforward
//! of an interface as possible, resulting in fast compilation times and clear code. In addition, the
//! lack of generics and Rust type-system "hacks" result in clear and concise errors that are easy to
//! parse and fix for the user.
//!
//! In terms of runtime performance, ultraviolet was designed from the start with performance in mind.
//! To do so, we provide two separate kinds of each type, each with nearly identical functionality,
//! one with usual scalar f32 values, and the other a 'wide' type which uses SIMD f32x4 vectors for
//! each value. This design is clear and explicit in intent, and it also allows code to
//! take full advantage of SIMD.
//!
//! The 'wide' types use an "SoA" (Structure of Arrays) architecture
//! such that each wide data structure actually contains the data for 4 or 8 of its associated data type and will do any operation
//! on all of the simd 'lanes' at the same time. For example, a `Vec3x8` is equivalent to 8 `Vec3`s all bundled together into one
//! data structure.
//!
//! Doing this is potentially *much* (factor of 10) faster than an standard "AoS" (Array of Structs) layout,
//! though it does depend on your workload and algorithm requirements. Algorithms must be carefully architected to take full advantage
//! of this, and doing so can be easier said than done, especially if your algorithm involves significant branching.
//!
//! `ultraviolet` was the first Rust math library to be designed in this "AoSoA" manner, though
//! `nalgebra` now supports it for several of their data structures as well.
//!
//! ## Benchmarks
//!
//! See [`mathbench-rs`](https://github.com/bitshifter/mathbench-rs) for latest benchmarks.
//!
//! ## Cargo Features
//!
//! To help further improve build times, `ultraviolet` puts various functionality under feature flags. For example, the 2d and 3d projective geometric algebras
//! as well as f64 and integer types are disabled by default. In order to enable them, enable the corresponding crate feature flags in your `Cargo.toml`. For example:
//!
//! ```toml
//! [dependencies]
//! ultraviolet = { version = "0.9", features = [ "f64", "int" ] }
//! ```
//!
//! Will enable the `f64` and `int` features. Here's a list of the available features:
//!
//! * `f64` – Enable `f64` bit wide floating point support. Naming convention is `D[Type]`, such as `DVec3x4` would be a collection of 4 3d vectors with `f64` precision each.
//! * `int` – Enable integer vector types.
//! * `bytemuck` – Enable casting of many types to byte arrays, for use with graphics APIs.
//! * `mint` – Enable interoperation with other math crates through the `mint` interface.
//! * `num-traits` – Enable [identity traits](https://docs.rs/num-traits/latest/num_traits/identities/index.html) for interoperation with other math crates.
//! * `serde` – Enable `Serialize` and `Deserialize` implementations for many scalar types.
//!
//! ## Crate Features
//!
//! This crate is currently being dogfooded in my ray tracer [`rayn`](https://github.com/termhn/rayn),
//! and is being used by various independent Rust game developers for various projects.
//! It does what those users have currently needed it to do.
//!
//! There are a couple relatively unique/novel features in this library, the most important being the use of the Geometric Algebra.
//!
//! Instead of implementing complex number algebra (for 2d rotations) and Quaternion algebra (for 3d rotations), we use
//! Rotors, a concept taken from Geometric Algebra, to represent 2d and 3d rotations.
//!
//! What this means for the programmer is that you will be using the `Rotor3` type in place of
//! a Quaternion, though you can expect it to do basically all the same things that a Quaternion does. In fact, Quaternions
//! are directly isomorphic to Rotors (meaning they are in essense the same thing, just formulated differently). The reason this decision was made was twofold:
//! first, the derivation of the math is actually quite simple to understand. All the derivations for the code implemented in the Rotor structs in this
//! library are written out in the `derivations` folder of the GitHub repo; I derived them manually as part of the implementation.
//!
//! On the other hand, Quaternions are often basically just seen as black boxes that we programmers use to do rotations because
//! they have some nice properties, but that we don't really understand. You can use Rotors this same way, but you can also easily
//! understand them. Second is that in some sense they can be seen as 'more correct' than Quaternions. Specifically, they
//! facilitate a more proper understanding of rotation as being something that occurs *within a plane* rather than something
//! that occurs *around an axis*, as it is generally thought. Finally, Rotors also generalize to 4 and even higher dimensions,
//! and if someone wants to they could implement a Rotor4 which retains all the properties of a Rotor3/Quaternion but does rotation
//! in 4 dimensions instead, something which simply is not possible to do with Quaternions.
//!
//! If it's missing something you need it to do, bug me on the [GitHub issue tracker](https://github.com/termhn/ultraviolet/issues) and/or Rust community discord server
//! (I'm Fusha there) and I'll try to add it for you, if I believe it fits with the vision of the lib :)

#![deny(
    rust_2018_compatibility,
    rust_2018_idioms,
    future_incompatible,
    nonstandard_style,
    unused,
    clippy::all
)]

extern crate alloc;
#[cfg(feature = "serde")]
extern crate serde;

mod util;

pub(crate) use util::Splat;

pub mod bivec;
#[cfg(feature = "int")]
pub mod conversion;
#[cfg(feature = "int")]
pub mod int;
pub mod interp;
pub mod mat;
pub mod projection;
pub mod rotor;
pub mod traits;
pub mod transform;
pub mod vec;

#[cfg(feature = "serde")]
mod impl_serde;
#[cfg(feature = "serde")]
pub use impl_serde::*;

#[cfg(feature = "mint")]
mod impl_mint;
#[cfg(feature = "mint")]
pub use impl_mint::*;

#[cfg(feature = "bytemuck")]
mod impl_bytemuck;

#[cfg(feature = "bytemuck")]
pub use impl_bytemuck::*;

pub use bivec::*;
#[cfg(feature = "int")]
pub use conversion::*;
#[cfg(feature = "int")]
pub use int::*;
pub use interp::*;
pub use mat::*;
pub use rotor::*;
pub use transform::*;
pub use vec::*;

pub use wide::f32x4;
pub use wide::f32x8;
pub use wide::f64x2;
pub use wide::f64x4;

pub use wide::f32x4 as m32x4;
pub use wide::f32x8 as m32x8;
pub use wide::f64x2 as m64x2;
pub use wide::f64x4 as m64x4;

pub(crate) use wide::{CmpGe, CmpLt};
