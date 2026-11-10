//! MS-Numpress encoding and decoding in safe Rust: linear prediction (m/z, retention
//! time), positive integer and short logged float (intensities).
//!
//! Ported from the reference implementation, ms-numpress `MSNumpress.cpp` (Copyright
//! 2013 Johan Teleman, Apache-2.0). Encoders produce its bytes exactly. Decoders check
//! every read, so corrupt input is an [`Error`], never a panic. Values are computed in
//! f64 as there; decoding can store them as f32 or f64 (see [`Float`]).
//!
//! Every function appends to its output buffer. An encoder that fails leaves the buffer
//! as it found it.

#![forbid(unsafe_code)]

mod half_bytes;
mod linear;
mod pic;
mod slof;

pub use linear::{decode_linear, encode_linear, optimal_linear_fixed_point};
pub use pic::{decode_pic, encode_pic};
pub use slof::{decode_slof, encode_slof, optimal_slof_fixed_point};

use std::fmt;

#[doc = include_str!("../README.md")]
#[cfg(doctest)]
pub struct ReadmeDoctests;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Error {
    /// The array is shorter than its header.
    TooShort(usize),
    Corrupt,
    /// A value does not fit the codec at the given fixed point (or is NaN).
    Overflow,
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TooShort(length) => write!(f, "numpress array of {length} bytes is too short"),
            Self::Corrupt => write!(f, "corrupt numpress array"),
            Self::Overflow => write!(f, "value cannot be numpress-encoded at this fixed point"),
        }
    }
}

impl std::error::Error for Error {}

/// The type decoded values are stored as.
pub trait Float: sealed::Sealed + Copy {
    fn from_f64(value: f64) -> Self;
}

impl Float for f32 {
    fn from_f64(value: f64) -> Self {
        value as f32
    }
}

impl Float for f64 {
    fn from_f64(value: f64) -> Self {
        value
    }
}

mod sealed {
    pub trait Sealed {}
    impl Sealed for f32 {}
    impl Sealed for f64 {}
}

/// A numpress codec, as named by mzML's binary-array compression CV terms.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Codec {
    Linear,
    Pic,
    Slof,
}

impl Codec {
    pub fn decode<T: Float>(self, data: &[u8], out: &mut Vec<T>) -> Result<(), Error> {
        match self {
            Self::Linear => decode_linear(data, out),
            Self::Pic => decode_pic(data, out),
            Self::Slof => decode_slof(data, out),
        }
    }

    /// Encodes at the reference's optimal fixed point for `values`.
    pub fn encode(self, values: &[f64], out: &mut Vec<u8>) -> Result<(), Error> {
        match self {
            Self::Linear => encode_linear(values, optimal_linear_fixed_point(values), out),
            Self::Pic => encode_pic(values, out),
            Self::Slof => encode_slof(values, optimal_slof_fixed_point(values), out),
        }
    }
}

const FIXED_POINT_BYTES: usize = 8;

/// The fixed point is stored big-endian.
fn read_fixed_point(data: &[u8]) -> Result<f64, Error> {
    let bytes = data
        .get(..FIXED_POINT_BYTES)
        .ok_or(Error::TooShort(data.len()))?;
    Ok(f64::from_be_bytes(bytes.try_into().expect("8 bytes")))
}

fn write_fixed_point(fixed_point: f64, out: &mut Vec<u8>) {
    out.extend_from_slice(&fixed_point.to_be_bytes());
}

fn all_or_nothing(
    out: &mut Vec<u8>,
    encode: impl FnOnce(&mut Vec<u8>) -> Result<(), Error>,
) -> Result<(), Error> {
    let start = out.len();
    let result = encode(out);
    if result.is_err() {
        out.truncate(start);
    }
    result
}
