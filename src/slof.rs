//! Short logged float (intensities): a fixed point, then little-endian u16 values `x`
//! standing for `exp(x / fixed_point) − 1`.
//!
//! `ln(v + 1)` is written as such, not `ln_1p`, to round exactly as the reference does.

use crate::{all_or_nothing, read_fixed_point, write_fixed_point, Error, Float, FIXED_POINT_BYTES};

/// The largest fixed point at which every value fits a u16.
pub fn optimal_slof_fixed_point(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut largest: f64 = 1.0;
    for value in values {
        largest = largest.max((value + 1.0).ln());
    }
    (u16::MAX as f64 / largest).floor()
}

/// The reference rejects only `> u16::MAX`; below `-1.5` (or NaN) its cast is undefined.
fn short_logged(value: f64, fixed_point: f64) -> Result<u16, Error> {
    let scaled = (value + 1.0).ln() * fixed_point;
    if !(scaled > -1.5 && scaled <= u16::MAX as f64) {
        return Err(Error::Overflow);
    }
    Ok((scaled + 0.5) as u16)
}

pub fn encode_slof(values: &[f64], fixed_point: f64, out: &mut Vec<u8>) -> Result<(), Error> {
    all_or_nothing(out, |out| {
        write_fixed_point(fixed_point, out);
        for &value in values {
            out.extend_from_slice(&short_logged(value, fixed_point)?.to_le_bytes());
        }
        Ok(())
    })
}

pub fn decode_slof<T: Float>(data: &[u8], out: &mut Vec<T>) -> Result<(), Error> {
    let fixed_point = read_fixed_point(data)?;
    let (pairs, odd_byte) = data[FIXED_POINT_BYTES..].as_chunks::<2>();
    if !odd_byte.is_empty() {
        return Err(Error::Corrupt);
    }
    out.extend(pairs.iter().map(|&pair| {
        let x = u16::from_le_bytes(pair);
        T::from_f64((x as f64 / fixed_point).exp() - 1.0)
    }));
    Ok(())
}
