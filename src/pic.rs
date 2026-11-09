//! Positive integer (intensities): each value rounded to an integer, stored as a
//! half-byte int. Each value is off by at most 0.5.

use crate::half_bytes::{HalfByteReader, HalfByteWriter};
use crate::{all_or_nothing, Error, Float};

/// The reference's range: rounds to `0..=i32::MAX`.
fn rounded(value: f64) -> Result<u32, Error> {
    if !(value >= -0.5 && value + 0.5 <= i32::MAX as f64) {
        return Err(Error::Overflow);
    }
    Ok((value + 0.5) as u32)
}

pub fn encode_pic(values: &[f64], out: &mut Vec<u8>) -> Result<(), Error> {
    all_or_nothing(out, |out| {
        let mut ints = HalfByteWriter::new(out);
        for &value in values {
            ints.push_int(rounded(value)?);
        }
        ints.finish();
        Ok(())
    })
}

pub fn decode_pic<T: Float>(data: &[u8], out: &mut Vec<T>) -> Result<(), Error> {
    let mut ints = HalfByteReader::new(data);
    while !ints.at_end() {
        out.push(T::from_f64(ints.next_int()? as f64));
    }
    Ok(())
}
