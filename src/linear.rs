//! Linear prediction (m/z, retention time): a fixed point, the first two scaled values
//! as little-endian u32, then each scaled value's difference from `2·prev − prev_prev`
//! as a half-byte int. Each value is off by at most `0.5 / fixed_point`.

use crate::half_bytes::{HalfByteReader, HalfByteWriter};
use crate::{all_or_nothing, read_fixed_point, write_fixed_point, Error, Float, FIXED_POINT_BYTES};

const FIRST_SEED: usize = FIXED_POINT_BYTES;
const SECOND_SEED: usize = FIRST_SEED + 4;
const DIFFERENCES: usize = SECOND_SEED + 4;

/// The largest fixed point at which every difference fits an i32.
pub fn optimal_linear_fixed_point(values: &[f64]) -> f64 {
    let largest = match values {
        [] => return 0.0,
        [only] => *only,
        [first, second, ..] => {
            let mut largest = first.max(*second);
            for window in values.windows(3) {
                let extrapolated = window[1] + (window[1] - window[0]);
                let difference = window[2] - extrapolated;
                largest = largest.max((difference.abs() + 1.0).ceil());
            }
            largest
        }
    };
    (i32::MAX as f64 / largest).floor()
}

fn scaled(value: f64, fixed_point: f64) -> Result<i64, Error> {
    let scaled = value * fixed_point + 0.5;
    if !(scaled >= i64::MIN as f64 && scaled < i64::MAX as f64) {
        return Err(Error::Overflow);
    }
    Ok(scaled as i64)
}

/// The first two values are stored as u32, so must scale into its range.
fn seed(value: f64, fixed_point: f64) -> Result<u32, Error> {
    u32::try_from(scaled(value, fixed_point)?).map_err(|_| Error::Overflow)
}

pub fn encode_linear(values: &[f64], fixed_point: f64, out: &mut Vec<u8>) -> Result<(), Error> {
    all_or_nothing(out, |out| {
        write_fixed_point(fixed_point, out);
        if values.is_empty() {
            return Ok(());
        }
        let first = seed(values[0], fixed_point)?;
        out.extend_from_slice(&first.to_le_bytes());
        if values.len() == 1 {
            return Ok(());
        }
        let second = seed(values[1], fixed_point)?;
        out.extend_from_slice(&second.to_le_bytes());

        let mut before_previous: i64;
        let mut previous = i64::from(first);
        let mut current = i64::from(second);
        let mut ints = HalfByteWriter::new(out);
        for &value in &values[2..] {
            before_previous = previous;
            previous = current;
            current = scaled(value, fixed_point)?;
            let extrapolated = 2 * i128::from(previous) - i128::from(before_previous);
            let difference =
                i32::try_from(i128::from(current) - extrapolated).map_err(|_| Error::Overflow)?;
            ints.push_int(difference as u32);
        }
        ints.finish();
        Ok(())
    })
}

pub fn decode_linear<T: Float>(data: &[u8], out: &mut Vec<T>) -> Result<(), Error> {
    if data.len() == FIXED_POINT_BYTES {
        return Ok(());
    }
    let fixed_point = read_fixed_point(data)?;
    let seed = |at: usize| -> Result<i64, Error> {
        let bytes = data.get(at..at + 4).ok_or(Error::TooShort(data.len()))?;
        Ok(u32::from_le_bytes(bytes.try_into().expect("4 bytes")) as i64)
    };
    let mut before_previous: i64;
    let mut previous = seed(FIRST_SEED)?;
    out.push(T::from_f64(previous as f64 / fixed_point));
    if data.len() == SECOND_SEED {
        return Ok(());
    }
    let mut current = seed(SECOND_SEED)?;
    out.push(T::from_f64(current as f64 / fixed_point));

    let mut ints = HalfByteReader::new(&data[DIFFERENCES..]);
    while !ints.at_end() {
        before_previous = previous;
        previous = current;
        let difference = ints.next_int()? as i32;
        current = previous
            .checked_sub(before_previous)
            .and_then(|step| previous.checked_add(step))
            .and_then(|extrapolated| extrapolated.checked_add(difference as i64))
            .ok_or(Error::Corrupt)?;
        out.push(T::from_f64(current as f64 / fixed_point));
    }
    Ok(())
}
