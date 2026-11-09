//! Against bytes and values produced by the reference C++, plus corrupt input.

mod common;

use common::{decoded_by_reference, encoded_by_reference, Random, CODECS};
use numpression::{Codec, Error};

type Table<T> = &'static [(&'static str, &'static [T], &'static [u8])];

fn encoded_cases() -> impl Iterator<Item = (Codec, &'static str, Vec<f64>, &'static [u8])> {
    let tag = |codec, table: Table<u64>| {
        table.iter().map(move |&(name, bits, bytes)| {
            let values = bits.iter().map(|&b| f64::from_bits(b)).collect();
            (codec, name, values, bytes)
        })
    };
    tag(Codec::Linear, encoded_by_reference::LINEAR)
        .chain(tag(Codec::Pic, encoded_by_reference::PIC))
        .chain(tag(Codec::Slof, encoded_by_reference::SLOF))
}

fn decoded_cases() -> impl Iterator<Item = (Codec, &'static str, &'static [u8], &'static [u32])> {
    let tag = |codec, table: &'static [(&'static str, &'static [u8], &'static [u32])]| {
        table
            .iter()
            .map(move |&(name, bytes, bits)| (codec, name, bytes, bits))
    };
    tag(Codec::Linear, decoded_by_reference::LINEAR)
        .chain(tag(Codec::Pic, decoded_by_reference::PIC))
        .chain(tag(Codec::Slof, decoded_by_reference::SLOF))
}

fn decoded_bits(codec: Codec, bytes: &[u8]) -> Result<Vec<u32>, Error> {
    let mut out: Vec<f32> = Vec::new();
    codec.decode(bytes, &mut out)?;
    Ok(out.iter().map(|v| v.to_bits()).collect())
}

#[test]
fn encodes_exactly_as_the_reference() {
    for (codec, name, values, expected) in encoded_cases() {
        let mut out = Vec::new();
        codec.encode(&values, &mut out).unwrap();
        assert_eq!(out, expected, "{codec:?} {name}");
    }
}

#[test]
fn decodes_exactly_as_the_reference() {
    for (codec, name, bytes, expected) in decoded_cases() {
        assert_eq!(
            decoded_bits(codec, bytes).as_deref(),
            Ok(expected),
            "{name}"
        );
    }
}

#[test]
fn encoding_and_decoding_append() {
    for (codec, name, values, expected) in encoded_cases() {
        let mut bytes = vec![0xaa];
        codec.encode(&values, &mut bytes).unwrap();
        assert_eq!(
            (bytes[0], &bytes[1..]),
            (0xaa, expected),
            "{codec:?} {name}"
        );

        let mut decoded = vec![-1.0f64];
        codec.decode(expected, &mut decoded).unwrap();
        assert_eq!(decoded.len(), 1 + values.len(), "{codec:?} {name}");
        assert_eq!(decoded[0], -1.0, "{codec:?} {name}");
    }
}

#[test]
fn failed_encoding_leaves_the_output_alone() {
    let unencodable: [(Codec, &[f64]); 6] = [
        (Codec::Linear, &[1.0, f64::NAN]),
        (Codec::Linear, &[-1.0, 1.0]),
        (Codec::Pic, &[3.0, -1.0]),
        (Codec::Pic, &[3.0, 1e10]),
        (Codec::Slof, &[3.0, f64::NAN]),
        (Codec::Slof, &[-0.5]),
    ];
    for (codec, values) in unencodable {
        let mut out = vec![7, 8, 9];
        assert_eq!(
            codec.encode(values, &mut out),
            Err(Error::Overflow),
            "{codec:?} {values:?}"
        );
        assert_eq!(out, [7, 8, 9], "{codec:?} {values:?}");
    }
}

#[test]
fn truncated_arrays_fail_or_decode_a_prefix() {
    for (codec, name, bytes, expected) in decoded_cases() {
        for length in 0..bytes.len() {
            if let Ok(bits) = decoded_bits(codec, &bytes[..length]) {
                assert!(expected.starts_with(&bits), "{name} cut to {length} bytes");
            }
        }
    }
}

#[test]
fn garbage_never_panics() {
    let mut random = Random::new();
    for length in 0..200 {
        for _ in 0..50 {
            let bytes: Vec<u8> = (0..length).map(|_| random.byte()).collect();
            for codec in CODECS {
                let _ = decoded_bits(codec, &bytes);
            }
        }
    }
}

#[test]
fn odd_slof_payload_is_corrupt() {
    let (_, _, bytes, _) = decoded_cases().find(|case| case.1 == "slof_1").unwrap();
    assert_eq!(
        decoded_bits(Codec::Slof, &bytes[..bytes.len() - 1]),
        Err(Error::Corrupt)
    );
}
