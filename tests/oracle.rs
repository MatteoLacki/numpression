//! Against `numpress-rs`, an independent port of the reference, on valid input only:
//! it panics on some corrupt input.

mod common;

use common::{values_for, Random, CODECS};
use numpression::{optimal_linear_fixed_point, optimal_slof_fixed_point, Codec};

fn oracle_encode(codec: Codec, values: &[f64]) -> Vec<u8> {
    let mut out = Vec::new();
    match codec {
        Codec::Linear => {
            out = numpress_rs::numpress_compress(values, numpress_rs::optimal_scaling(values))
                .unwrap();
        }
        Codec::Pic => {
            numpress_rs::encode_pic(values, &mut out).unwrap();
        }
        Codec::Slof => {
            let fixed_point = numpress_rs::optimal_slof_fixed_point(values);
            numpress_rs::encode_slof(values, &mut out, fixed_point).unwrap();
        }
    }
    out
}

fn oracle_decode(codec: Codec, bytes: &[u8]) -> Vec<f64> {
    let mut out = Vec::new();
    match codec {
        Codec::Linear => numpress_rs::decode_linear(bytes, &mut out),
        Codec::Pic => numpress_rs::decode_pic(bytes, &mut out),
        Codec::Slof => numpress_rs::decode_slof(bytes, &mut out),
    }
    .unwrap();
    out
}

#[test]
fn fixed_points_match_the_oracle() {
    let mut random = Random::new();
    for _ in 0..1000 {
        let count = random.below(400);
        let mz = random.mz(count);
        let intensities = random.intensities(count);
        assert_eq!(
            optimal_linear_fixed_point(&mz).to_bits(),
            numpress_rs::optimal_scaling(&mz).to_bits()
        );
        if count > 0 {
            assert_eq!(
                optimal_slof_fixed_point(&intensities).to_bits(),
                numpress_rs::optimal_slof_fixed_point(&intensities).to_bits()
            );
        }
    }
}

#[test]
fn bytes_and_values_match_the_oracle() {
    let mut random = Random::new();
    for _ in 0..1000 {
        let count = random.below(400);
        for codec in CODECS {
            // The oracle gives an empty slof array fixed point 65535, the reference 0
            // (covered by the reference vectors).
            if codec == Codec::Slof && count == 0 {
                continue;
            }
            let values = values_for(codec, &mut random, count);
            let mut ours = Vec::new();
            codec.encode(&values, &mut ours).unwrap();
            assert_eq!(ours, oracle_encode(codec, &values), "{codec:?} encode");

            let mut decoded: Vec<f64> = Vec::new();
            codec.decode(&ours, &mut decoded).unwrap();
            let oracle: Vec<u64> = oracle_decode(codec, &ours)
                .iter()
                .map(|v| v.to_bits())
                .collect();
            let decoded: Vec<u64> = decoded.iter().map(|v| v.to_bits()).collect();
            assert_eq!(decoded, oracle, "{codec:?} decode");
        }
    }
}
