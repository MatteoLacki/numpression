//! Decoding what was encoded stays within each codec's error bound.

mod common;

use common::{values_for, Random, CODECS};
use numpression::{optimal_linear_fixed_point, optimal_slof_fixed_point, Codec};

/// Float slack on top of the codec's own rounding.
const SLACK: f64 = 1e-9;

#[test]
fn round_trips_within_the_codec_bounds() {
    let mut random = Random::new();
    for _ in 0..500 {
        let count = random.below(400);
        for codec in CODECS {
            let values = values_for(codec, &mut random, count);
            let mut bytes = Vec::new();
            codec.encode(&values, &mut bytes).unwrap();
            let mut decoded: Vec<f64> = Vec::new();
            codec.decode(&bytes, &mut decoded).unwrap();
            assert_eq!(decoded.len(), values.len(), "{codec:?}");

            for (&value, &back) in values.iter().zip(&decoded) {
                let (error, bound) = match codec {
                    Codec::Linear => (
                        (back - value).abs(),
                        0.5 / optimal_linear_fixed_point(&values),
                    ),
                    Codec::Pic => ((back - value).abs(), 0.5),
                    Codec::Slof => (
                        ((back + 1.0).ln() - (value + 1.0).ln()).abs(),
                        0.5 / optimal_slof_fixed_point(&values),
                    ),
                };
                assert!(
                    error <= bound * (1.0 + SLACK) + SLACK,
                    "{codec:?}: {value} came back as {back}"
                );
            }
        }
    }
}
