#![allow(dead_code)]

pub mod decoded_by_reference;
pub mod encoded_by_reference;

use numpression::Codec;

pub const CODECS: [Codec; 3] = [Codec::Linear, Codec::Pic, Codec::Slof];

/// Deterministic xorshift, so failures reproduce.
pub struct Random(u64);

impl Random {
    pub fn new() -> Self {
        Self(0x2545_f491_4f6c_dd1d)
    }

    pub fn next_u64(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    pub fn byte(&mut self) -> u8 {
        (self.next_u64() >> 24) as u8
    }

    /// Uniform in `[0, 1)`.
    pub fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    pub fn below(&mut self, bound: usize) -> usize {
        (self.unit() * bound as f64) as usize
    }

    /// Sorted m/z-like values in `[100, 2000)`.
    pub fn mz(&mut self, count: usize) -> Vec<f64> {
        let mut values: Vec<f64> = (0..count).map(|_| 100.0 + self.unit() * 1900.0).collect();
        values.sort_by(f64::total_cmp);
        values
    }

    /// Intensity-like values spread over orders of magnitude in `[0, 1e7)`.
    pub fn intensities(&mut self, count: usize) -> Vec<f64> {
        (0..count)
            .map(|_| 10f64.powf(self.unit() * 7.0) - 1.0)
            .collect()
    }
}

pub fn values_for(codec: Codec, random: &mut Random, count: usize) -> Vec<f64> {
    match codec {
        Codec::Linear => random.mz(count),
        Codec::Pic | Codec::Slof => random.intensities(count),
    }
}
