//! Integers stored as variable runs of half-bytes, high half of each byte first.
//!
//! An integer is a head half-byte, then its remaining half-bytes lowest first. A head
//! `h` <= 8 says `h` leading zero half-bytes were dropped, `h` > 8 that `h - 8` leading
//! 0xf half-bytes were (at most 7, so the head never exceeds 0xf).

use crate::Error;

pub(crate) struct HalfByteReader<'a> {
    data: &'a [u8],
    byte: usize,
    in_low_half: bool,
}

impl<'a> HalfByteReader<'a> {
    pub(crate) fn new(data: &'a [u8]) -> Self {
        Self {
            data,
            byte: 0,
            in_low_half: false,
        }
    }

    /// No integer left: past the last byte, or at a last low half-byte of zero padding.
    pub(crate) fn at_end(&self) -> bool {
        let last = self.data.len().saturating_sub(1);
        self.byte >= self.data.len()
            || (self.byte == last && self.in_low_half && self.data[last] & 0xf == 0)
    }

    fn next_half_byte(&mut self) -> Result<u8, Error> {
        let byte = *self.data.get(self.byte).ok_or(Error::Corrupt)?;
        if self.in_low_half {
            self.byte += 1;
            self.in_low_half = false;
            Ok(byte & 0xf)
        } else {
            self.in_low_half = true;
            Ok(byte >> 4)
        }
    }

    pub(crate) fn next_int(&mut self) -> Result<u32, Error> {
        let head = self.next_half_byte()? as usize;
        let (dropped, mut value) = match head {
            0..=8 => (head, 0u32),
            _ => {
                let ones = head - 8;
                let leading_ones = (0..ones).fold(0u32, |v, i| v | (0xf000_0000 >> (4 * i)));
                (ones, leading_ones)
            }
        };
        if dropped == 8 {
            return Ok(value);
        }
        let pending_low_half = usize::from(!self.in_low_half);
        if self.byte + ((8 - dropped) - pending_low_half) / 2 >= self.data.len() {
            return Err(Error::Corrupt);
        }
        for i in dropped..8 {
            let half_byte = self.next_half_byte()? as u32;
            value |= half_byte << ((i - dropped) * 4);
        }
        Ok(value)
    }
}

pub(crate) struct HalfByteWriter<'a> {
    out: &'a mut Vec<u8>,
    pending_high_half: Option<u8>,
}

impl<'a> HalfByteWriter<'a> {
    pub(crate) fn new(out: &'a mut Vec<u8>) -> Self {
        Self {
            out,
            pending_high_half: None,
        }
    }

    fn push_half_byte(&mut self, half_byte: u8) {
        let half_byte = half_byte & 0xf;
        match self.pending_high_half.take() {
            None => self.pending_high_half = Some(half_byte),
            Some(high) => self.out.push(high << 4 | half_byte),
        }
    }

    pub(crate) fn push_int(&mut self, value: u32) {
        let (head, dropped) = match value & 0xf000_0000 {
            0 => {
                let zeros = value.leading_zeros() / 4;
                (zeros, zeros)
            }
            0xf000_0000 => {
                let ones = ((!value).leading_zeros() / 4).min(7);
                (ones + 8, ones)
            }
            _ => (0, 0),
        };
        self.push_half_byte(head as u8);
        for i in dropped..8 {
            self.push_half_byte((value >> (4 * (i - dropped))) as u8);
        }
    }

    /// Pads a last lone half-byte with a zero low half.
    pub(crate) fn finish(self) {
        if let Some(high) = self.pending_high_half {
            self.out.push(high << 4);
        }
    }
}
