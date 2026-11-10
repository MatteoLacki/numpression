# How the MS-Numpress codecs work

MS-Numpress (Teleman et al., *Mol Cell Proteomics* 2014) compresses the two kinds of
arrays a mass spectrum is made of: smooth, increasing coordinates (m/z, retention time)
and noisy, positive magnitudes (intensities). All three codecs are lossy in a controlled
way: each turns floats into small integers, then stores the integers in as few bytes as
their size allows. This document walks through the byte layout of each, the error it
makes, and how to call it from this crate. The code blocks are compiled and run by
`cargo test`, so the bytes shown are the bytes the crate produces.

## Two building blocks

### The fixed point

Linear and slof scale every value by a *fixed point* `f` and round to an integer:
`i = trunc(x·f + 0.5)`, which is rounding to nearest for non-negative `x·f`. Decoding
divides back, `x' = i / f`, so the error is at most `0.5 / f`. A larger `f` means smaller
error but larger integers. Each codec has an *optimal* fixed point: the largest `f` at
which its integers still fit their storage.

`f` itself is stored in the first 8 bytes of the array as a **big-endian** f64 (the
reference writes it that way on every platform). Everything after it is little-endian.

### Half-byte integers

Linear and pic store a sequence of 32-bit integers in variable-length runs of half-bytes
(4 bits, "nibbles"), so small numbers cost little. An integer is written as a *head*
nibble followed by some of its own nibbles, **lowest first**:

- If the top nibbles are zeros: count them (`n`, 0 to 8). Head = `n`, then the remaining
  `8 − n` nibbles.
- If the top nibbles are `0xf` (a small negative number in two's complement): count them,
  capped at 7 (`n`). Head = `8 + n`, then the remaining `8 − n` nibbles.
- Otherwise: head = `0`, then all 8 nibbles.

| Value | Hex | Head | Nibbles after it | Total nibbles |
|---|---|---|---|---|
| 0 | `00000000` | 8 | none | 1 |
| 5 | `00000005` | 7 | 5 | 2 |
| 300 | `0000012c` | 5 | c, 2, 1 | 4 |
| −3 | `fffffffd` | 15 | d | 2 |
| −1 | `ffffffff` | 15 | f | 2 |
| 0x12345678 | `12345678` | 0 | 8, 7, 6, 5, 4, 3, 2, 1 | 9 |

The nibbles of consecutive integers are packed two per byte, **high half first**. A last
lone nibble gets a zero low half as padding. The decoder recognises the padding because
it is the low half of the final byte and is zero; a real integer can never start there
and end in the same byte, because a head of 0 is followed by 8 more nibbles.

## Positive integer (pic): intensities, absolute error ≤ 0.5

Each value is rounded to an integer and written as a half-byte integer. There is no
header and no fixed point. Values must round into `0..=i32::MAX`.

Ion counts are often small integers, so they cost 2 to 4 nibbles (1 to 2 bytes) each,
against 4 bytes for an f32. The error is absolute, 0.5, so pic suits detector counts and
is a poor choice for intensities scaled to small fractions.

```rust
use numpression::{decode_pic, encode_pic};

let intensities = [0.0, 1.0, 15.0, 300.0];
let mut bytes = Vec::new();
encode_pic(&intensities, &mut bytes)?;
// nibbles 8 | 7 1 | 7 f | 5 c 2 1, packed two per byte, the last one padded
assert_eq!(bytes, [0x87, 0x17, 0xf5, 0xc2, 0x10]);

let mut decoded: Vec<f64> = Vec::new();
decode_pic(&bytes, &mut decoded)?;
assert_eq!(decoded, intensities);
# Ok::<(), numpression::Error>(())
```

## Linear prediction (linear): m/z and retention time, error ≤ 0.5 / f

Sorted m/z values change slowly: each one is close to the straight line through the two
before it. Linear prediction stores only how far each value misses that line.

1. Scale and round every value: `i_k = trunc(x_k · f + 0.5)`.
2. Store the first two integers in full, as little-endian u32.
3. For every later value, predict `p_k = 2·i_(k−1) − i_(k−2)` and store the residual
   `r_k = i_k − p_k` as a half-byte integer (its two's-complement bits as a u32).

Layout: `f` (8 bytes, big-endian) | `i_0` (4 bytes) | `i_1` (4 bytes) | residuals. Decoding
runs the same recurrence forwards: `i_k = 2·i_(k−1) − i_(k−2) + r_k`, then `x_k = i_k / f`.
The integers are reconstructed exactly, so the only error is the initial rounding,
`0.5 / f` per value.

```rust
use numpression::{decode_linear, encode_linear};

let mz = [100.0, 101.5, 103.25];
let mut bytes = Vec::new();
encode_linear(&mz, 100.0, &mut bytes)?;
// i = 10000, 10150, 10325; prediction 2·10150 − 10000 = 10300, residual 25 = 0x19
assert_eq!(
    bytes,
    [
        0x40, 0x59, 0, 0, 0, 0, 0, 0, // f = 100.0, big-endian
        0x10, 0x27, 0, 0,             // i_0 = 10000
        0xa6, 0x27, 0, 0,             // i_1 = 10150
        0x69, 0x10,                   // residual: head 6, nibbles 9 1, padding
    ]
);

let mut decoded: Vec<f64> = Vec::new();
decode_linear(&bytes, &mut decoded)?;
assert_eq!(decoded, mz);
# Ok::<(), numpression::Error>(())
```

**Optimal fixed point.** `optimal_linear_fixed_point` returns
`floor((2³¹ − 1) / M)`, where `M` is the largest of `x_0`, `x_1` and, over all later
values, `ceil(|x_k − (2·x_(k−1) − x_(k−2))| + 1)`. So the seeds and every residual fit a
signed 32-bit integer. For a spectrum with peaks up to m/z 2000 and gaps of a few
hundred Th, `f` lands around 10⁶ to 10⁷: an error of 10⁻⁷ to 10⁻⁶ Th, well below
instrument accuracy. The more regular the spacing, the smaller the residuals and the
fewer nibbles each value takes. Profile data and dense centroided spectra compress best.

Limits: the first two values must scale into `0..=u32::MAX` (they are stored as u32), and
every residual must fit an i32. Otherwise `encode_linear` returns `Error::Overflow`. A
single negative value still encodes, with a negative fixed point.

## Short logged float (slof): intensities, relative error ≈ 0.5 / f

Intensities span many orders of magnitude, and what matters is their relative size.
Slof stores each one as a fixed-point logarithm in exactly 2 bytes:

`u = trunc(ln(x + 1) · f + 0.5)` as a little-endian u16; decoding gives
`x' = exp(u / f) − 1`.

Layout: `f` (8 bytes, big-endian) | one u16 per value. The error is at most `0.5 / f` in
`ln(x + 1)`, so about `0.5 / f` *relative* error in `x + 1`.

`optimal_slof_fixed_point` returns `floor(65535 / max(1, max ln(x + 1)))`, the largest
`f` at which every `u` fits a u16. For intensities up to 10⁷, `ln(10⁷ + 1) ≈ 16.1`, so
`f = 4070` and the relative error is about 1.2·10⁻⁴.

```rust
use numpression::{decode_slof, encode_slof};

let intensities = [0.0, 1.0, 1000.0];
let mut bytes = Vec::new();
encode_slof(&intensities, 1000.0, &mut bytes)?;
// u = 0, round(ln 2 · 1000) = 693, round(ln 1001 · 1000) = 6909
assert_eq!(
    bytes,
    [0x40, 0x8f, 0x40, 0, 0, 0, 0, 0, 0x00, 0x00, 0xb5, 0x02, 0xfd, 0x1a]
);

let mut decoded: Vec<f64> = Vec::new();
decode_slof(&bytes, &mut decoded)?;
for (x, back) in intensities.iter().zip(&decoded) {
    assert!(((back + 1.0).ln() - (x + 1.0).ln()).abs() <= 0.5 / 1000.0);
}
# Ok::<(), numpression::Error>(())
```

Values must satisfy `ln(x + 1) · f ≤ 65535`; negative values above −1 encode, but only
while the rounded logarithm stays representable.

## Choosing a codec

| Array | Codec | Bytes per value (typical) | Error |
|---|---|---|---|
| m/z, retention time | `Linear` | 1–3 | absolute, ≤ 0.5 / f |
| integer-like counts | `Pic` | 1–2 | absolute, ≤ 0.5 |
| intensities over many decades | `Slof` | 2 | relative, ≈ 0.5 / f |

All three are usually followed by zlib in mzML. Numpress leaves a lot of
low-entropy structure (repeated heads, small residuals), which zlib squeezes further.

## Calling it from Rust

`Codec` picks a codec at runtime (as when reading the compression CV term of an mzML
array) and always encodes at the optimal fixed point:

```rust
use numpression::Codec;

let mz = [445.12, 446.5, 512.25, 1001.75];
let mut bytes = Vec::new();
Codec::Linear.encode(&mz, &mut bytes)?;

let mut decoded: Vec<f32> = Vec::new(); // f32 or f64
Codec::Linear.decode(&bytes, &mut decoded)?;
assert_eq!(decoded.len(), mz.len());
# Ok::<(), numpression::Error>(())
```

The free functions (`encode_linear`, `encode_pic`, `encode_slof`, their `decode_`
counterparts, `optimal_linear_fixed_point`, `optimal_slof_fixed_point`) give control over
the fixed point, for example to fix it across a whole file, or to cap the error:
`encode_linear(&mz, 1e6, &mut out)` guarantees at most 5·10⁻⁷ Th per value, or returns
`Error::Overflow` if the residuals do not fit at that resolution.

Conventions shared by every function:

- **Appending.** Output goes to the end of the `Vec` passed in, so several arrays can be
  written into one buffer. An encoder that fails truncates the buffer back to where it
  started.
- **Values are computed in f64.** Decoding into `Vec<f32>` rounds once, at the end.
- **Errors.** `Error::TooShort` for an array shorter than its header, `Error::Corrupt` for
  bytes that are not a valid encoding (truncated integers, a slof payload of odd length, a
  linear sequence that overflows), `Error::Overflow` for values an encoder cannot hold
  (including NaN). No input makes a decoder panic.

### In an mzML reader or writer

An mzML `<binary>` element holds base64 text. When writing: encode with numpress, then
zlib if wanted, then base64, and name the result with the matching CV term (for example
MS:1002746, "MS-Numpress linear prediction compression followed by zlib compression").
When reading, undo it in reverse: base64, inflate if the term says zlib, then numpress.
Numpress arrays carry no 32- or 64-bit float CV term; the codec defines the type.
The CV terms per codec are listed in the [README](README.md).
