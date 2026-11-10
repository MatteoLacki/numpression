# numpression

[MS-Numpress](https://github.com/ms-numpress/ms-numpress) encoding and decoding for mass
spectrometry data, in safe Rust with no dependencies.

| Codec | For | Error per value |
|---|---|---|
| `Linear` (linear prediction) | m/z, retention time | at most `0.5 / fixed_point` |
| `Pic` (positive integer) | intensities | at most 0.5 |
| `Slof` (short logged float) | intensities | at most `0.5 / fixed_point` in `ln(x + 1)` |

- Encoders produce the same bytes as the reference C++ (`MSNumpress.cpp`).
- Decoders bounds-check every read: a corrupt array returns an `Error`, never panics.
- `#![forbid(unsafe_code)]`.

How each codec works, byte by byte, with worked examples: [THEORY.md](THEORY.md).

## Usage

```toml
[dependencies]
numpression = "0.1"
```

```rust
use numpression::Codec;

let mz = [445.12, 446.5, 512.25, 1001.75];
let mut bytes = Vec::new();
Codec::Linear.encode(&mz, &mut bytes)?; // at the optimal fixed point

let mut decoded: Vec<f64> = Vec::new(); // or Vec<f32>
Codec::Linear.decode(&bytes, &mut decoded)?;
# Ok::<(), numpression::Error>(())
```

Every function appends to its output buffer; an encoder that fails (a value the codec
cannot hold, or NaN) leaves the buffer as it found it. To choose the fixed point
yourself, use `encode_linear` / `encode_slof` with `optimal_linear_fixed_point` /
`optimal_slof_fixed_point` as a starting point.

In mzML, the codecs are named by these binary-array compression CV terms:

| Codec | Alone | Followed by zlib |
|---|---|---|
| `Linear` | MS:1002312 | MS:1002746 |
| `Pic` | MS:1002313 | MS:1002747 |
| `Slof` | MS:1002314 | MS:1002748 |

## Tests

`cargo test` checks the encoders against bytes produced by the reference C++, the codecs
against [`numpress-rs`](https://crates.io/crates/numpress-rs) on random data, error bounds
on round trips, corrupt input, and fragments of 30 peptides in simulated DIA scans.

## DIA-NN interop (`interop/`, not part of the crate)

`interop/run_diann.sh <diann-linux> <fasta> <workdir> [python]` predicts a DIA-NN
library from the FASTA, simulates a DIA run from it (`make_synthetic_dia.py`, needs numpy,
pandas, pyarrow), re-encodes the run as 32-bit zlib, linear + pic and linear + slof
(`renumpress`), searches each with DIA-NN and compares the reports
(`compare_reports.py`). With DIA-NN 2.7.0 and a 200-protein FASTA, all 3000 spiked
precursors are identified in every encoding.

## License

Apache-2.0. Ported from ms-numpress `MSNumpress.cpp`, Copyright 2013 Johan Teleman,
Apache-2.0.
