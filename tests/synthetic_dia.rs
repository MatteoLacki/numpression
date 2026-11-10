//! 30 tryptic peptides in a simulated DIA run: every MS2 scan is encoded and decoded
//! with linear (m/z) plus pic or slof (intensity), and at each peptide's apex scan its
//! fragments must still be found at the right m/z with the right relative intensities.
//! Fragments are b/y ions with a fixed, made-up intensity pattern; noise peaks on top.

mod common;

use common::Random;
use numpression::{optimal_linear_fixed_point, Codec};

const PEPTIDES: [&str; 30] = [
    "LVNELTEFAK",
    "AEFVEVTK",
    "YLYEIAR",
    "HPYFYAPELLYYANK",
    "QTALVELLK",
    "DDPHACYSTVFDK",
    "SLHTLFGDELCK",
    "EYEATLEECCAK",
    "LGEYGFQNALIVR",
    "KVPQVSTPTLVEVSR",
    "VFDEFKPLVEEPQNLIK",
    "TCVADESHAGCEK",
    "ETYGDMADCCEK",
    "NECFLSHK",
    "FKDLGEEHFK",
    "LCVLHEK",
    "SHCIAEVEK",
    "ADDKETCFAEEGK",
    "AVMDDFAAFVEK",
    "MPCTEDYLSLILNR",
    "GLVLIAFSQYLQQCPFDEHVK",
    "RHPEYAVSVLLR",
    "DAIPENLPPLTADFAEDK",
    "TPVSEKVTK",
    "CCTESLVNR",
    "QNCDQFEK",
    "IGDYAGIK",
    "LAQGVSGAVR",
    "VVSTQTALA",
    "ELLESYIDGR",
];

const PROTON: f64 = 1.007_276_466_8;
const WATER: f64 = 18.010_564_684;
const CYCLE_S: f64 = 1.5;
const SCANS_PER_CYCLE: usize = 25;
const WINDOW_START: f64 = 400.0;
const WINDOW_WIDTH: f64 = 24.0;
const SIGMA_S: f64 = 4.0;
const NOISE_PEAKS: usize = 150;

/// Monoisotopic residue masses; C carries carbamidomethyl.
fn residue_mass(aa: u8) -> f64 {
    match aa {
        b'G' => 57.021_464,
        b'A' => 71.037_114,
        b'S' => 87.032_028,
        b'P' => 97.052_764,
        b'V' => 99.068_414,
        b'T' => 101.047_679,
        b'C' => 160.030_649,
        b'L' | b'I' => 113.084_064,
        b'N' => 114.042_927,
        b'D' => 115.026_943,
        b'Q' => 128.058_578,
        b'K' => 128.094_963,
        b'E' => 129.042_593,
        b'M' => 131.040_485,
        b'H' => 137.058_912,
        b'F' => 147.068_414,
        b'R' => 156.101_111,
        b'Y' => 163.063_329,
        b'W' => 186.079_313,
        _ => panic!("unknown residue {}", aa as char),
    }
}

struct Peptide {
    precursor_mz: f64,
    apex_s: f64,
    abundance: f64,
    fragment_mz: Vec<f64>,
    relative_intensity: Vec<f64>,
}

/// Singly charged b2.. and y2.. ions. y ions fall off from the precursor end, b ions
/// are weaker; a made-up but peptide-specific pattern.
fn peptide(sequence: &str, number: usize, random: &mut Random) -> Peptide {
    let residues = sequence.as_bytes();
    let mut masses = Vec::with_capacity(residues.len());
    for &aa in residues {
        masses.push(residue_mass(aa));
    }
    let total: f64 = masses.iter().sum();
    let mut fragment_mz = Vec::new();
    let mut relative_intensity = Vec::new();
    let mut b = 0.0;
    for cut in 1..residues.len() {
        b += masses[cut - 1];
        let y = total - b;
        let from_end = (residues.len() - cut) as f64;
        if cut >= 2 {
            fragment_mz.push(b + PROTON);
            relative_intensity.push(0.3 * (0.8f64).powf(cut as f64) + 0.05);
        }
        if from_end >= 2.0 {
            fragment_mz.push(y + WATER + PROTON);
            relative_intensity.push((0.88f64).powf(cut as f64) * (1.0 - 0.02 * from_end).max(0.2));
        }
    }
    Peptide {
        precursor_mz: lowest_charge_in_windows(total + WATER),
        apex_s: 30.0 + 10.0 * number as f64,
        abundance: 1e4 * 10f64.powf(3.0 * random.unit()),
        fragment_mz,
        relative_intensity,
    }
}

/// The precursor m/z at the lowest charge (2 or 3) that falls inside the windows.
fn lowest_charge_in_windows(neutral_mass: f64) -> f64 {
    for charge in 2..=3 {
        let mz = (neutral_mass + charge as f64 * PROTON) / charge as f64;
        if window_of(mz).is_some() {
            return mz;
        }
    }
    panic!("no charge puts mass {neutral_mass} inside the windows");
}

fn window_of(mz: f64) -> Option<usize> {
    let window = ((mz - WINDOW_START) / WINDOW_WIDTH).floor();
    if window >= 0.0 && (window as usize) < SCANS_PER_CYCLE {
        return Some(window as usize);
    }
    None
}

/// One MS2 scan: fragments of every peptide in `window`, scaled by its elution profile,
/// plus noise peaks; sorted by m/z.
fn scan(
    peptides: &[Peptide],
    window: usize,
    time_s: f64,
    random: &mut Random,
) -> (Vec<f64>, Vec<f64>) {
    let mut peaks: Vec<(f64, f64)> = Vec::new();
    for peptide in peptides {
        if window_of(peptide.precursor_mz) != Some(window) {
            continue;
        }
        let elution = (-0.5 * ((time_s - peptide.apex_s) / SIGMA_S).powi(2)).exp();
        if elution < 0.01 {
            continue;
        }
        for k in 0..peptide.fragment_mz.len() {
            let intensity = peptide.abundance * elution * peptide.relative_intensity[k];
            peaks.push((peptide.fragment_mz[k], intensity));
        }
    }
    for _ in 0..NOISE_PEAKS {
        peaks.push((
            150.0 + 1350.0 * random.unit(),
            100.0 + 2000.0 * random.unit(),
        ));
    }
    peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut mz = Vec::with_capacity(peaks.len());
    let mut intensity = Vec::with_capacity(peaks.len());
    for (m, i) in peaks {
        mz.push(m);
        intensity.push(i);
    }
    (mz, intensity)
}

fn round_trip(codec: Codec, values: &[f64]) -> Vec<f64> {
    let mut bytes = Vec::new();
    codec.encode(values, &mut bytes).unwrap();
    let mut decoded = Vec::new();
    codec.decode(&bytes, &mut decoded).unwrap();
    assert_eq!(decoded.len(), values.len());
    decoded
}

fn pearson(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let mean_x = x.iter().sum::<f64>() / n;
    let mean_y = y.iter().sum::<f64>() / n;
    let (mut xy, mut xx, mut yy) = (0.0, 0.0, 0.0);
    for k in 0..x.len() {
        xy += (x[k] - mean_x) * (y[k] - mean_y);
        xx += (x[k] - mean_x) * (x[k] - mean_x);
        yy += (y[k] - mean_y) * (y[k] - mean_y);
    }
    xy / (xx * yy).sqrt()
}

/// Index of the decoded peak nearest `target`.
fn nearest(mz: &[f64], target: f64) -> usize {
    let at = mz.partition_point(|&m| m < target);
    let mut best = at.min(mz.len() - 1);
    if at > 0 && (mz[at - 1] - target).abs() < (mz[best] - target).abs() {
        best = at - 1;
    }
    best
}

#[test]
fn peptides_survive_numpress_at_their_apex() {
    let mut random = Random::new();
    let mut peptides = Vec::new();
    for (number, sequence) in PEPTIDES.iter().enumerate() {
        peptides.push(peptide(sequence, number, &mut random));
    }
    let cycles = ((peptides.last().unwrap().apex_s + 30.0) / CYCLE_S) as usize;
    let scan_dt = CYCLE_S / SCANS_PER_CYCLE as f64;

    for intensity_codec in [Codec::Pic, Codec::Slof] {
        let mut checked = 0;
        for (number, peptide) in peptides.iter().enumerate() {
            let window = window_of(peptide.precursor_mz)
                .unwrap_or_else(|| panic!("{} outside the windows", PEPTIDES[number]));
            let apex_cycle = ((peptide.apex_s - (window as f64 + 0.5) * scan_dt) / CYCLE_S).round();
            assert!(apex_cycle >= 0.0 && (apex_cycle as usize) < cycles);
            let time_s = apex_cycle * CYCLE_S + (window as f64 + 0.5) * scan_dt;
            let (mz, intensity) = scan(&peptides, window, time_s, &mut random);

            let decoded_mz = round_trip(Codec::Linear, &mz);
            let decoded_intensity = round_trip(intensity_codec, &intensity);
            let mz_bound = 0.5 / optimal_linear_fixed_point(&mz) + 1e-9;

            let mut found = Vec::with_capacity(peptide.fragment_mz.len());
            for &fragment in &peptide.fragment_mz {
                let at = nearest(&decoded_mz, fragment);
                assert!(
                    (decoded_mz[at] - fragment).abs() <= mz_bound,
                    "{intensity_codec:?} {}: fragment {fragment} came back as {}",
                    PEPTIDES[number],
                    decoded_mz[at]
                );
                found.push(decoded_intensity[at]);
            }
            let r = pearson(&found, &peptide.relative_intensity);
            assert!(
                r > 0.99,
                "{intensity_codec:?} {}: r = {r}",
                PEPTIDES[number]
            );
            checked += 1;
        }
        assert_eq!(checked, PEPTIDES.len());
    }
}
