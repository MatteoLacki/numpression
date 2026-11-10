//! Rewrite every binary array of an mzML with another encoding and rebuild the spectrum
//! index. Usage: renumpress <in.mzML> <out.mzML> <zlib32|pic|slof>
//!
//! zlib32: 32-bit float + zlib (DIA-NN's recommended msconvert output).
//! pic / slof: numpression linear for m/z, pic / slof for intensity, then zlib.
//! Expects one element per line, as pmsms2mzml and make_synthetic_dia.py write.
//! Index offsets point at the start of the `<spectrum` line.

use base64::{engine::general_purpose::STANDARD, Engine};
use flate2::{read::ZlibDecoder, write::ZlibEncoder, Compression};
use numpression::Codec;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};

const ZLIB: &str = "MS:1000574";
const NO_COMPRESSION: &str = "MS:1000576";
const FLOAT_32: &str = "MS:1000521";
const FLOAT_64: &str = "MS:1000523";
const LINEAR_ZLIB: &str = "MS:1002746";
const PIC_ZLIB: &str = "MS:1002747";
const SLOF_ZLIB: &str = "MS:1002748";
const LINEAR: &str = "MS:1002312";
const PIC: &str = "MS:1002313";
const SLOF: &str = "MS:1002314";
const MZ_ARRAY: &str = "MS:1000514";
const INTENSITY_ARRAY: &str = "MS:1000515";

#[derive(Clone, Copy, PartialEq)]
enum Mode {
    Zlib32,
    Pic,
    Slof,
}

#[derive(Default)]
struct ArrayIn {
    header_indent: String,
    cv_indent: String,
    kept: Vec<String>,
    zlib: bool,
    numpress: Option<Codec>,
    f64: bool,
    is_mz: bool,
    is_intensity: bool,
}

#[derive(Default)]
struct Stats {
    arrays: u64,
    values: u64,
    same_codec: u64,
    identical_numpress_bytes: u64,
}

fn accession(line: &str) -> Option<&str> {
    let start = line.find("accession=\"")? + "accession=\"".len();
    let end = start + line[start..].find('"')?;
    Some(&line[start..end])
}

fn indent(line: &str) -> String {
    line[..line.len() - line.trim_start().len()].to_string()
}

/// The decoded values and the bytes they were decoded from (after zlib).
fn decode(array: &ArrayIn, base64: &str) -> (Vec<f64>, Vec<u8>) {
    let raw = STANDARD.decode(base64).expect("base64");
    let bytes = if array.zlib {
        let mut inflated = Vec::new();
        ZlibDecoder::new(raw.as_slice())
            .read_to_end(&mut inflated)
            .expect("zlib");
        inflated
    } else {
        raw
    };
    let mut values = Vec::new();
    match array.numpress {
        Some(codec) => codec.decode(&bytes, &mut values).expect("numpress"),
        None if array.f64 => {
            for chunk in bytes.chunks_exact(8) {
                values.push(f64::from_le_bytes(chunk.try_into().unwrap()));
            }
        }
        None => {
            for chunk in bytes.chunks_exact(4) {
                values.push(f32::from_le_bytes(chunk.try_into().unwrap()) as f64);
            }
        }
    }
    (values, bytes)
}

fn zlib(bytes: &[u8]) -> Vec<u8> {
    let mut encoder = ZlibEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(bytes).unwrap();
    encoder.finish().unwrap()
}

/// The output codec of an array and its CV accessions, in pmsms2mzml's order.
fn target(mode: Mode, array: &ArrayIn) -> (Option<Codec>, &'static [&'static str]) {
    match mode {
        Mode::Zlib32 => (None, &[FLOAT_32, ZLIB]),
        _ if array.is_mz => (Some(Codec::Linear), &[LINEAR_ZLIB, ZLIB]),
        Mode::Pic => (Some(Codec::Pic), &[PIC_ZLIB, ZLIB]),
        Mode::Slof => (Some(Codec::Slof), &[SLOF_ZLIB, ZLIB]),
    }
}

fn cv_name(accession: &str) -> &'static str {
    match accession {
        FLOAT_32 => "32-bit float",
        ZLIB => "zlib compression",
        LINEAR_ZLIB => "MS-Numpress linear prediction compression",
        PIC_ZLIB => "MS-Numpress positive integer compression",
        SLOF_ZLIB => "MS-Numpress short logged float compression",
        _ => unreachable!(),
    }
}

fn emit(out: &mut BufWriter<File>, written: &mut u64, line: &str) {
    out.write_all(line.as_bytes()).unwrap();
    out.write_all(b"\n").unwrap();
    *written += line.len() as u64 + 1;
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 4 {
        eprintln!("usage: renumpress <in.mzML> <out.mzML> <zlib32|pic|slof>");
        std::process::exit(2);
    }
    let mode = match args[3].as_str() {
        "zlib32" => Mode::Zlib32,
        "pic" => Mode::Pic,
        "slof" => Mode::Slof,
        other => panic!("unknown mode {other}"),
    };
    let input = BufReader::with_capacity(1 << 24, File::open(&args[1]).unwrap());
    let mut out = BufWriter::with_capacity(1 << 24, File::create(&args[2]).unwrap());

    let mut written: u64 = 0;
    let mut spectrum_offsets: Vec<u64> = Vec::new();
    let mut index_list_offset: u64 = 0;
    let mut in_spectrum_index = false;
    let mut next_index_entry = 0usize;
    let mut array: Option<ArrayIn> = None;
    let mut stats = Stats::default();

    for line in input.lines() {
        let line = line.unwrap();
        let trimmed = line.trim_start();

        if trimmed.starts_with("<binaryDataArray ") {
            array = Some(ArrayIn {
                header_indent: indent(&line),
                ..Default::default()
            });
            continue;
        }
        if let Some(current) = array.as_mut() {
            if trimmed.starts_with("<cvParam") {
                current.cv_indent = indent(&line);
                match accession(trimmed).unwrap_or("") {
                    ZLIB => current.zlib = true,
                    NO_COMPRESSION | FLOAT_32 => {}
                    FLOAT_64 => current.f64 = true,
                    LINEAR => current.numpress = Some(Codec::Linear),
                    PIC => current.numpress = Some(Codec::Pic),
                    SLOF => current.numpress = Some(Codec::Slof),
                    LINEAR_ZLIB => {
                        current.numpress = Some(Codec::Linear);
                        current.zlib = true;
                    }
                    PIC_ZLIB => {
                        current.numpress = Some(Codec::Pic);
                        current.zlib = true;
                    }
                    SLOF_ZLIB => {
                        current.numpress = Some(Codec::Slof);
                        current.zlib = true;
                    }
                    MZ_ARRAY => {
                        current.is_mz = true;
                        current.kept.push(line.clone());
                    }
                    INTENSITY_ARRAY => {
                        current.is_intensity = true;
                        current.kept.push(line.clone());
                    }
                    _ => current.kept.push(line.clone()),
                }
                continue;
            }
            if trimmed.starts_with("<binary>") || trimmed.starts_with("<binary/>") {
                assert!(current.is_mz || current.is_intensity, "unknown array type");
                let payload = trimmed
                    .strip_prefix("<binary>")
                    .and_then(|s| s.strip_suffix("</binary>"))
                    .unwrap_or("");
                let (values, input_bytes) = decode(current, payload);
                let (codec, accessions) = target(mode, current);
                let encoded = match codec {
                    Some(codec) => {
                        let mut bytes = Vec::new();
                        codec.encode(&values, &mut bytes).expect("encode");
                        if current.numpress == Some(codec) {
                            stats.same_codec += 1;
                            stats.identical_numpress_bytes += u64::from(bytes == input_bytes);
                        }
                        zlib(&bytes)
                    }
                    None => {
                        let mut floats = Vec::with_capacity(4 * values.len());
                        for &value in &values {
                            floats.extend_from_slice(&(value as f32).to_le_bytes());
                        }
                        zlib(&floats)
                    }
                };
                let base64 = STANDARD.encode(&encoded);
                stats.arrays += 1;
                stats.values += values.len() as u64;

                let header = format!(
                    "{}<binaryDataArray encodedLength=\"{}\">",
                    current.header_indent,
                    base64.len()
                );
                emit(&mut out, &mut written, &header);
                for accession in accessions {
                    let cv = format!(
                        "{}<cvParam cvRef=\"MS\" accession=\"{accession}\" name=\"{}\" value=\"\"/>",
                        current.cv_indent,
                        cv_name(accession)
                    );
                    emit(&mut out, &mut written, &cv);
                }
                for kept in &current.kept {
                    emit(&mut out, &mut written, kept);
                }
                let binary = format!("{}<binary>{base64}</binary>", indent(&line));
                emit(&mut out, &mut written, &binary);
                array = None;
                continue;
            }
            panic!("unexpected line inside binaryDataArray: {trimmed}");
        }

        if trimmed.starts_with("<spectrum ") {
            spectrum_offsets.push(written);
        } else if trimmed.starts_with("<indexList ") {
            index_list_offset = written;
        } else if trimmed.starts_with("<index name=\"spectrum\"") {
            in_spectrum_index = true;
        } else if trimmed.starts_with("</index>") {
            in_spectrum_index = false;
        } else if in_spectrum_index && trimmed.starts_with("<offset ") {
            let close = line.find("\">").unwrap() + 2;
            let rewritten = format!(
                "{}{}</offset>",
                &line[..close],
                spectrum_offsets[next_index_entry]
            );
            next_index_entry += 1;
            emit(&mut out, &mut written, &rewritten);
            continue;
        } else if trimmed.starts_with("<indexListOffset>") {
            let rewritten = format!(
                "{}<indexListOffset>{index_list_offset}</indexListOffset>",
                indent(&line)
            );
            emit(&mut out, &mut written, &rewritten);
            continue;
        }
        emit(&mut out, &mut written, &line);
    }
    out.flush().unwrap();
    assert_eq!(
        next_index_entry,
        spectrum_offsets.len(),
        "index entries != spectra"
    );
    eprintln!(
        "spectra {} arrays {} values {} | arrays already in the target numpress codec {}, re-encoded byte-identical {}",
        spectrum_offsets.len(),
        stats.arrays,
        stats.values,
        stats.same_codec,
        stats.identical_numpress_bytes
    );
}
