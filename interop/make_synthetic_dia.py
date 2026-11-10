"""A synthetic DIA run (indexed mzML, 32-bit float + zlib) from a DIA-NN library parquet.

Usage: make_synthetic_dia.py <lib.parquet> <out.mzML> <spiked.tsv> [--spiked N]

N random target precursors elute as Gaussians at their predicted RT (mapped onto a
30 min gradient); MS2 windows carry their library fragments, MS1 their isotopes; m/z
jitter, intensity noise and random noise peaks on top. Needs numpy, pandas, pyarrow.
"""

import argparse
import base64
import zlib

import numpy as np
import pandas as pd

GRADIENT_MIN = 30.0
CYCLE_S = 1.5
WINDOWS = [(400.0 + 24.0 * i, 424.0 + 24.0 * i) for i in range(25)]
SIGMA_S = 4.0
PPM = 2e-6
NOISE_MS2, NOISE_MS1 = 150, 300

parser = argparse.ArgumentParser()
parser.add_argument("library")
parser.add_argument("mzml")
parser.add_argument("spiked_tsv")
parser.add_argument("--spiked", type=int, default=3000)
parser.add_argument("--seed", type=int, default=7)
args = parser.parse_args()
rng = np.random.default_rng(args.seed)

lib = pd.read_parquet(args.library)
lib = lib[lib["Decoy"] == 0]
precursors = lib.drop_duplicates("Precursor.Id")[
    ["Precursor.Id", "Precursor.Mz", "Precursor.Charge", "RT"]
]
lo, hi = np.percentile(precursors["RT"], [0.5, 99.5])
in_range = precursors[
    precursors["Precursor.Mz"].between(WINDOWS[0][0], WINDOWS[-1][1])
    & precursors["Precursor.Charge"].between(2, 4)
    & precursors["RT"].between(lo, hi)
]
spiked = in_range.sample(args.spiked, random_state=args.seed).reset_index(drop=True)
spiked["rt_s"] = 60.0 * (1.0 + (GRADIENT_MIN - 2.0) * (spiked["RT"] - lo) / (hi - lo))
spiked["abundance"] = np.exp(rng.normal(np.log(3e5), 1.0, len(spiked)))
spiked["window"] = ((spiked["Precursor.Mz"] - WINDOWS[0][0]) // 24).astype(int)
spiked[["Precursor.Id", "Precursor.Mz", "Precursor.Charge", "rt_s", "abundance"]].to_csv(
    args.spiked_tsv, sep="\t", index=False
)

fragments = {
    pid: (g["Product.Mz"].to_numpy(), g["Relative.Intensity"].to_numpy())
    for pid, g in lib[lib["Precursor.Id"].isin(spiked["Precursor.Id"])].groupby("Precursor.Id")
}
pids = spiked["Precursor.Id"].to_numpy()
pmz = spiked["Precursor.Mz"].to_numpy()
pz = spiked["Precursor.Charge"].to_numpy()
prt = spiked["rt_s"].to_numpy()
pab = spiked["abundance"].to_numpy()
pwin = spiked["window"].to_numpy()


def profile(t):
    return np.exp(-0.5 * ((t - prt) / SIGMA_S) ** 2)


def jitter(mz):
    return mz * (1.0 + rng.normal(0.0, PPM, len(mz)))


def noisy(intensity):
    return intensity * np.clip(1.0 + rng.normal(0.0, 0.1, len(intensity)), 0.05, None)


def noise_peaks(count, lo_mz, hi_mz):
    return rng.uniform(lo_mz, hi_mz, count), np.exp(rng.normal(np.log(2e3), 1.0, count))


def ms1(t):
    p = profile(t)
    live = np.nonzero(p > 0.01)[0]
    mzs, ints = [], []
    for i in live:
        lam = pmz[i] * pz[i] / 1800.0
        ratios = np.array([1.0, lam, lam * lam / 2.0])
        mzs.append(pmz[i] + np.arange(3) * 1.003355 / pz[i])
        ints.append(pab[i] * p[i] * ratios / ratios.sum())
    nmz, nint = noise_peaks(NOISE_MS1, 350.0, 1250.0)
    return jitter(np.concatenate(mzs + [nmz])), noisy(np.concatenate(ints + [nint]))


def ms2(t, window):
    p = profile(t)
    live = np.nonzero((p > 0.01) & (pwin == window))[0]
    mzs, ints = [], []
    for i in live:
        fmz, frel = fragments[pids[i]]
        mzs.append(fmz)
        ints.append(pab[i] * p[i] * frel)
    nmz, nint = noise_peaks(NOISE_MS2, 150.0, 1500.0)
    return jitter(np.concatenate(mzs + [nmz])), noisy(np.concatenate(ints + [nint]))


def array_xml(values, accession, name, unit):
    data = base64.b64encode(zlib.compress(values.astype("<f4").tobytes())).decode()
    return (
        f'            <binaryDataArray encodedLength="{len(data)}">\n'
        '              <cvParam cvRef="MS" accession="MS:1000521" name="32-bit float" value=""/>\n'
        '              <cvParam cvRef="MS" accession="MS:1000574" name="zlib compression" value=""/>\n'
        f'              <cvParam cvRef="MS" accession="{accession}" name="{name}" value="" {unit}/>\n'
        f"              <binary>{data}</binary>\n"
        "            </binaryDataArray>\n"
    )


def spectrum_xml(index, t, mz, inten, window):
    order = np.argsort(mz)
    mz, inten = mz[order], inten[order]
    level = 1 if window is None else 2
    kind = (
        '<cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum" value=""/>'
        if level == 1
        else '<cvParam cvRef="MS" accession="MS:1000580" name="MSn spectrum" value=""/>'
    )
    s = [
        f'        <spectrum index="{index}" id="scan={index + 1}" defaultArrayLength="{len(mz)}">\n',
        f"          {kind}\n",
        f'          <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="{level}"/>\n',
        '          <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" value=""/>\n',
        '          <cvParam cvRef="MS" accession="MS:1000130" name="positive scan" value=""/>\n',
        '          <scanList count="1">\n',
        '            <cvParam cvRef="MS" accession="MS:1000795" name="no combination" value=""/>\n',
        "            <scan>\n",
        f'              <cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="{t / 60.0:.6f}" unitCvRef="UO" unitAccession="UO:0000031" unitName="minute"/>\n',
        "            </scan>\n",
        "          </scanList>\n",
    ]
    if level == 2:
        lo_mz, hi_mz = WINDOWS[window]
        centre = (lo_mz + hi_mz) / 2.0
        s += [
            '          <precursorList count="1">\n',
            "            <precursor>\n",
            "              <isolationWindow>\n",
            f'                <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="{centre}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>\n',
            f'                <cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="{centre - lo_mz}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>\n',
            f'                <cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="{hi_mz - centre}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>\n',
            "              </isolationWindow>\n",
            '              <selectedIonList count="1">\n',
            "                <selectedIon>\n",
            f'                  <cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="{centre}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>\n',
            "                </selectedIon>\n",
            "              </selectedIonList>\n",
            "              <activation>\n",
            '                <cvParam cvRef="MS" accession="MS:1000422" name="beam-type collision-induced dissociation" value=""/>\n',
            '                <cvParam cvRef="MS" accession="MS:1000045" name="collision energy" value="27" unitCvRef="UO" unitAccession="UO:0000266" unitName="electronvolt"/>\n',
            "              </activation>\n",
            "            </precursor>\n",
            "          </precursorList>\n",
        ]
    s += [
        '          <binaryDataArrayList count="2">\n',
        array_xml(mz, "MS:1000514", "m/z array", 'unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"'),
        array_xml(inten, "MS:1000515", "intensity array", 'unitCvRef="MS" unitAccession="MS:1000131" unitName="number of detector counts"'),
        "          </binaryDataArrayList>\n",
        "        </spectrum>\n",
    ]
    return "".join(s)


cycles = int(GRADIENT_MIN * 60.0 / CYCLE_S)
scan_dt = CYCLE_S / (1 + len(WINDOWS))
count = cycles * (1 + len(WINDOWS))
head = f"""<?xml version="1.0" encoding="utf-8"?>
<indexedmzML xmlns="http://psi.hupo.org/ms/mzml" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.2_idx.xsd">
  <mzML xmlns="http://psi.hupo.org/ms/mzml" id="synthetic_dia" version="1.1.0">
    <cvList count="2">
      <cv id="MS" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology" version="4.1.41" URI="https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo"/>
      <cv id="UO" fullName="Unit Ontology" version="09:04:2014" URI="https://raw.githubusercontent.com/bio-ontology-research-group/unit-ontology/master/unit.obo"/>
    </cvList>
    <fileDescription>
      <fileContent>
        <cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum" value=""/>
        <cvParam cvRef="MS" accession="MS:1000580" name="MSn spectrum" value=""/>
        <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" value=""/>
      </fileContent>
    </fileDescription>
    <softwareList count="1">
      <software id="synthetic" version="0">
        <cvParam cvRef="MS" accession="MS:1000799" name="custom unreleased software tool" value="make_synthetic_dia"/>
      </software>
    </softwareList>
    <instrumentConfigurationList count="1">
      <instrumentConfiguration id="IC">
        <cvParam cvRef="MS" accession="MS:1000031" name="instrument model" value=""/>
      </instrumentConfiguration>
    </instrumentConfigurationList>
    <dataProcessingList count="1">
      <dataProcessing id="synthetic">
        <processingMethod order="0" softwareRef="synthetic">
          <cvParam cvRef="MS" accession="MS:1000544" name="Conversion to mzML" value=""/>
        </processingMethod>
      </dataProcessing>
    </dataProcessingList>
    <run id="synthetic_dia" defaultInstrumentConfigurationRef="IC">
      <spectrumList count="{count}" defaultDataProcessingRef="synthetic">
"""
with open(args.mzml, "wb") as out:
    out.write(head.encode())
    offsets = []
    index = 0
    for c in range(cycles):
        for w in [None] + list(range(len(WINDOWS))):
            t = c * CYCLE_S + (0 if w is None else (w + 1) * scan_dt)
            mz, inten = ms1(t) if w is None else ms2(t, w)
            offsets.append(out.tell() + 8)  # the '<' of '<spectrum', after its indent
            out.write(spectrum_xml(index, t, mz, inten, w).encode())
            index += 1
    out.write(b"      </spectrumList>\n    </run>\n  </mzML>\n")
    index_list = out.tell() + 2
    out.write(b'  <indexList count="2">\n    <index name="spectrum">\n')
    for i, o in enumerate(offsets):
        out.write(f'      <offset idRef="scan={i + 1}">{o}</offset>\n'.encode())
    out.write(b'    </index>\n    <index name="chromatogram">\n    </index>\n  </indexList>\n')
    out.write(f"  <indexListOffset>{index_list}</indexListOffset>\n</indexedmzML>\n".encode())
print(f"{count} spectra, {args.spiked} spiked precursors, RT range {prt.min():.0f}-{prt.max():.0f} s")
