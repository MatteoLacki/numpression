#!/usr/bin/env bash
# DIA-NN interop check: predict a library from a FASTA, simulate a DIA run from it,
# re-encode the run with numpression, search every encoding and compare the reports.
#
# Usage: run_diann.sh <diann-linux> <fasta> <workdir> [python]
# python needs numpy, pandas and pyarrow. A FASTA of ~200 proteins takes a few minutes.
set -euo pipefail

diann=$(realpath "$1")
fasta=$(realpath "$2")
mkdir -p "$3"
work=$(realpath "$3")
python=${4:-python3}
here=$(cd "$(dirname "$0")" && pwd)

cargo build --release --manifest-path "$here/renumpress/Cargo.toml"
renumpress="$here/renumpress/target/release/renumpress"

cd "$work"
"$diann" --fasta "$fasta" --fasta-search --predictor --gen-spec-lib \
    --out-lib lib.parquet --threads 16 > predict.log 2>&1
"$diann" --lib lib.predicted.speclib --gen-spec-lib --out-lib lib.parquet \
    --threads 16 > export.log 2>&1

"$python" "$here/make_synthetic_dia.py" lib.parquet synthetic.mzML spiked.tsv
for mode in zlib32 pic slof; do
    "$renumpress" synthetic.mzML "synthetic_$mode.mzML" "$mode"
done

for name in synthetic_zlib32 synthetic_pic synthetic_slof; do
    mkdir -p "run_$name"
    "$diann" --f "$work/$name.mzML" --lib lib.predicted.speclib --threads 16 \
        --temp "run_$name" --out "run_$name/report.parquet" > "run_$name/log.txt" 2>&1
done

"$python" "$here/compare_reports.py" spiked.tsv \
    run_synthetic_zlib32/report.parquet run_synthetic_pic/report.parquet \
    run_synthetic_slof/report.parquet
