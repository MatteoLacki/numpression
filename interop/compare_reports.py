"""Compare DIA-NN reports of the same synthetic run in different encodings.

Usage: compare_reports.py <spiked.tsv> <reference report.parquet> <report.parquet>...
Per report: IDs at 1% run-level q-value, how many were spiked, overlap with the
reference, and precursor quantity and RT differences against it.
"""

import sys

import pandas as pd

spiked = set(pd.read_csv(sys.argv[1], sep="\t")["Precursor.Id"])
paths = sys.argv[2:]
reports = {p: pd.read_parquet(p) for p in paths}
reference = paths[0]
ids = {p: set(r.loc[r["Q.Value"] <= 0.01, "Precursor.Id"]) for p, r in reports.items()}

print(f"{len(spiked)} spiked precursors; reference {reference}")
print(f"{'IDs@1%':>7} {'spiked':>7} {'other':>6} {'shared':>7} {'+':>4} {'-':>4} {'|dQ| med':>9} {'|dQ| p99':>9} {'max |dRT| s':>11}  report")
base = reports[reference].set_index("Precursor.Id")
for p, r in reports.items():
    m = base[["Precursor.Quantity", "RT"]].join(
        r.set_index("Precursor.Id")[["Precursor.Quantity", "RT"]], rsuffix="_n", how="inner"
    )
    m = m[(m["Precursor.Quantity"] > 0) & (m["Precursor.Quantity_n"] > 0)]
    dq = (m["Precursor.Quantity_n"] / m["Precursor.Quantity"] - 1).abs()
    drt = (m["RT_n"] - m["RT"]).abs().max() * 60
    s = ids[p]
    print(
        f"{len(s):7d} {len(s & spiked):7d} {len(s - spiked):6d} {len(s & ids[reference]):7d}"
        f" {len(s - ids[reference]):4d} {len(ids[reference] - s):4d}"
        f" {dq.median():9.1e} {dq.quantile(0.99):9.1e} {drt:11.2f}  {p}"
    )
