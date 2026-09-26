"""Test every eligible candidate peak before selecting a locus for follow-up."""

import argparse
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from multiome_core import h5_for, read_barcodes


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--h5-root", type=Path, required=True)
    p.add_argument("--tables", type=Path, required=True)
    a = p.parse_args()
    pairs = pd.read_csv(a.tables / "matched_pairs.tsv.gz", sep="\t")
    candidate = pd.read_csv(a.tables / "candidate_peak_gene.tsv.gz", sep="\t")
    mapping = pd.read_csv(a.tables / "consensus_peak_map.tsv.gz", sep="\t").set_index("peak")
    anchors = sorted(candidate.peak.unique())
    by_library = []
    for gsm, sample_pairs in pairs.groupby("gsm"):
        barcodes = sorted(set(sample_pairs.high_barcode) | set(sample_pairs.low_barcode))
        with h5py.File(h5_for(a.h5_root, gsm)) as h5:
            _, P, _, peak_names = read_barcodes(h5, barcodes)
        bidx = {b: i for i, b in enumerate(barcodes)}
        pidx = {str(p): i for i, p in enumerate(peak_names)}
        local = [
            pidx.get(anchor if gsm == "GSM6339597" else mapping.at[anchor, gsm + "_peak"], -1)
            for anchor in anchors
        ]
        valid = np.asarray(local) >= 0
        available = np.asarray(anchors)[valid]
        X = P[:, np.asarray(local)[valid]].tocsr()
        for (gate, contrast), sub in sample_pairs.groupby(["gate", "contrast"]):
            hi = np.asarray([bidx[b] for b in sub.high_barcode])
            lo = np.asarray([bidx[b] for b in sub.low_barcode])
            high, low = X[hi], X[lo]
            h = np.asarray(high.sum(axis=0)).ravel().astype(int)
            l = np.asarray(low.sum(axis=0)).ravel().astype(int)
            both = np.asarray(high.multiply(low).sum(axis=0)).ravel().astype(int)
            plus, minus = h - both, l - both
            disc = plus + minus
            row = pd.DataFrame(
                {
                    "gsm": gsm,
                    "gate": gate,
                    "contrast": contrast,
                    "peak": available,
                    "n_pairs": len(sub),
                    "high_open": h,
                    "low_open": l,
                    "high_only": plus,
                    "low_only": minus,
                    "delta_pp": 100 * (h - l) / len(sub),
                }
            )
            by_library.append(row)
        print(gsm, "peaks", len(available), flush=True)
    lib = pd.concat(by_library, ignore_index=True)
    lib.to_csv(
        a.tables / "candidate_peak_effects_by_library.tsv.gz",
        sep="\t",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    rows = []
    for (gate, contrast, peak), sub in lib.groupby(["gate", "contrast", "peak"]):
        if len(sub) < 3:
            continue
        h, l = int(sub.high_open.sum()), int(sub.low_open.sum())
        plus, minus = int(sub.high_only.sum()), int(sub.low_only.sum())
        disc = plus + minus
        pval = min(1.0, 2 * stats.binom.cdf(min(plus, minus), disc, 0.5)) if disc else 1.0
        rows.append(
            {
                "gate": gate,
                "contrast": contrast,
                "peak": peak,
                "n_libraries": len(sub),
                "n_pairs": int(sub.n_pairs.sum()),
                "high_open": h,
                "low_open": l,
                "high_only": plus,
                "low_only": minus,
                "delta_pp": 100 * (h - l) / sub.n_pairs.sum(),
                "positive_libraries": int((sub.delta_pp > 0).sum()),
                "p_pair_binomial": pval,
            }
        )
    out = pd.DataFrame(rows)
    out["q_candidate_peaks"] = np.nan
    for _, sub in out.groupby(["gate", "contrast"]):
        test = (sub.high_open + sub.low_open) >= 20
        if test.any():
            out.loc[sub.index[test], "q_candidate_peaks"] = multipletests(
                sub.loc[test, "p_pair_binomial"], method="fdr_bh"
            )[1]
    out.to_csv(
        a.tables / "candidate_peak_effects_pooled.tsv.gz",
        sep="\t",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    primary = out[(out.gate == "TSS_ge_3") & (out.contrast == "2plus_vs_1")]
    print("primary ATAC peak q<.05", int((primary.q_candidate_peaks < 0.05).sum()))


if __name__ == "__main__":
    main()
