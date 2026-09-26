"""Retain fixed matched pairs after removing per-library Scrublet score tails."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from multiome_core import paired_stats


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tables", type=Path, required=True)
    p.add_argument("--scrublet", type=Path, required=True)
    a = p.parse_args()
    d = pd.read_csv(a.tables / "myod1_locus_pair_scores.tsv.gz", sep="\t")
    d = d[d.gate == "TSS_ge_3"]
    s = pd.read_csv(a.scrublet, sep="\t")
    score = s.set_index(["gsm", "barcode"]).scrublet_score
    rows = []
    for tail in [0, 1, 2.5, 5, 10]:
        if tail == 0:
            keep = d
        else:
            cutoff = s.groupby("gsm").scrublet_score.quantile(1 - tail / 100)
            hi = np.array([score.at[(g, b)] <= cutoff.at[g] for g, b in zip(d.gsm, d.high_barcode)])
            lo = np.array([score.at[(g, b)] <= cutoff.at[g] for g, b in zip(d.gsm, d.low_barcode)])
            keep = d[hi & lo]
        for (contrast, region_set), sub in keep.groupby(["contrast", "region_set"]):
            h, l = sub.high_open.mean(), sub.low_open.mean()
            sample = sub.groupby("gsm").difference.mean()
            rows.append(
                {
                    "contrast": contrast,
                    "region_set": region_set,
                    "scrublet_tail_removed_pct": tail,
                    "n_pairs": len(sub),
                    "high_open_pct": 100 * h,
                    "low_open_pct": 100 * l,
                    "fold_open": h / l if l else np.nan,
                    "delta_pp": 100 * (h - l),
                    "positive_libraries": int((sample > 0).sum()),
                    "p_pair": paired_stats(sub.difference)["p_pair"],
                }
            )
    out = pd.DataFrame(rows)
    out.to_csv(a.tables / "doublet_tail_sensitivity.tsv", sep="\t", index=False)
    print(
        out[
            out.region_set.isin(
                ["common_all4_19", "positive_link_nonMRF_5", "positive_link_nonMRF_tss2_6"]
            )
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
