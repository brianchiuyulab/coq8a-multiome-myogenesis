"""Display MYOD1 region effects for pre- and post-linkage region sets.

The 32/19-peak sets are COQ-outcome independent. The linked sets were learned
from same-nucleus peak/RNA associations in this dataset and remain exploratory.
"""

import argparse
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from multiome_core import h5_for, paired_stats, read_barcodes


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--h5-root", type=Path, required=True)
    p.add_argument("--tables", type=Path, required=True)
    a = p.parse_args()
    candidates = pd.read_csv(a.tables / "candidate_peak_gene.tsv.gz", sep="\t")
    myod = candidates[candidates.gene == "MYOD1"]
    links = pd.read_csv(a.tables / "peak_gene_links.tsv", sep="\t")
    links = links[
        (links.gene == "MYOD1")
        & (links.q_all_links < 0.05)
        & (links.partial_r > 0)
        & (links.n_libraries == 4)
    ]
    links_tss2 = pd.read_csv(a.tables / "peak_gene_links_tss2.tsv", sep="\t")
    links_tss2 = links_tss2[
        (links_tss2.gene == "MYOD1")
        & (links_tss2.q_all_links < 0.05)
        & (links_tss2.partial_r > 0)
        & (links_tss2.n_libraries == 4)
    ]
    region_sets = {
        "all_candidate_32": set(myod.peak),
        "common_all4_19": set(myod.loc[myod.n_libraries == 4, "peak"]),
        "positive_link_10": set(links.peak),
        "positive_link_nonMRF_5": set(links.loc[links.mrf_max_score < 0.95, "peak"]),
        "positive_link_nonMRF_tss2_6": set(links_tss2.loc[links_tss2.mrf_max_score < 0.95, "peak"]),
    }
    assert list(map(len, region_sets.values())) == [32, 19, 10, 5, 6]
    mapping = pd.read_csv(a.tables / "consensus_peak_map.tsv.gz", sep="\t").set_index("peak")
    pairs = pd.read_csv(a.tables / "matched_pairs.tsv.gz", sep="\t")
    rows, raw = [], []
    for gsm, gsm_pairs in pairs.groupby("gsm"):
        barcodes = sorted(set(gsm_pairs.high_barcode) | set(gsm_pairs.low_barcode))
        with h5py.File(h5_for(a.h5_root, gsm)) as h5:
            _, P, _, peak_names = read_barcodes(h5, barcodes)
        bidx = {b: i for i, b in enumerate(barcodes)}
        pidx = {str(p): i for i, p in enumerate(peak_names)}
        local_sets = {}
        for name, anchors in region_sets.items():
            local = []
            for anchor in anchors:
                target = anchor if gsm == "GSM6339597" else mapping.at[anchor, gsm + "_peak"]
                if isinstance(target, str) and target in pidx:
                    local.append(pidx[target])
            local_sets[name] = sorted(set(local))
        for (gate, contrast), ps in gsm_pairs.groupby(["gate", "contrast"]):
            hi = [bidx[b] for b in ps.high_barcode]
            lo = [bidx[b] for b in ps.low_barcode]
            for region_set, cols in local_sets.items():
                if not cols:
                    continue
                high = np.asarray(P[hi][:, cols].mean(axis=1)).ravel()
                low = np.asarray(P[lo][:, cols].mean(axis=1)).ravel()
                d = high - low
                stats = paired_stats(d)
                rows.append(
                    {
                        "gsm": gsm,
                        "source": ps.source.iloc[0],
                        "stage": ps.stage.iloc[0],
                        "gate": gate,
                        "contrast": contrast,
                        "region_set": region_set,
                        "n_regions": len(cols),
                        "high_open_pct": 100 * high.mean(),
                        "low_open_pct": 100 * low.mean(),
                        "fold_open": high.mean() / low.mean() if low.mean() else np.nan,
                        "delta_pp": 100 * d.mean(),
                        **stats,
                    }
                )
                raw.extend(
                    {
                        "gsm": gsm,
                        "gate": gate,
                        "contrast": contrast,
                        "region_set": region_set,
                        "pair": int(pid),
                        "high_barcode": hb,
                        "low_barcode": lb,
                        "high_open": h,
                        "low_open": l,
                        "difference": delta,
                    }
                    for pid, hb, lb, h, l, delta in zip(
                        ps.pair, ps.high_barcode, ps.low_barcode, high, low, d
                    )
                )
        print(gsm, "locus complete", flush=True)
    sample = pd.DataFrame(rows)
    sample.to_csv(a.tables / "myod1_locus_by_library.tsv", sep="\t", index=False)
    detail = pd.DataFrame(raw)
    detail.to_csv(
        a.tables / "myod1_locus_pair_scores.tsv.gz",
        sep="\t",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    summary = []
    for (gate, contrast, region_set), sub in detail.groupby(["gate", "contrast", "region_set"]):
        ss = sample[
            (sample.gate == gate)
            & (sample.contrast == contrast)
            & (sample.region_set == region_set)
        ]
        high, low = sub.high_open.mean(), sub.low_open.mean()
        summary.append(
            {
                "gate": gate,
                "contrast": contrast,
                "region_set": region_set,
                "n_regions": int(ss.n_regions.max()),
                "high_open_pct": 100 * high,
                "low_open_pct": 100 * low,
                "fold_open": high / low if low else np.nan,
                "delta_pp": 100 * (high - low),
                "positive_libraries": int((ss.delta_pp > 0).sum()),
                **paired_stats(sub.difference),
            }
        )
    summary = pd.DataFrame(summary)
    summary.to_csv(a.tables / "myod1_locus_summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
