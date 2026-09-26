"""Validate temporal membership, paired summaries and direct H5 spot checks."""

import argparse
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import stats

from multiome_core import h5_for


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ["tables", "temporal"]:
        p.add_argument("--" + arg, type=Path, required=True)
    p.add_argument(
        "--h5-root",
        type=Path,
        help="Include direct raw-H5 spot checks when inputs are available",
    )
    a = p.parse_args()
    members = pd.read_csv(a.temporal / "temporal_set_memberships.tsv.gz", sep="\t")
    definitions = pd.read_csv(a.temporal / "temporal_set_definitions.tsv", sep="\t")
    annotations = pd.read_csv(
        a.temporal / "target_temporal_annotations.tsv.gz", sep="\t"
    )
    candidates = pd.read_csv(a.tables / "candidate_peak_gene.tsv.gz", sep="\t")
    assert set(members.peak) <= set(candidates.loc[candidates.n_libraries == 4, "peak"])
    assert not members.duplicated(["set_id", "peak"]).any()
    for row in definitions.itertuples():
        assert row.n_peaks == int((members.set_id == row.set_id).sum())
    mid = set(members.loc[members.set_id == "P2_O50_Middle_24to48h", "peak"])
    annotated = annotations[(annotations.promoter_kb == 2) & annotations.peak.isin(mid)]
    assert (annotated.q_dynamic < 0.05).all() and (annotated.rep_delta72 > 0).all()
    assert ((annotated.half_rise_hour > 24) & (annotated.half_rise_hour <= 48)).all()
    scores = pd.read_csv(a.temporal / "matched_pair_temporal_scores.tsv.gz", sep="\t")
    effects = pd.read_csv(a.temporal / "temporal_effect_grid.tsv", sep="\t")
    checked = 0
    for r in effects.query(
        "scope=='programme' and stage_analysis=='pooled'"
    ).itertuples():
        x = scores[
            (scores.set_id == r.set_id)
            & (scores.gate == r.gate)
            & (scores.contrast == r.contrast)
        ]
        h, l = x.high_open.to_numpy(), x.low_open.to_numpy()
        difference = h - l
        assert len(x) == r.n_pairs
        assert np.isclose(h.mean() / l.mean(), r.fold_open, rtol=1e-11)
        t = difference.mean() / (difference.std(ddof=1) / np.sqrt(len(x)))
        pvalue = 2 * stats.t.sf(abs(t), len(x) - 1)
        assert np.isclose(pvalue, r.p_pair, atol=1e-12)
        checked += 1
    if a.h5_root is None:
        print(
            json.dumps(
                dict(
                    status="passed",
                    programme_summaries_recomputed=checked,
                    mode="released_tables_only",
                    external_middle_peaks=len(mid),
                )
            )
        )
        return
    pairs = pd.read_csv(a.tables / "matched_pairs.tsv.gz", sep="\t")
    pairs = pairs[(pairs.gate == "TSS_ge_3") & (pairs.contrast == "3plus_vs_1")]
    mapping = pd.read_csv(a.tables / "consensus_peak_map.tsv.gz", sep="\t").set_index(
        "peak"
    )
    direct_checks = 0
    for gsm, lib in pairs.groupby("gsm"):
        with h5py.File(h5_for(a.h5_root, gsm)) as h5:
            matrix = h5["matrix"]
            names = np.char.decode(matrix["features"]["name"][:])
            barcodes = np.char.decode(matrix["barcodes"][:])
            fidx = {str(x): i for i, x in enumerate(names)}
            bidx = {str(x): i for i, x in enumerate(barcodes)}
            target = (
                list(mid)
                if gsm == "GSM6339597"
                else mapping.loc[sorted(mid), gsm + "_peak"].tolist()
            )
            target_ids = {fidx[x] for x in target}
            for pair in [lib.iloc[0], lib.iloc[-1]]:
                stored = scores[
                    (scores.set_id == "P2_O50_Middle_24to48h")
                    & (scores.gsm == gsm)
                    & (scores.gate == "TSS_ge_3")
                    & (scores.contrast == "3plus_vs_1")
                    & (scores.pair == pair.pair)
                ].iloc[0]
                for group in ["high", "low"]:
                    barcode = pair[group + "_barcode"]
                    col = bidx[barcode]
                    left, right = matrix["indptr"][col : col + 2]
                    ids = matrix["indices"][left:right]
                    values = matrix["data"][left:right]
                    nonzero = set(ids[values > 0].tolist())
                    opened = len(nonzero & target_ids) / len(target_ids)
                    assert np.isclose(opened, stored[group + "_open"], atol=1e-12)
                    direct_checks += 1
    report = dict(
        status="passed",
        programme_summaries_recomputed=checked,
        direct_H5_nucleus_scores_checked=direct_checks,
        external_middle_peaks=len(mid),
        main_matched_pairs=len(pairs),
        libraries=pairs.gsm.nunique(),
        source_lines=pairs.source.nunique(),
    )
    (a.temporal / "validation.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report))


if __name__ == "__main__":
    main()
