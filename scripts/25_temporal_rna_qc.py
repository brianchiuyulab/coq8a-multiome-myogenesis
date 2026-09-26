"""Audit RNA features/mitochondrial counts and challenge fixed-pair temporal effects."""

import argparse
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from multiome_core import h5_for, read_barcodes, paired_stats


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ["h5-root", "tables", "temporal", "scrublet"]:
        p.add_argument("--" + arg, type=Path, required=True)
    a = p.parse_args()
    pairs = pd.read_csv(a.tables / "matched_pairs.tsv.gz", sep="\t")
    records = []
    for gsm, ps in pairs.groupby("gsm"):
        barcodes = sorted(set(ps.high_barcode) | set(ps.low_barcode))
        with h5py.File(h5_for(a.h5_root, gsm)) as h5:
            rna, _, genes, _ = read_barcodes(h5, barcodes)
        total = np.asarray(rna.sum(axis=1)).ravel()
        mito = np.char.startswith(genes, "MT-")
        records.append(
            pd.DataFrame(
                {
                    "gsm": gsm,
                    "barcode": barcodes,
                    "n_rna_genes": rna.getnnz(axis=1),
                    "pct_mito_umi": 100
                    * np.asarray(rna[:, mito].sum(axis=1)).ravel()
                    / total,
                }
            )
        )
    cells = pd.concat(records, ignore_index=True)
    cells.to_csv(a.temporal / "matched_nucleus_rna_qc.tsv.gz", sep="\t", index=False)
    scrublet = pd.read_csv(a.scrublet, sep="\t")
    cells = cells.merge(scrublet, on=["gsm", "barcode"], validate="one_to_one")
    assert cells.scrublet_score.notna().all()
    rules = {
        "Baseline": np.ones(len(cells), dtype=bool),
        "RNA_genes_ge500_mito_le5pct": (cells.n_rna_genes >= 500)
        & (cells.pct_mito_umi <= 5),
        "RNA_genes_ge500_mito_le10pct": (cells.n_rna_genes >= 500)
        & (cells.pct_mito_umi <= 10),
        "RNA_genes_ge500_mito_le15pct": (cells.n_rna_genes >= 500)
        & (cells.pct_mito_umi <= 15),
        "RNA_genes_ge500_mito_le20pct": (cells.n_rna_genes >= 500)
        & (cells.pct_mito_umi <= 20),
        "Remove_predicted_doublets": ~cells.scrublet_predicted_doublet.astype(bool),
    }
    for tail in [1, 2.5, 5, 10]:
        cut = scrublet.groupby("gsm").scrublet_score.quantile(1 - tail / 100)
        rules[f"Remove_scrublet_top{tail}pct"] = cells.scrublet_score <= cells.gsm.map(
            cut
        )
    scores = pd.read_csv(a.temporal / "matched_pair_temporal_scores.tsv.gz", sep="\t")
    scores = scores[scores.set_id.str.startswith("P2_O50_")]
    joined = scores.merge(
        pairs[["gsm", "gate", "contrast", "pair", "high_barcode", "low_barcode"]],
        on=["gsm", "gate", "contrast", "pair"],
        validate="many_to_one",
    )
    rows = []
    for name, mask in rules.items():
        passing = set(zip(cells.loc[mask, "gsm"], cells.loc[mask, "barcode"]))
        keep = np.array(
            [
                (g, h) in passing and (g, l) in passing
                for g, h, l in zip(joined.gsm, joined.high_barcode, joined.low_barcode)
            ]
        )
        for stage in ["pooled", "stem", "differentiated"]:
            sub = joined[keep]
            if stage != "pooled":
                sub = sub[sub.stage == stage]
            for key, frame in sub.groupby(["set_id", "gate", "contrast"]):
                high, low = frame.high_open.mean(), frame.low_open.mean()
                s = paired_stats(frame.high_open - frame.low_open)
                directions = (
                    frame.assign(delta=frame.high_open - frame.low_open)
                    .groupby("gsm")
                    .delta.mean()
                    > 1e-10
                ).sum()
                rows.append(
                    dict(
                        qc_rule=name,
                        stage_analysis=stage,
                        set_id=key[0],
                        gate=key[1],
                        contrast=key[2],
                        n_pairs=len(frame),
                        fold_open=high / low if low else np.nan,
                        delta_pp=100 * (high - low),
                        p_pair=s["p_pair"],
                        positive_libraries=int(directions),
                    )
                )
    result = pd.DataFrame(rows)
    all_settings = joined[["set_id", "gate", "contrast"]].drop_duplicates()
    full = (
        pd.DataFrame({"qc_rule": list(rules)})
        .merge(
            pd.DataFrame({"stage_analysis": ["pooled", "stem", "differentiated"]}),
            how="cross",
        )
        .merge(all_settings, how="cross")
    )
    result = full.merge(
        result,
        on=["qc_rule", "stage_analysis", "set_id", "gate", "contrast"],
        how="left",
        validate="one_to_one",
    )
    result["n_pairs"] = result.n_pairs.fillna(0).astype(int)
    result.to_csv(
        a.temporal / "rna_qc_sensitivity.tsv", sep="\t", index=False, na_rep="NA"
    )


if __name__ == "__main__":
    main()
