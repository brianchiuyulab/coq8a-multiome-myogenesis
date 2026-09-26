"""Audit peak-level BH adjustment under explicit candidate scopes."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


def main(root, gene):
    tables = root / "results/tables"
    effects = pd.read_csv(tables / "candidate_peak_effects_pooled.tsv.gz", sep="\t")
    effects = effects.loc[
        effects.gate.eq("TSS_ge_3")
        & effects.contrast.eq("3plus_vs_1")
        & effects.n_libraries.eq(4)
    ].copy()
    if effects.peak.duplicated().any():
        raise ValueError("Expected one row per peak for the fixed contrast")
    effects["FC"] = effects.high_open / effects.low_open.replace(0, np.nan)
    annotation = pd.read_csv(
        root / "results/enhancer_regions/peak_selectivity.tsv", sep="\t"
    )
    strict = annotation.loc[
        annotation.strong_enhancer & annotation.n_other_strong.eq(0), "peak"
    ]
    enhancers = annotation.loc[annotation.strong_enhancer, "peak"]
    neighbors = pd.read_csv(tables / "candidate_peak_gene.tsv.gz", sep="\t")
    gene_peaks = neighbors.loc[neighbors.gene.eq(gene), "peak"]
    scopes = {
        "myogenesis_candidates": effects.peak,
        "HSMM_strong_enhancers": enhancers,
        "muscle_selective_enhancers": strict,
        "gene_selective_scope_posthoc": strict[strict.isin(gene_peaks)],
    }
    output = []
    for name, peaks in scopes.items():
        frame = effects.loc[effects.peak.isin(peaks)].copy()
        if frame.empty or not np.isfinite(frame.p_pair_binomial).all():
            raise ValueError(f"Empty or nonfinite hypothesis family: {name}")
        frame["q_BH"] = multipletests(frame.p_pair_binomial, method="fdr_bh")[1]
        frame["scope"] = name
        frame["n_hypotheses"] = len(frame)
        frame = frame.loc[frame.peak.isin(gene_peaks)]
        output.append(frame)
    result = pd.concat(output, ignore_index=True)
    result = result.merge(
        neighbors.loc[neighbors.gene.eq(gene), ["peak", "gene", "nearest_tss_bp"]],
        on="peak",
        validate="many_to_one",
    )
    out = root / "results/candidate_comparison"
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{gene}_fdr_scope_audit.tsv"
    result.to_csv(path, sep="\t", index=False, na_rep="NA")
    print(path)
    print(result.groupby("scope").n_hypotheses.first().to_string())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument("--gene", default="CAV3")
    args = parser.parse_args()
    main(args.root, args.gene)
