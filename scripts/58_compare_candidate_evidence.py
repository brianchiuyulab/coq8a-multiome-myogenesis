"""Combine existing human associations with an independent differentiation reference."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def main(a):
    a.out.mkdir(parents=True, exist_ok=True)
    evidence = pd.read_csv(
        a.root / "results/unstratified/complete_followup_evidence.tsv", sep="\t"
    )
    annotation = pd.read_csv(
        a.root / "results/external_enhancers/external_peak_annotations.tsv", sep="\t"
    )
    selectivity = pd.read_csv(
        a.root / "results/enhancer_regions/peak_selectivity.tsv", sep="\t"
    )
    mouse = pd.read_csv(
        a.root / "results/external_c2c12/human_to_mouse_ATAC_annotation.tsv", sep="\t"
    )
    time = pd.read_csv(
        a.root / "results/external_c2c12/all_gene_time_contrasts.tsv.gz", sep="\t"
    )
    time["symbol_match"] = time.gene.str.upper()
    unique = time.groupby("symbol_match").ID.nunique()
    time = time[time.symbol_match.isin(unique[unique.eq(1)].index)]
    lib = pd.read_csv(
        a.root / "results/tables/candidate_peak_effects_by_library.tsv.gz", sep="\t"
    )
    lib = lib[(lib.gate == "TSS_ge_3") & (lib.contrast == "3plus_vs_1")]
    inventory = pd.read_csv(
        a.root / "results/design/target_sample_inventory.tsv", sep="\t"
    )
    lib = lib.merge(
        inventory[["gsm", "source", "stage"]], on="gsm", validate="many_to_one"
    )
    lib["high_fraction"] = lib.high_open / lib.n_pairs
    lib["low_fraction"] = lib.low_open / lib.n_pairs
    source = (
        lib.groupby(["peak", "source"])[["high_fraction", "low_fraction"]]
        .mean()
        .reset_index()
    )
    source["FC_equal_stage"] = source.high_fraction / source.low_fraction.replace(
        0, np.nan
    )
    source["delta_pp"] = 100 * (source.high_fraction - source.low_fraction)
    source.to_csv(
        a.out / "peak_effects_by_source.tsv.gz",
        sep="\t",
        index=False,
        na_rep="NA",
        compression={"method": "gzip", "mtime": 0},
    )
    wide = source.pivot(
        index="peak", columns="source", values="FC_equal_stage"
    ).add_prefix("ATAC_FC_")
    evidence = evidence.merge(wide, on="peak", how="left")
    evidence = evidence.merge(
        annotation[["peak", "strong_enhancer", "MYOD_bound_enhancer"]],
        on="peak",
        how="left",
    )
    evidence = evidence.merge(
        selectivity[["peak", "n_other_strong"]], on="peak", how="left"
    )
    evidence = evidence.merge(mouse, on="peak", how="left")
    for t in [12, 24, 48, 60, 96]:
        sub = time[time.time_h.eq(t)][["symbol_match", "FC", "pvalue", "padj"]].rename(
            columns={
                "symbol_match": "gene",
                "FC": f"mouse_RNA_FC_{t}h",
                "pvalue": f"mouse_RNA_p_{t}h",
                "padj": f"mouse_RNA_q_{t}h",
            }
        )
        evidence = evidence.merge(sub, on="gene", how="left", validate="many_to_one")
    evidence.to_csv(
        a.out / "complete_existing_followup_with_external_evidence.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )
    positive = evidence[
        (evidence.fold_open > 1)
        & (evidence.p_pair_binomial < 0.05)
        & (evidence.RNA_ratio > 1)
        & (evidence.p_peak_RNA < 0.05)
        & (evidence.difference > 0)
        & (evidence.p_pair < 0.05)
    ].copy()
    positive.to_csv(
        a.out / "nominal_positive_chain_with_external_evidence.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )
    focal = {
        "MYOD1": "chr11:17649919-17650798",
        "CAV3": "chr3:8733438-8733968",
        "DMD": "chrX:31220862-31221781",
    }
    selected = pd.concat(
        [
            evidence[(evidence.gene == g) & (evidence.peak == p)]
            for g, p in focal.items()
        ]
    )
    cols = [
        "gene",
        "peak",
        "fold_open",
        "p_pair_binomial",
        "ATAC_FC_line1",
        "ATAC_FC_line2",
        "RNA_ratio",
        "p_peak_RNA",
        "difference",
        "p_pair",
        "strong_enhancer",
        "n_other_strong",
        "mapping_pass",
        "GM_overlap",
        "DM60h_overlap",
        "mouse_RNA_FC_24h",
        "mouse_RNA_q_24h",
        "mouse_RNA_FC_48h",
        "mouse_RNA_q_48h",
    ]
    selected[cols].to_csv(
        a.out / "focal_comparison.tsv", sep="\t", index=False, na_rep="NA"
    )
    print(selected[cols].to_string(index=False))
    corr = pd.read_csv(
        a.root / "results/unstratified/correlations_combined.tsv.gz", sep="\t"
    )
    closing = corr[
        corr.peak.isin(["chr11:1880073-1880969", "chr11:1882161-1883026"])
        & corr.model.eq("technical_state_coq")
    ]
    closing.to_csv(a.out / "closing_block_RNA_links.tsv", sep="\t", index=False)
    checks = []
    for file in sorted(a.cache.glob("GSM*_meta.tsv")):
        d = pd.read_csv(file, sep="\t")
        checks.append(
            dict(
                gsm=file.name.split("_")[0],
                n_QC=len(d),
                n_COQ8A_ge3=int((d.COQ8A_umi >= 3).sum()),
                n_COQ8A_eq1=int((d.COQ8A_umi == 1).sum()),
            )
        )
    pd.DataFrame(checks).merge(
        inventory[["gsm", "n_main_matched_pairs"]], on="gsm"
    ).to_csv(a.out / "COQ8A_group_inventory.tsv", sep="\t", index=False)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument("--cache", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    main(p.parse_args())
