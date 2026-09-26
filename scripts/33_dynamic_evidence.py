"""Assemble peak effects, cis links, and fixed gene-level RNA associations.

Every distance-based target remains in the table. No composite score or
significance-based target reassignment is applied. Missing RNA statistics
indicate a gene outside the original 221-gene RNA screen.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tables", type=Path, required=True)
    parser.add_argument("--temporal", type=Path, required=True)
    args = parser.parse_args()
    peaks = pd.read_csv(args.temporal / "bidirectional_peak_effects.tsv", sep="\t")
    links = pd.read_csv(args.temporal / "dynamic_cis_links.tsv", sep="\t")
    links = links.rename(columns={"positive_libraries": "link_positive_libraries"})
    original = pd.read_csv(args.tables / "candidate_peak_gene.tsv.gz", sep="\t")
    membership = set(zip(original.peak, original.gene))
    links["original_221_neighborhood"] = [
        (p, g) in membership for p, g in zip(links.peak, links.gene)
    ]
    evidence = links.merge(
        peaks, on=["gate", "contrast", "phase", "peak"], validate="many_to_one"
    )
    rna = pd.read_csv(args.tables / "gene_effects_pooled.tsv", sep="\t")
    rna = rna.loc[
        rna.modality.eq("RNA"),
        [
            "gate",
            "contrast",
            "gene",
            "difference",
            "p_pair",
            "q_221",
            "positive_libraries",
        ],
    ]
    rna = rna.rename(
        columns={c: "rna_" + c for c in rna if c not in ["gate", "contrast", "gene"]}
    )
    evidence = evidence.merge(
        rna, on=["gate", "contrast", "gene"], how="left", validate="many_to_one"
    )
    evidence["link_status"] = np.where(
        evidence.p_link.isna(),
        "insufficient_detection",
        np.where(
            evidence.q_all_testable_dynamic_cis_links.lt(0.05),
            "FDR_positive",
            np.where(evidence.p_link.lt(0.05), "nominal_only", "not_significant"),
        ),
    )
    evidence["atac_rna_same_direction"] = np.where(
        evidence.rna_difference.isna(),
        None,
        evidence.delta_pp * evidence.rna_difference > 0,
    )
    # Count directional persistence without choosing thresholds after seeing outcomes.
    sensitivity = peaks.groupby("peak").agg(
        sensitivity_positive_settings=("delta_pp", lambda x: int((x > 0).sum())),
        sensitivity_negative_settings=("delta_pp", lambda x: int((x < 0).sum())),
        sensitivity_nominal_settings=("p_exact", lambda x: int((x < 0.05).sum())),
    )
    evidence = evidence.merge(sensitivity, on="peak", validate="many_to_one")
    evidence.sort_values(
        ["gate", "contrast", "phase", "p_exact", "peak", "gene"]
    ).to_csv(
        args.temporal / "dynamic_evidence.tsv.gz",
        sep="\t",
        index=False,
        na_rep="NA",
        compression={"method": "gzip", "mtime": 0},
    )
    main = evidence.loc[
        evidence.gate.eq("TSS_ge_3") & evidence.contrast.eq("3plus_vs_1")
    ]
    near = main.loc[main.original_221_neighborhood].copy()
    near.to_csv(
        args.temporal / "dynamic_neighborhood_evidence_main.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )
    print("Full rows:", len(evidence), "; main original-neighborhood rows:", len(near))
    columns = [
        "phase",
        "peak",
        "gene",
        "fold_open",
        "p_exact",
        "q_within_phase_peaks",
        "partial_r",
        "p_link",
        "q_within_phase_links",
        "rna_difference",
        "rna_p_pair",
        "rna_q_221",
        "n_gene_detected",
        "n_eligible_libraries",
    ]
    print(
        near.loc[
            near.gene.isin(
                ["CSRP3", "GNAO1", "MYOD1", "MYOG", "LDB3", "TNNC1", "SPARC", "CDH13"]
            )
        ]
        .sort_values(["gene", "p_exact"])[columns]
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
